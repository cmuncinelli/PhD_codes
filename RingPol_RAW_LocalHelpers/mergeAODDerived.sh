#!/bin/bash

# ==============================================================================
# ALICE O2 DERIVED AOD MERGER
# ==============================================================================
# Companion script to runAnalysisInternalMultithreading.sh.
# Finds derived AO2D files produced by that script, builds the merger input
# list on the fly, and calls o2-aod-merger.
#
# Usage:
#   ./mergeAODDerived.sh <results_dir> [suffix] [max_folder_size_bytes]
# Example:
#   ./mergeAODDerived.sh ./results "-JustTOFBase-FullQA_permissivePt"
#
# Arguments:
#   results_dir          Path to the results/ directory written by the analysis
#                        launcher (contains the derived_aods/ sub-folder).
#   suffix               Optional. Config suffix used in the launcher (e.g.
#                        "-JustTOFBase-FullQA_permissivePt"). When provided,
#                        only files matching AO2D_Derived<suffix>_*.root are
#                        merged. When omitted, ALL AO2D_Derived_*.root files in
#                        the folder are merged.
#   max_folder_size      Optional. Maximum uncompressed DF-folder size in bytes
#                        passed to o2-aod-merger (default: 10000000000 = 10 GB).
#                        Set large to keep every DF folder separate (safest for
#                        index validation). Set small to stress-test remapping.
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# 1. ARGUMENT PARSING
# ------------------------------------------------------------------------------
RESULTS_DIR="${1:-}"
SUFFIX="${2:-}"          # Optional: e.g. "-JustTOFBase-FullQA_permissivePt"
MAX_SIZE="${3:-10000000000}"   # Default: 10 GB -- keeps DF folders separate

if [[ -z "$RESULTS_DIR" ]]; then
    echo "Usage: $0 <results_dir> [suffix] [max_folder_size_bytes]"
    exit 1
fi

RESULTS_DIR=$(realpath "$RESULTS_DIR")
DERIVED_AOD_DIR="${RESULTS_DIR}/derived_aods"
MERGE_OUTPUT_DIR="${RESULTS_DIR}/merged_aods"

if [[ ! -d "$DERIVED_AOD_DIR" ]]; then
    echo "ERROR: derived_aods directory not found under ${RESULTS_DIR}"
    exit 1
fi

mkdir -p "$MERGE_OUTPUT_DIR"

# ------------------------------------------------------------------------------
# 1.5. LOGGING SETUP
# ------------------------------------------------------------------------------
# The log file name includes the suffix when provided.
CLEAN_SUFFIX="${SUFFIX#-}"
if [[ -n "$CLEAN_SUFFIX" ]]; then
    MERGE_LOG="${MERGE_OUTPUT_DIR}/merge_${CLEAN_SUFFIX}.log"
else
    MERGE_LOG="${MERGE_OUTPUT_DIR}/merge.log"
fi

# ------------------------------------------------------------------------------
# 2. COLLECT INPUT FILES
# ------------------------------------------------------------------------------
# Build the glob pattern: if SUFFIX is given, match only that suffix;
# otherwise match everything produced by the launcher.
if [[ -n "$SUFFIX" ]]; then
    PATTERN="AO2D_Derived${SUFFIX}_*.root"
else
    PATTERN="AO2D_Derived_*.root"
fi

echo "========================================================"
echo "  ALICE O2 Derived AOD Merger"
echo "  Derived AOD directory : ${DERIVED_AOD_DIR}"
echo "  File pattern          : ${PATTERN}"
echo "  Max DF folder size    : ${MAX_SIZE} bytes"
echo "  Output directory      : ${MERGE_OUTPUT_DIR}"
echo "  Output log            : ${MERGE_LOG}"
echo "========================================================"

# Build the input list file (one absolute path per line, no 'file:' prefix --
# o2-aod-merger reads bare paths).
MERGER_INPUT="${MERGE_OUTPUT_DIR}/merger_input${SUFFIX}.txt"
> "$MERGER_INPUT"   # Truncate / create

# Use find + sort to guarantee a stable, numerically ordered list.
# The batch ID suffix (_00, _01, ...) is numeric, so a plain sort works.
while IFS= read -r -d '' F; do
    realpath "$F" >> "$MERGER_INPUT"
done < <(find "$DERIVED_AOD_DIR" -maxdepth 1 -name "$PATTERN" -print0 | sort -z)

NUM_FILES=$(wc -l < "$MERGER_INPUT")
if [[ "$NUM_FILES" -eq 0 ]]; then
    echo "ERROR: No files matching '${PATTERN}' found in ${DERIVED_AOD_DIR}"
    rm "$MERGER_INPUT"
    exit 1
fi

echo ""
echo "  Found ${NUM_FILES} file(s) to merge:"
cat "$MERGER_INPUT" | sed 's/^/    /'
echo ""

# ------------------------------------------------------------------------------
# 3. RUN THE MERGER
# ------------------------------------------------------------------------------
# Derive a meaningful output name from the suffix (strip leading dash if any).
CLEAN_SUFFIX="${SUFFIX#-}"   # "-JustTOFBase" -> "JustTOFBase", "" -> ""
if [[ -n "$CLEAN_SUFFIX" ]]; then
    OUTPUT_FILE="${MERGE_OUTPUT_DIR}/AO2D_DerivedMerged_${CLEAN_SUFFIX}.root"
else
    OUTPUT_FILE="${MERGE_OUTPUT_DIR}/AO2D_DerivedMerged.root"
fi

echo "  Running o2-aod-merger..."
echo "  Output: ${OUTPUT_FILE}"
echo ""

time o2-aod-merger \
    --input  "$MERGER_INPUT" \
    --output "$OUTPUT_FILE" \
    --max-size "$MAX_SIZE" \
    > "$MERGE_LOG" 2>&1

EXIT_CODE=$?

echo ""
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "========================================================"
    echo "  Merge complete."
    echo "  Output  : ${OUTPUT_FILE}"
    echo "  Input list kept at: ${MERGER_INPUT}"
    echo "========================================================"
else
    echo "ERROR: o2-aod-merger exited with code ${EXIT_CODE}."
    echo "Input list preserved at: ${MERGER_INPUT}"
    exit $EXIT_CODE
fi

# ------------------------------------------------------------------------------
# 4. OPTIONAL INDEX VALIDATION (ROOT macro)
# ------------------------------------------------------------------------------
# Runs a quick check that all fIndexRingCollisions values in O2ringjet and
# O2ringlav0 are within the valid range of O2ringcollision rows, for every
# DF folder in the merged file.  Requires ROOT to be in PATH.
#
# Comment out this block if you do not need the extra validation step.

if command -v root &>/dev/null; then
    echo ""
    echo "  Running index validation on merged file..."
    VALIDATION_MACRO="${MERGE_OUTPUT_DIR}/check_indices.C"

    # Write the ROOT macro on the fly.
    # Uses only the first DF folder found in the file for a quick sanity check.
    cat > "$VALIDATION_MACRO" << 'EOF'
// check_indices.C -- quick index validation for RingCollision-derived tables
// Run with: root -l -b -q 'check_indices.C("path/to/file.root")'
void check_indices(const char* fname)
{
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: cannot open %s\n", fname);
        return;
    }

    // Iterate over all DF_* keys in the file
    bool allOK = true;
    TIter next(f->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (!name.BeginsWith("DF_")) continue;

        TDirectoryFile* df = (TDirectoryFile*)f->Get(name);
        if (!df) continue;

        TTree* collisions = (TTree*)df->Get("O2ringcollision");
        if (!collisions) {
            printf("WARN : %s -- O2ringcollision not found (skip)\n", name.Data());
            continue;
        }
        long long nColl = collisions->GetEntries();

        // Tables that carry fIndexRingCollisions:
        const char* childTrees[] = {"O2ringjet", "O2ringlav0", "O2ringleadp", nullptr};
        for (int t = 0; childTrees[t] != nullptr; ++t) {
            TTree* child = (TTree*)df->Get(childTrees[t]);
            if (!child) continue; // table may be empty in some DFs -- that is OK

            int idx = -1;
            if (child->SetBranchAddress("fIndexRingCollisions", &idx) < 0) {
                printf("WARN : %s/%s -- branch fIndexRingCollisions not found\n",
                       name.Data(), childTrees[t]);
                continue;
            }
            for (long long i = 0; i < child->GetEntries(); ++i) {
                child->GetEntry(i);
                if (idx < 0 || idx >= (int)nColl) {
                    printf("BAD  : %s/%s row %lld -> idx=%d (nColl=%lld)\n",
                           name.Data(), childTrees[t], i, idx, nColl);
                    allOK = false;
                }
            }
        }
        printf("OK   : %s (nColl=%lld)\n", name.Data(), nColl);
    }
    f->Close();
    printf("\nIndex validation: %s\n", allOK ? "ALL VALID" : "ERRORS FOUND");
}
EOF

    root -l -b -q "${VALIDATION_MACRO}(\"${OUTPUT_FILE}\")" \
        2>&1 | grep -v "^Processing\|^Info\|Warning"

    rm -f "$VALIDATION_MACRO"
else
    echo "  (Skipping index validation: ROOT not found in PATH)"
fi