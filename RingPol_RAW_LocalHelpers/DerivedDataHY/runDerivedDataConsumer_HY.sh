#!/usr/bin/env bash

###############################################################################
# Script name:
#   runDerivedDataConsumer_HY.sh
#
# Purpose:
#   Runs the ALICE O2 derived-data consumer task on AO2D files downloaded
#   from Hyperloop, processing them in batches to avoid overwhelming the DPL
#   shared-memory scheduler on large datasets (400+ files).
#
#   For each batch the consumer produces a single AnalysisResults.root which
#   is renamed to ConsumerResults_<SUFFIX>_<BATCH_ID>.root and kept in the
#   temporary staging area. After all batches finish, hadd merges the
#   surviving results into a single final file and the staging area
#   (including all intermediate per-batch files) is deleted on exit.
#
#   Input  : <WORK_DIR>/AO2Ds/AO2D_*.root
#   Output : <WORK_DIR>/results_consumer/ConsumerResults_<CONS_SUFFIX>.root
#   Logs   : <WORK_DIR>/results_consumer/logs/batch_<BATCH_ID>_<CONS_SUFFIX>.log
#   Config : <WORK_DIR>/results_consumer/used_configs/consumer-config_<CONS_SUFFIX>.json
#
# Usage:
#   ./runDerivedDataConsumer_HY.sh <WORK_DIR> <consumer_config.json>
#
# Arguments:
#   $1 : WORK_DIR          -- absolute path to the per-wagon working directory
#                             (e.g. /home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC_min3ITS)
#   $2 : consumer_config   -- path to the consumer dpl-config JSON file
#                             (e.g. dpl-config-DerivedConsumer-JustLambda.json)
#
# Example:
#   ./runDerivedDataConsumer_HY.sh \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC_min3ITS \
#       /home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC_min3ITS/dpl-config-DerivedConsumer-JustLambda.json
#
###############################################################################

# ==============================================================================
# TUNING KNOBS
# ==============================================================================
FILES_PER_BATCH=80          # Number of AO2D files per batch. Reduce if DPL
                            # still reports scheduling stalls on this machine.
SHM_SIZE="64000000000"      # Shared memory: 64 GB
MEM_RATE_LIMIT="4000000000" # AOD memory rate limit: 4 GB/s
THREADS=32                  # Pipeline threads for the consumer task
READERS=5                   # Parallel AOD reader threads

# ==============================================================================
# 1. ARGUMENT PARSING AND VALIDATION
# ==============================================================================
WORK_DIR="$1"
JSON_CONSUMER_CONFIG="$2"

if [ -z "$WORK_DIR" ] || [ -z "$JSON_CONSUMER_CONFIG" ]; then
  echo "Usage: $0 <WORK_DIR> <consumer_config.json>"
  echo ""
  echo "  WORK_DIR        : per-wagon working directory containing AO2Ds/"
  echo "  consumer_config : dpl-config JSON for the consumer task"
  exit 1
fi

WORK_DIR="$(realpath "$WORK_DIR")"
JSON_CONSUMER_CONFIG="$(realpath "$JSON_CONSUMER_CONFIG")"

if [ ! -d "$WORK_DIR" ]; then
  echo "Error: WORK_DIR does not exist: ${WORK_DIR}"
  exit 1
fi

if [ ! -f "$JSON_CONSUMER_CONFIG" ]; then
  echo "Error: Consumer config not found: ${JSON_CONSUMER_CONFIG}"
  exit 1
fi

AOD_DIR="${WORK_DIR}/AO2Ds"

if [ ! -d "$AOD_DIR" ]; then
  echo "Error: AO2Ds/ directory not found in: ${WORK_DIR}"
  echo "Run download_hyperloop.sh first."
  exit 1
fi

# ==============================================================================
# 2. DYNAMIC NAMING
#
# Output files are named after the consumer config suffix only.
#
# Example:
#   dpl-config-DerivedConsumer-JustLambda.json  ->  CONS_SUFFIX = JustLambda
#   Final output: ConsumerResults_JustLambda.root
# ==============================================================================
CONS_BASENAME=$(basename "$JSON_CONSUMER_CONFIG" .json)
CONS_SUFFIX="${CONS_BASENAME#dpl-config-DerivedConsumer-}"

# ==============================================================================
# 3. DIRECTORY SETUP
# ==============================================================================
CONSUMER_RESULTS="${WORK_DIR}/results_consumer"
CONSUMER_LOGS="${CONSUMER_RESULTS}/logs"
CONSUMER_CONFIGS="${CONSUMER_RESULTS}/used_configs"
CONSUMER_TEMP="${WORK_DIR}/temp_consumer_stage"

mkdir -p "$CONSUMER_RESULTS"
mkdir -p "$CONSUMER_LOGS"
mkdir -p "$CONSUMER_CONFIGS"

# Remove any leftover staging area from a previous run
rm -rf "$CONSUMER_TEMP"
mkdir -p "$CONSUMER_TEMP"

# ==============================================================================
# 4. SIGNAL HANDLING
# ==============================================================================
cleanup() {
  if [ -d "$CONSUMER_TEMP" ]; then
    echo "  [Cleanup] Removing staging area: ${CONSUMER_TEMP}"
    rm -rf "$CONSUMER_TEMP" 2>/dev/null
  fi
}

handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  echo "    Stopping child processes..."
  pkill -P $$ 2>/dev/null
  exit 130
}

trap cleanup EXIT
trap handle_interrupt INT TERM

# ==============================================================================
# 5. BUILD FULL SORTED FILE LIST AND SPLIT INTO BATCHES
#
# Files are sorted by numeric index (AO2D_1.root, AO2D_2.root, ...) using
# version sort, then split into fixed-size chunks written to the staging area.
# ==============================================================================
FULL_LIST="${CONSUMER_TEMP}/all_inputs.txt"

find "$AOD_DIR" -maxdepth 1 -name "AO2D_*.root" \
  | sort -t_ -k2 -V \
  > "$FULL_LIST"

TOTAL_FILES=$(wc -l < "$FULL_LIST")

if [ "$TOTAL_FILES" -eq 0 ]; then
  echo "Error: No AO2D_*.root files found in: ${AOD_DIR}"
  exit 1
fi

# Split into batch files: batch_00, batch_01, ...
# Each batch file contains at most FILES_PER_BATCH paths (no "file:" prefix
# yet; that is added per-batch when building the O2 input list below).
split -d -l "$FILES_PER_BATCH" "$FULL_LIST" "${CONSUMER_TEMP}/batch_"

NUM_BATCHES=$(find "$CONSUMER_TEMP" -maxdepth 1 -name "batch_*" | wc -l)

echo "========================================================"
echo "  HY Derived Data Consumer (batched)"
echo "  WORK_DIR         : ${WORK_DIR}"
echo "  Consumer suffix  : ${CONS_SUFFIX}"
echo "  Total AOD files  : ${TOTAL_FILES}"
echo "  Files per batch  : ${FILES_PER_BATCH}"
echo "  Number of batches: ${NUM_BATCHES}"
echo "  Final output     : ${CONSUMER_RESULTS}/ConsumerResults_${CONS_SUFFIX}.root"
echo "========================================================"

# ==============================================================================
# 6. BATCH LOOP
# ==============================================================================
MERGE_LIST="${CONSUMER_TEMP}/merge_list.txt"
> "$MERGE_LIST"

FAILED_BATCHES=()

for BATCH_FILE in $(find "$CONSUMER_TEMP" -maxdepth 1 -name "batch_*" | sort); do

  BATCH_ID="${BATCH_FILE##*batch_}"
  BATCH_SIZE=$(wc -l < "$BATCH_FILE")

  echo "--------------------------------------------------------"
  echo "  Batch ${BATCH_ID} / $(( NUM_BATCHES - 1 ))  (${BATCH_SIZE} files)"
  echo "--------------------------------------------------------"

  # Per-batch working directory: O2 scratch files stay contained here
  BATCH_WORK_DIR="${CONSUMER_TEMP}/work_${BATCH_ID}"
  mkdir -p "$BATCH_WORK_DIR"

  # Build the O2 input list for this batch (prepend "file:" prefix)
  BATCH_INPUT_LIST="${BATCH_WORK_DIR}/inputs.txt"
  awk '{print "file:"$0}' "$BATCH_FILE" > "$BATCH_INPUT_LIST"

  BATCH_LOG="${CONSUMER_LOGS}/batch_${BATCH_ID}_${CONS_SUFFIX}.log"
  # Per-batch result files live in the temp staging area alongside the batch
  # split files. They are deleted after hadd; only the merged final output
  # lands in results_consumer/.
  BATCH_OUTPUT="${CONSUMER_TEMP}/ConsumerResults_${CONS_SUFFIX}_${BATCH_ID}.root"

  # Run the consumer from inside the batch work directory
  cd "$BATCH_WORK_DIR" || exit 1

  echo "  [O2] Launching consumer for batch ${BATCH_ID}..."

  time \
  o2-analysis-lf-lambdajetpolarizationionsderived \
      -b \
      --configuration "json://${JSON_CONSUMER_CONFIG}" \
      --pipeline "lambdajetpolarizationionsderived:${THREADS}" \
      --readers "$READERS" \
      --aod-file "@${BATCH_INPUT_LIST}" \
      --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
      --shm-segment-size "$SHM_SIZE" \
      > "$BATCH_LOG" 2>&1

  EXIT_CODE=$?

  cd "$WORK_DIR" || exit 1

  if [ $EXIT_CODE -eq 0 ]; then
    if [ -f "${BATCH_WORK_DIR}/AnalysisResults.root" ]; then
      mv "${BATCH_WORK_DIR}/AnalysisResults.root" "$BATCH_OUTPUT"
      echo "$BATCH_OUTPUT" >> "$MERGE_LIST"
      echo "  [OK] Batch ${BATCH_ID} done -> $(basename "$BATCH_OUTPUT")"
    else
      echo "  [WARN] Batch ${BATCH_ID}: O2 exited 0 but AnalysisResults.root not found."
      echo "         Check log: ${BATCH_LOG}"
      FAILED_BATCHES+=("$BATCH_ID")
    fi
  else
    echo "  [FAIL] Batch ${BATCH_ID} exited with code ${EXIT_CODE}. Skipping."
    echo "         Check log: ${BATCH_LOG}"
    FAILED_BATCHES+=("$BATCH_ID")
  fi

  # Remove the per-batch work directory regardless of outcome; logs are safe
  # in CONSUMER_LOGS and batch results (if any) are in CONSUMER_TEMP.
  rm -rf "$BATCH_WORK_DIR"

done

# ==============================================================================
# 7. REPORT FAILED BATCHES
# ==============================================================================
if [ ${#FAILED_BATCHES[@]} -gt 0 ]; then
  echo ""
  echo "========================================================"
  echo "  WARNING: The following batches failed or produced no output:"
  for B in "${FAILED_BATCHES[@]}"; do
    echo "    batch_${B}  ->  log: ${CONSUMER_LOGS}/batch_${B}_${CONS_SUFFIX}.log"
  done
  echo "  The final merged file will be based on the successful batches only."
  echo "========================================================"
fi

# ==============================================================================
# 8. MERGE BATCH RESULTS
# ==============================================================================
SUCCESSFUL=$(wc -l < "$MERGE_LIST")
FINAL_OUTPUT="${CONSUMER_RESULTS}/ConsumerResults_${CONS_SUFFIX}.root"

echo ""
echo "========================================================"
echo "  Merging ${SUCCESSFUL} / ${NUM_BATCHES} batch result(s)..."
echo "========================================================"

if [ "$SUCCESSFUL" -eq 0 ]; then
  echo "Error: No batches produced output. Nothing to merge."
  exit 1
fi

if [ "$SUCCESSFUL" -eq 1 ]; then
  # Only one batch succeeded: rename instead of running hadd on a single file
  mv "$(cat "$MERGE_LIST")" "$FINAL_OUTPUT"
else
  hadd -f "$FINAL_OUTPUT" $(cat "$MERGE_LIST")

  if [ $? -ne 0 ]; then
    echo "Error: hadd failed. Intermediate batch files NOT deleted."
    exit 1
  fi
fi

# Delete intermediate per-batch result files (keep only the merged output)
while IFS= read -r BATCH_RESULT || [ -n "$BATCH_RESULT" ]; do
  [ -f "$BATCH_RESULT" ] && rm -f "$BATCH_RESULT"
done < "$MERGE_LIST"

# Save a copy of the config used for this run
cp "$JSON_CONSUMER_CONFIG" "${CONSUMER_CONFIGS}/consumer-config_${CONS_SUFFIX}.json"

echo ""
echo "SUCCESS."
echo "  Output  : ${FINAL_OUTPUT}"
echo "  Logs    : ${CONSUMER_LOGS}/"
echo "  Config  : ${CONSUMER_CONFIGS}/consumer-config_${CONS_SUFFIX}.json"
if [ ${#FAILED_BATCHES[@]} -gt 0 ]; then
  echo "  NOTE    : ${#FAILED_BATCHES[@]} batch(es) failed -- results are partial."
fi
