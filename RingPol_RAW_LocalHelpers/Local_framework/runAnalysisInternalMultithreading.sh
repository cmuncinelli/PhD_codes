#!/bin/bash

# ==============================================================================
# ALICE O2 ANALYSIS LAUNCHER: INTERNAL PARALLELISM MODE
# ==============================================================================
# Hardware Target: AMD Threadripper 3990X (64C/128T) | 256GB RAM | HDD -> NVMe
# Strategy:
#   1. Split input into large batches (~60GB data).
#   2. Copy batch sequentially to NVMe (Max HDD throughput).
#   3. Run ONE massive O2 pipeline with high internal parallelism.
#   4. Clean up immediately to free space for next batch.
# ==============================================================================
# Run as: ./runAnalysisInternalMultithreading.sh ./input_data_storage.txt ./dpl-config-JustTOFBase-FullQA_permissivePt.json ./aod-writer.json

# --- INPUT ARGUMENTS ---
INPUT_LIST="$1"
JSON_CONFIG_PATH="$2"
AOD_WRITER_JSON="$3"

# --- TUNING KNOBS (Adjust based on performance monitoring) ---
FILES_PER_BATCH=20             # 20 files * 3GB = ~60GB batch size
SHM_SIZE="128000000000"        # 128GB Shared Memory to buffer a lot of data
MEM_RATE_LIMIT="6000000000"    # 6GB/s Rate Limit

# --- PATH SETUP ---
WORK_DIR=$(pwd)
TEMP_BASE="${WORK_DIR}/temp_staging_area"
RESULTS_DIR="${WORK_DIR}/results"
LOGS_DIR="${RESULTS_DIR}/logs"
DERIVED_AOD_DIR="${RESULTS_DIR}/derived_aods"

# ==============================================================================
# 1. SIGNAL HANDLING
# ==============================================================================
# This function runs ONLY when the script exits (normally or forced)
cleanup() {
    if [ -d "$TEMP_BASE" ]; then
        echo "  [Cleanup] Removing temp staging area..."
        rm -rf "$TEMP_BASE" 2>/dev/null # We redirect to /dev/null to hide errors if folder is already gone
        echo "  [Cleanup] Staging area removed."
    fi
}

# This function runs ONLY on Ctrl+C
handle_interrupt() {
    echo ""
    echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
    echo "    Stopping current jobs..."
    pkill -P $$ 2>/dev/null # Attempt to kill child processes (the O2 pipeline) immediately
    exit 130 # Force exit. This triggers the 'cleanup' function via the EXIT trap.
}

trap cleanup EXIT # 1. Trap EXIT: Always cleans up, no matter how we leave.
trap handle_interrupt INT TERM # 2. Trap INT/TERM: Calls handle_interrupt, which then EXITS.

# ==============================================================================
# 2. INITIALIZATION & DYNAMIC NAMING
# ==============================================================================
if [[ -z "$INPUT_LIST" ]] || [[ -z "$JSON_CONFIG_PATH" ]] || [[ -z "$AOD_WRITER_JSON" ]]; then
    echo "Usage: $0 <input_list.txt> <config.json> <aod-writer.json>"
    exit 1
fi

INPUT_LIST=$(realpath "$INPUT_LIST")
JSON_CONFIG_PATH=$(realpath "$JSON_CONFIG_PATH")
AOD_WRITER_JSON=$(realpath "$AOD_WRITER_JSON")

# --- DYNAMIC FILENAME LOGIC ---
# 1. Get basename (e.g., "dpl-config-ITSEnforced.json" -> "dpl-config-ITSEnforced")
CONFIG_BASENAME=$(basename "$JSON_CONFIG_PATH" .json)

# 2. Strip "dpl-config" from the start.
# If filename is "dpl-config-ITSEnforced", SUFFIX becomes "-ITSEnforced"
SUFFIX="${CONFIG_BASENAME#dpl-config}"

echo "========================================================"
echo "  Run Configuration"
echo "  Basename: ${CONFIG_BASENAME}"
echo "  Suffix:   ${SUFFIX}"
echo "========================================================"

# Clean start
if [ -d "$TEMP_BASE" ]; then
    rm -rf "$TEMP_BASE" 2>/dev/null
fi
mkdir -p "$TEMP_BASE/batches"
mkdir -p "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR/used_configs" # To store used dpl-configs for reference
mkdir -p "$LOGS_DIR"
mkdir -p "$DERIVED_AOD_DIR"

echo "  Splitting input list..."
grep -vE '^\s*$' "$INPUT_LIST" | split -d -l "$FILES_PER_BATCH" - "$TEMP_BASE/batches/batch_"
NUM_BATCHES=$(find "$TEMP_BASE/batches" -name "batch_*" | wc -l)
echo "  Total Batches: $NUM_BATCHES"

# File to keep track of outputs for merging
MERGE_LIST_FILE="${RESULTS_DIR}/merge_list.txt"
> "$MERGE_LIST_FILE"

# ==============================================================================
# 3. MAIN LOOP
# ==============================================================================
# Iterate over batches sorted numerically
for BATCH_FILE in $(ls "$TEMP_BASE/batches/batch_"* | sort); do
    
    # Extract Batch ID (e.g., batch_01 -> 01)
    BATCH_ID="${BATCH_FILE##*_}"
    
    echo "--------------------------------------------------------"
    echo " Processing Batch $BATCH_ID / $((NUM_BATCHES-1))"
    echo "--------------------------------------------------------"

    # Define per-batch workspace
    BATCH_WORK_DIR="$TEMP_BASE/work_$BATCH_ID"
    mkdir -p "$BATCH_WORK_DIR"
    LOCAL_INPUT_LIST="$BATCH_WORK_DIR/local_list.txt"
    BATCH_LOG="$LOGS_DIR/batch_$BATCH_ID.log"

    # ---------------------------------------------------------
    # A. STAGE-IN (Sequential HDD -> NVMe)
    # ---------------------------------------------------------
    echo "  [I/O] Staging files to NVMe..."
    while read -r LINE; do
        # Handle 'file:' prefix if present, or raw paths
        CLEAN_PATH="${LINE#file:}"
        FILENAME=$(basename "$CLEAN_PATH")
        
        if [ -f "$CLEAN_PATH" ]; then
            cp "$CLEAN_PATH" "$BATCH_WORK_DIR/$FILENAME"
            echo "file:$BATCH_WORK_DIR/$FILENAME" >> "$LOCAL_INPUT_LIST"
        else
            echo "    WARNING: Source file missing: $CLEAN_PATH"
        fi
    done < "$BATCH_FILE"

    # ---------------------------------------------------------
    # B. EXECUTE PIPELINE
    # ---------------------------------------------------------
    echo "  [O2] Launching internal pipeline..."
    echo "       (Check $BATCH_LOG for live details)"
    
    cd "$BATCH_WORK_DIR" || exit 1
    
    # Common Configuration options
    # Notice you need to pass the --aod-writer argument to all tasks, otherwise we will not write the appropriate tables!
    OPTION="-b --configuration json://$JSON_CONFIG_PATH --aod-writer-json "${AOD_WRITER_JSON}""

    # Updates based on the workflow_dump.log:
    # 1. Event Selection -> eventselection-run3
    # 2. Propagation     -> propagation-service (The heavy hitter!)
    # 3. PID TPC         -> pid-tpc-service
    # 4. Strangeness     -> strangenesstofpid (Added 8 threads to prevent bottleneck)
    # 5. Lambda Jet      -> lambdajetpolarizationions
    # Total Explicit Workers: ~104
    # Estimated RAM Usage: ~120 GB
    time \
    o2-analysis-event-selection-service ${OPTION} \
        --pipeline eventselection-run3:16 | \
    o2-analysis-multcenttable ${OPTION} | \
    o2-analysis-propagationservice ${OPTION} \
        --pipeline propagation-service:32 | \
    o2-analysis-pid-tpc-service ${OPTION} \
        --pipeline pid-tpc-service:16 | \
    o2-analysis-ft0-corrected-table ${OPTION} | \
    o2-analysis-pid-tof-base ${OPTION} \
        --pipeline tof-event-time:6 | \
    o2-analysis-lf-strangenesstofpid ${OPTION} \
        --pipeline strangenesstofpid:8 | \
    o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
        --pipeline lambdajetpolarizationions:24 \
        --aod-file "@${LOCAL_INPUT_LIST}" \
        --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
        --shm-segment-size "$SHM_SIZE" \
        > "$BATCH_LOG" 2>&1
    # Notice that the instructions on which tables to write are actually contained in $JSON_CONFIG_PATH (the config json).
    # The aod-writer-json just helps with formatting! For more info, read the "internal-dpl-aod-writer" key in the .json.

    EXIT_CODE=$?

    # ---------------------------------------------------------
    # C. STAGE-OUT (Results & Derived Data)
    # ---------------------------------------------------------
    if [ $EXIT_CODE -eq 0 ]; then
        echo "  [OK] O2 finished successfully."
        
        # 1. Handle Main Results with SUFFIX
        if [ -f "AnalysisResults.root" ]; then
            # Name format: AnalysisResults-ITSEnforced_00.root
            TARGET_RES="$RESULTS_DIR/AnalysisResults${SUFFIX}_${BATCH_ID}.root"
            mv "AnalysisResults.root" "$TARGET_RES"
            echo "    Saved AnalysisResults to ${TARGET_RES}."
            echo "$TARGET_RES" >> "$MERGE_LIST_FILE"
        else
            echo "    ERROR: AnalysisResults.root not found!"
        fi

        # 2. Handle Derived AODs with SUFFIX
        # Note: Derived AODs usually have fixed names defined in the JSON writer,
        # but we are overriding that anyways via the wildcard "AO2D_*.root".
        DERIVED_FILE=$(find . -maxdepth 1 -name "AO2D_LambdaJetsRing_DerivedTest.root" | head -n 1)
        if [ -n "$DERIVED_FILE" ]; then
            # Name format: AO2D_Derived-ITSEnforced_00.root
            echo "    Moving file $DERIVED_FILE"
            mv "$DERIVED_FILE" "$DERIVED_AOD_DIR/AO2D_Derived${SUFFIX}_${BATCH_ID}.root"
            echo "    Saved derived AOD to \"$DERIVED_AOD_DIR/AO2D_Derived${SUFFIX}_${BATCH_ID}.root\"."
        else
            echo "    Warning: No derived AOD found."
        fi
    else
        echo "  [FAIL] Pipeline failed with code $EXIT_CODE"
        echo "  Check log: $BATCH_LOG"
        echo "  Skipping to next batch..."
        # This is handled by the local cleanup step below
    fi

    # ---------------------------------------------------------
    # D. LOCAL CLEANUP
    # ---------------------------------------------------------
    cd "$WORK_DIR" || exit 1
    rm -rf "$BATCH_WORK_DIR"
    echo "  [Clean] Batch workspace removed."
    echo ""

done

# ==============================================================================
# 4. FINAL MERGE & CLEANUP
# ==============================================================================
echo "========================================================"
echo " Final Step: Merging Results"
echo "========================================================"

# Final Output Name with SUFFIX (e.g., AnalysisResults-ITSEnforced.root)
FINAL_OUTPUT="$RESULTS_DIR/AnalysisResultsMerged${SUFFIX}.root"

if [ -s "$MERGE_LIST_FILE" ]; then
    hadd -f -j 16 "$FINAL_OUTPUT" @"$MERGE_LIST_FILE"
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "SUCCESS: Merge complete."
        echo "Output: $FINAL_OUTPUT"
        
        # --- Delete intermediate files ---
        echo "  [Cleanup] Deleting intermediate batch files..."
        xargs rm < "$MERGE_LIST_FILE"
        rm "$MERGE_LIST_FILE"
        
        # --- Copy Config for reproducibility ---
        cp "$JSON_CONFIG_PATH" "$RESULTS_DIR/used_configs/dpl-config${SUFFIX}.json"
    else
        echo "ERROR: hadd merge failed. Intermediate files NOT deleted."
        exit 1
    fi
else
    echo "ERROR: No results generated to merge."
    exit 1
fi