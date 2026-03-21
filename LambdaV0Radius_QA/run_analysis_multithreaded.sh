#!/bin/bash

# ==============================================================================
# ALICE O2 MULTI-THREADED STAGED (HDD -> SSD Cache) ANALYSIS LAUNCHER
# ==============================================================================
# Usage: ./run_analysis_multithreaded.sh <input_list.txt> <config.json> [num_jobs]
# Example: ./run_analysis_multithreaded.sh input_data_storage.txt dpl-config.json 16
# The input_list.txt file should have the prefix file:/ for each file as in
# gen_input_paths.sh!
# ==============================================================================

INPUT_LIST="$1"
JSON_CONFIG_PATH="$2"
MAX_JOBS="${3:-16}" # Number of CPU threads to use

# --- TUNING PARAMETERS ---
SHM_SIZE="4000000000"        # 4GB Shared Memory (Enough for most AO2Ds)
MEM_RATE_LIMIT="500000000"   # 500MB Rate Limit
FILES_PER_CHUNK=5            # Set this to 5 to save SSD space (5 * 20 jobs = 100 active files, or about 60 GB)

if [[ -z "$INPUT_LIST" ]] || [[ -z "$JSON_CONFIG_PATH" ]]; then
    echo "Usage: $0 <input_list.txt> <config.json> [num_jobs]"
    exit 1
fi

INPUT_LIST=$(realpath "$INPUT_LIST")
JSON_CONFIG_PATH=$(realpath "$JSON_CONFIG_PATH")
WORK_DIR=$(pwd)
TEMP_BASE="${WORK_DIR}/temp_staging_area"
RESULTS_DIR="${WORK_DIR}/results"
LOGS_DIR="${RESULTS_DIR}/logs"  # <--- NEW: Permanent Log Directory

# --- SAFETY TRAP ---
# This ensures the temp folder is ALWAYS deleted, even if you Ctrl+C or the script crashes.
cleanup() {
    # Does not delete logs though! Be aware of that!
    if [ -d "$TEMP_BASE" ]; then
        echo ""
        echo "  [Cleanup] Removing temp staging area..."
        rm -rf "$TEMP_BASE"
    fi
}
# Trigger cleanup on Exit (0/1), Interrupt (Ctrl+C), or Terminate signals
trap cleanup EXIT INT TERM

echo "========================================================"
echo "   Starting O2 Staged Analysis"
echo "   Strategy: Copy Batch -> Process -> Delete"
echo "   Threads: $MAX_JOBS"
echo "   Batch Size: $FILES_PER_CHUNK files"
echo "   Logs: $LOGS_DIR"
echo "========================================================"

# --- 1. CLEANUP & SETUP ---

if [ -d "$TEMP_BASE" ]; then
    rm -rf "$TEMP_BASE"
fi

# FIX: Clean old logs to prevent appending to previous runs
if [ -d "$LOGS_DIR" ]; then
    echo "   Cleaning old logs..."
    rm -rf "$LOGS_DIR"
fi

mkdir -p "$TEMP_BASE/lists"
mkdir -p "$TEMP_BASE/outputs"
mkdir -p "$RESULTS_DIR"
mkdir -p "$LOGS_DIR" # Ensure log dir exists

# --- 2. SPLIT INPUT ---

echo "   Splitting input list..."
grep -vE '^\s*$' "$INPUT_LIST" | split -l "$FILES_PER_CHUNK" - "$TEMP_BASE/lists/batch_"
NUM_BATCHES=$(ls "$TEMP_BASE/lists/batch_"* | wc -l)
echo "   Created $NUM_BATCHES batches."

# --- 3. THE WORKER ---

run_staged_job() {
    BATCH_LIST="$1"
    JOB_ID="$2"
    JSON="$3"
    BASE_TMP="$4"
    SHM="$5"
    MEM_LIMIT="$6"
    LOG_DIR="$7" # Pass the log dir

    # Define Log File in the PERMANENT directory
    LOG_FILE="${LOG_DIR}/job_${JOB_ID}.log"

    # Define Isolated SSD Workspace
    JOB_DIR="${BASE_TMP}/worker_${JOB_ID}"
    DATA_DIR="${JOB_DIR}/data"
    mkdir -p "$DATA_DIR"
    
    LOCAL_LIST_FILE="${JOB_DIR}/local_input.txt"
    touch "$LOCAL_LIST_FILE"
    
    # STAGE-IN
    while IFS= read -r RAW_LINE; do
        HDD_PATH="${RAW_LINE#"file:"}"
        FILENAME=$(basename "$HDD_PATH")
        if [ -f "$HDD_PATH" ]; then
            cp "$HDD_PATH" "${DATA_DIR}/${FILENAME}"
            echo "file:${DATA_DIR}/${FILENAME}" >> "$LOCAL_LIST_FILE"
        else
            echo "Warning: Source file not found: $HDD_PATH" >> "$LOG_FILE"
        fi
    done < "$BATCH_LIST"

    # EXECUTE O2
    cd "$JOB_DIR" || return 1
    OPTION="-b --configuration json://${JSON}"

    # We redirect output to the persistent log file
    o2-analysis-lf-strangenesstofpid ${OPTION} | \
    o2-analysis-lf-v0mlscoresconverter ${OPTION} | \
    o2-analysis-lf-asymmetric-rapidity-test ${OPTION} \
        --aod-file "@${LOCAL_LIST_FILE}" \
        --aod-memory-rate-limit "$MEM_LIMIT" \
        --shm-segment-size "$SHM" \
        >> "$LOG_FILE" 2>&1

    RC=$?

    # STAGE-OUT
    if [ $RC -eq 0 ] && [ -f "AnalysisResults.root" ]; then
        mv "AnalysisResults.root" "${BASE_TMP}/outputs/AnalysisResults_${JOB_ID}.root"
        STATUS="OK"
    else
        STATUS="Error!!!"
        echo "   !!! JOB FAILED. See log: $LOG_FILE" # Print to main screen
    fi

    # CLEANUP DATA ONLY
    rm -rf "$JOB_DIR"
    echo "   [Batch $JOB_ID] $STATUS"
}

export -f run_staged_job

# --- 4. EXECUTE ---

echo "  Processing queue... (Logs are in $LOGS_DIR)"
find "$TEMP_BASE/lists" -name "batch_*" | parallel --progress --eta -j "$MAX_JOBS" \
    run_staged_job {} {#} "$JSON_CONFIG_PATH" "$TEMP_BASE" "$SHM_SIZE" "$MEM_RATE_LIMIT" "$LOGS_DIR"
# # --bar provides a cleaner visual bar than --progress
# find "$TEMP_BASE/lists" -name "batch_*" | parallel --bar -j "$MAX_JOBS" \
#     run_staged_job {} {#} "$JSON_CONFIG_PATH" "$TEMP_BASE" "$SHM_SIZE" "$MEM_RATE_LIMIT" "$LOGS_DIR"

# --- 5. MERGE ---

echo "Merging Results..."
TARGET_FILE="${RESULTS_DIR}/AnalysisResults.root"

if ls "$TEMP_BASE/outputs"/AnalysisResults_*.root 1> /dev/null 2>&1; then
    hadd -f -k -j "$MAX_JOBS" "$TARGET_FILE" "$TEMP_BASE/outputs"/AnalysisResults_*.root
    mv "$PWD/dpl-config.json" results/Asymmetry.json
    echo "  Done. Final file: $TARGET_FILE"
else
    echo "!!! No results found to merge."
    echo "   Check the log files in $LOGS_DIR for details."
    exit 1
fi