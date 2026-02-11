#!/bin/bash

# ==============================================================================
# ALICE O2 MULTI-THREADED STAGED (HDD -> SSD Cache) ANALYSIS LAUNCHER
# ==============================================================================
# Usage: ./run_analysis_multithreaded.sh <input_list.txt> <config.json> [num_jobs]
# Example: ./run_analysis_multithreaded.sh input_data_storage.txt dpl-config.json 16
# The input_list.txt file should have the prefix file:/ for each file as in
# gen_input_paths.sh!
# ==============================================================================
# Run to get the dpl-config dummy file: o2-analysis-lf-asymmetric-rapidity-test | o2-analysis-lf-strangenesstofpid | o2-analysis-lf-v0mlscoresconverter --aod-file /storage2/cicero/asymmetric_rapidity_V0tests_LHC23_pass4_Thin_small/AO2D_1.root -b

INPUT_LIST="$1"
JSON_CONFIG_PATH="$2"
AOD_WRITER_JSON="$3"
MAX_JOBS="${4:-4}" # Number of CPU threads to use. Uses only 4 due to small jarvis4 available disk space.

# --- TUNING PARAMETERS ---
SHM_SIZE="12000000000"        # 12GB Shared Memory (Increased to deal with RAW AO2Ds and many pipelines)
MEM_RATE_LIMIT="600000000"   # 600MB Rate Limit
FILES_PER_CHUNK=5            # Set this to 5 to save SSD space (5 * 4 jobs = 20 active files of 3GB each (raw data), or about 60 GB)

if [[ -z "$INPUT_LIST" ]] || [[ -z "$JSON_CONFIG_PATH" ]] || [[ -z "$AOD_WRITER_JSON" ]]; then
    echo "Usage: $0 <input_list.txt> <config.json> <aod-writer.json> [num_jobs]"
    exit 1
fi

INPUT_LIST=$(realpath "$INPUT_LIST")
JSON_CONFIG_PATH=$(realpath "$JSON_CONFIG_PATH")
AOD_WRITER_JSON=$(realpath "$AOD_WRITER_JSON")
WORK_DIR=$(pwd)
TEMP_BASE="${WORK_DIR}/temp_staging_area"
RESULTS_DIR="${WORK_DIR}/results"
LOGS_DIR="${RESULTS_DIR}/logs"  # <--- NEW: Permanent Log Directory
DERIVED_AOD_DIR="${RESULTS_DIR}/derived_aods" # New folder, to save derived aod data!

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
mkdir -p "$DERIVED_AOD_DIR"

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
    AOD_WRITER_JSON="$4"
    BASE_TMP="$5"
    SHM="$6"
    MEM_LIMIT="$7"
    LOG_DIR="$8" # Pass the log dir
    DERIVED_AOD_DIR="$9"

    # Define Log File in the PERMANENT directory
    LOG_FILE="${LOG_DIR}/job_${JOB_ID}.log"

    # Define Isolated SSD Workspace
    JOB_DIR="${BASE_TMP}/worker_${JOB_ID}"
    DATA_DIR="${JOB_DIR}/data"
    mkdir -p "$DATA_DIR"
    
    LOCAL_LIST_FILE="${JOB_DIR}/local_input.txt"
    touch "$LOCAL_LIST_FILE"
    
    # ============================
    # STAGE-IN (Serialized HDD I/O)
    # ============================
    # Made something like an OMP_CRITICAL block, so that multiple parallel workers can't rush to copy files.
    # The throughtput is larger when running sequential reads instead of MAX_JOBS calls of cp, which will read
    # as if random access reads (thrashing), not sequential reads. This should increase the speed of copy before 
    # processing truly starts. At this stage of parallelization, it seemed like a true bottleneck to be solved.
    LOCKFILE="${BASE_TMP}/hdd_copy.lock"
    {
        echo "[Job ${JOB_ID}] Waiting for HDD lock..." >> "$LOG_FILE"
        flock -x 200
        echo "[Job ${JOB_ID}] Acquired HDD lock. Starting copy..." >> "$LOG_FILE"

        while IFS= read -r RAW_LINE; do
            HDD_PATH="${RAW_LINE#"file:"}"
            FILENAME=$(basename "$HDD_PATH")

            if [ -f "$HDD_PATH" ]; then
                # # Lower disk priority to reduce interference
                # ionice -c2 -n7 cp "$HDD_PATH" "${DATA_DIR}/${FILENAME}"
                cp "$HDD_PATH" "${DATA_DIR}/${FILENAME}"
                echo "file:${DATA_DIR}/${FILENAME}" >> "$LOCAL_LIST_FILE"
            else
                echo "Warning: Source file not found: $HDD_PATH" >> "$LOG_FILE"
            fi
        done < "$BATCH_LIST"
        echo "[Job ${JOB_ID}] Copy finished. Releasing HDD lock." >> "$LOG_FILE"
    } 200>"$LOCKFILE"

    # EXECUTE O2
    cd "$JOB_DIR" || return 1
    OPTION="-b --configuration json://${JSON}"

    o2-analysis-event-selection-service ${OPTION} | \
    o2-analysis-multcenttable ${OPTION} | \
    o2-analysis-propagationservice ${OPTION} | \
    o2-analysis-pid-tpc-service ${OPTION} | \
    o2-analysis-ft0-corrected-table ${OPTION} | \
    o2-analysis-pid-tof-base ${OPTION} | \
    o2-analysis-lf-strangenesstofpid ${OPTION} | \
    o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
        --aod-writer-json "${AOD_WRITER_JSON}" \
        --aod-file "@${LOCAL_LIST_FILE}" \
        --aod-memory-rate-limit "$MEM_LIMIT" \
        --shm-segment-size "$SHM" \
        >> "$LOG_FILE" 2>&1

    RC=$?

    # STAGE-OUT
    if [ $RC -eq 0 ]; then
        STATUS="OK"

        # --- Physics output ---
        if [ -f "AnalysisResults.root" ]; then
            mv "AnalysisResults.root" \
               "${BASE_TMP}/outputs/AnalysisResults_${JOB_ID}.root"
        else
            echo "Warning: AnalysisResults.root missing (job ${JOB_ID})" >> "$LOG_FILE"
        fi

        # --- Derived AOD output ---
        DERIVED_AOD_NAME="AO2D_LambdaJetsRing_DerivedTest.root"
        if [ -f "${DERIVED_AOD_NAME}" ]; then
            mv "${DERIVED_AOD_NAME}" \
               "${DERIVED_AOD_DIR}/AO2D_LambdaJetsRing_job${JOB_ID}.root"
        else
            echo "Warning: Derived AOD missing (job ${JOB_ID})" >> "$LOG_FILE"
        fi
    else
        STATUS="Error!!!"
        echo "   !!! JOB FAILED. See log: $LOG_FILE"
    fi

    # CLEANUP DATA ONLY
    rm -rf "$JOB_DIR"
    echo "   [Batch $JOB_ID] $STATUS"
}

export -f run_staged_job

# --- 4. EXECUTE ---

echo "  Processing queue... (Logs are in $LOGS_DIR)"
find "$TEMP_BASE/lists" -name "batch_*" | parallel --progress --eta -j "$MAX_JOBS" \
    run_staged_job {} {#} "$JSON_CONFIG_PATH" "$AOD_WRITER_JSON" \
    "$TEMP_BASE" "$SHM_SIZE" "$MEM_RATE_LIMIT" "$LOGS_DIR" "$DERIVED_AOD_DIR"
# # --bar provides a cleaner visual bar than --progress
# find "$TEMP_BASE/lists" -name "batch_*" | parallel --bar -j "$MAX_JOBS" \
#     run_staged_job {} {#} "$JSON_CONFIG_PATH" "$TEMP_BASE" "$SHM_SIZE" "$MEM_RATE_LIMIT" "$LOGS_DIR"

# --- 5. MERGE ---

echo "   Merging Results..."

# --- DYNAMIC FILENAME LOGIC ---
# 1. Get basename (e.g., "dpl-config-ITSEnforced")
CONFIG_BASENAME=$(basename "$JSON_CONFIG_PATH" .json)

echo "  Using basename ${CONFIG_BASENAME}"

# 2. Strip "dpl-config" from the start.
# If filename is "dpl-config-ITSEnforced", SUFFIX becomes "-ITSEnforced"
SUFFIX="${CONFIG_BASENAME#dpl-config}"

# 3. Construct target filename
TARGET_FILE="${RESULTS_DIR}/AnalysisResults${SUFFIX}.root"

if ls "$TEMP_BASE/outputs"/AnalysisResults_*.root 1> /dev/null 2>&1; then
    hadd -f -k -j "$MAX_JOBS" "$TARGET_FILE" "$TEMP_BASE/outputs"/AnalysisResults_*.root
    # Also copy the JSON with a matching name for records:
    cp "$JSON_CONFIG_PATH" "${RESULTS_DIR}/dpl-config${SUFFIX}.json" 2>/dev/null
    echo "  Done. Final file: $TARGET_FILE"

    #####################################################################################################
    ### NOTICE THERE IS NO NEED TO HADD THE DERIVED DATA!!! O2 ALREADY EXPECTS AN OUTPUT LIKE THAT!!! ###
    #####################################################################################################
else
    echo "!!! No results found to merge."
    echo "   Check the log files in $LOGS_DIR for details."
    exit 1
fi