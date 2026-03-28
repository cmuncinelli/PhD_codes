#!/bin/bash
###############################################################################
# Script name:
#   runProducerAsHY.sh
#
# Purpose:
#   Runs the TableProducer locally, mimicking a Hyperloop (HY) download.
#   It dynamically creates a "wagon" structure in the specified dataset folder,
#   stages data from a slow LAN to the local fast SSD in batches, processes it 
#   with massive internal parallelism, and outputs structured, derived AO2Ds.
#
# Hardware Target: 2x AMD EPYC 9634 84-Core Processor, 512 GB RAM, 
# SSD is Micron_7450_MTFDKBG3T8TFR, with a whopping 5000 MB/s read speed! (Thus I increased the number of workers)
# Strategy:
#   1. Stage LAN data sequentially to SSD/NVMe (overcoming I/O bottlenecks).
#   2. Blast through the batch using O2 with scaled-up threads & large SHM.
#   3. Save derived AO2Ds and QA directly into the HY-style folder structure.
#   4. Auto-generate input_data_storage.txt for immediate consumer use.
#
# Usage:
#   ./runProducerAsHY.sh [--TOF | --no-TOF] <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG_JSON> <AOD_WRITER_JSON>
###############################################################################

# --- INPUT ARGUMENTS & FLAG PARSING ---
USE_TOF="auto"
LOG_OUTPUT=false
OUTPUT_DATASET_DIR=""
INPUT_LIST=""
JSON_CONFIG_PATH=""
AOD_WRITER_JSON=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --TOF)
      USE_TOF=true
      shift
      ;;
    --no-TOF)
      USE_TOF=false
      shift
      ;;
    --do-log)
      LOG_OUTPUT=true
      shift
      ;;
    *)
      if [[ -z "$OUTPUT_DATASET_DIR" ]]; then OUTPUT_DATASET_DIR="$1"
      elif [[ -z "$INPUT_LIST" ]]; then INPUT_LIST="$1"
      elif [[ -z "$JSON_CONFIG_PATH" ]]; then JSON_CONFIG_PATH="$1"
      elif [[ -z "$AOD_WRITER_JSON" ]]; then AOD_WRITER_JSON="$1"
      fi
      shift
      ;;
  esac
done

if [[ -z "$OUTPUT_DATASET_DIR" ]] || [[ -z "$INPUT_LIST" ]] || [[ -z "$JSON_CONFIG_PATH" ]] || [[ -z "$AOD_WRITER_JSON" ]]; then
    echo "Usage: $0 [--TOF | --no-TOF] [--no-log] <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG_JSON> <AOD_WRITER_JSON>"
    exit 1
fi

OUTPUT_DATASET_DIR=$(realpath "$OUTPUT_DATASET_DIR")
INPUT_LIST=$(realpath "$INPUT_LIST")
JSON_CONFIG_PATH=$(realpath "$JSON_CONFIG_PATH")
AOD_WRITER_JSON=$(realpath "$AOD_WRITER_JSON")

# ==============================================================================
# 1. DYNAMIC HYPERLOOP STRUCTURE SETUP
# ==============================================================================
# Extract WAGON_SHORTNAME (e.g., "dpl-config-ITSandTPC.json" -> "ITSandTPC")
CONFIG_BASENAME=$(basename "$JSON_CONFIG_PATH" .json)
WAGON_SHORTNAME="${CONFIG_BASENAME#dpl-config-}"
WAGON_DIR="${OUTPUT_DATASET_DIR}/${WAGON_SHORTNAME}"

echo "========================================================"
echo "  Local Hyperloop Producer Initializing"
echo "  Dataset Dir : $OUTPUT_DATASET_DIR"
echo "  Wagon Name  : $WAGON_SHORTNAME"
echo "========================================================"

# Create the HY-style folder tree
mkdir -p "$WAGON_DIR/AO2Ds"
mkdir -p "$WAGON_DIR/results_producer/logs"
mkdir -p "$WAGON_DIR/results_producer/used_configs"

# Working/Staging area setup (kept inside the wagon folder to isolate runs)
TEMP_BASE="${WAGON_DIR}/temp_staging_area"
if [ -d "$TEMP_BASE" ]; then
    rm -rf "$TEMP_BASE" 2>/dev/null
fi
mkdir -p "$TEMP_BASE/batches"

# Copy the config for reproducibility
cp "$JSON_CONFIG_PATH" "$WAGON_DIR/results_producer/used_configs/"

# --- TOF CONFIGURATION LOGIC ---
if [[ "$USE_TOF" == "auto" ]]; then
    if grep -E -q '"processDataWithTOF"\s*:\s*("true"|true)' "$JSON_CONFIG_PATH"; then
        USE_TOF=true
        TOF_SOURCE="JSON Config"
    else
        USE_TOF=false
        TOF_SOURCE="JSON Config"
    fi
else
    TOF_SOURCE="Manual CLI Flag"
fi

echo "  TOF Pipeline: ${USE_TOF^^} (Source: $TOF_SOURCE)"
if [ "$LOG_OUTPUT" = false ]; then
    echo "  File Logging: DISABLED (Sending all output to /dev/null)"
fi

# --- TUNING KNOBS ---
# jarvis15 NUMA sweet spot: 20 files (~60 GB) perfectly balances Node 0's 257 GB RAM limit.
# It allows Linux to hold the entire LAN-to-SSD transfer in the Page Cache for blazing-fast 
# in-memory reads, without competing with our 128 GB SHM and forcing a cross-socket memory spill.
# (remember that in run_all_producers.sh I am running with "numactl --cpunodebind=0 --preferred=0" options!)
FILES_PER_BATCH=20             # debug at 10 is fine
SHM_SIZE="128000000000"        # 128 GB
MEM_RATE_LIMIT="8000000000"    # 8 GB/s, scaled for larger batches

# --- DETECT HARDWARE ---
TOTAL_CORES=$(lscpu | awk '/^Core\(s\) per socket:/ {cores=$4} /^Socket\(s\):/ {sockets=$2} END {print cores*sockets}')
TOTAL_THREADS=$(nproc)

# --- LIMIT: USE ONLY HALF OF CORES (IMPORTANT) ---
MAX_CORES_ALLOWED=$((TOTAL_CORES / 2))
MAX_THREADS_ALLOWED=$((MAX_CORES_ALLOWED * 2))   # assuming SMT (2 threads/core)

# # --- PIPELINE PARALLELISM ---
# # Rule of thumb: pid-tpc >= propagation (it is the CPU bottleneck).
# # Eventsel is lightweight (flag checking only), keep it moderate.
# # Lambda is the analysis task itself; keep it high to exploit parallelism.
# # Changes from previous iteration:
# PIPE_EVENTSEL=10
# PIPE_PROPAGATION=22
# PIPE_PIDTPC=28        # This is a real bottleneck in the analysis
# PIPE_MULTCENT=12      # Another bottleneck we had
# PIPE_TOF=8
# PIPE_STRANGE=8
# PIPE_LAMBDA=20

# --- PIPELINE PARALLELISM (REVISED DOWN) ---
# The sys > user timing proved the OS was overwhelmed by ZMQ inter-process 
# communication overhead. We are reducing the total process count drastically.
PIPE_EVENTSEL=4
PIPE_PROPAGATION=8
PIPE_PIDTPC=10
PIPE_MULTCENT=4
PIPE_TOF=4
PIPE_STRANGE=4
PIPE_LAMBDA=8

USED_THREADS=$((PIPE_EVENTSEL + PIPE_PROPAGATION + PIPE_PIDTPC + PIPE_MULTCENT + PIPE_TOF + PIPE_STRANGE + PIPE_LAMBDA))

# --- I/O TUNING (REVISED DOWN) ---
# Readers scale with FILES_PER_BATCH. Cap safeguard below handles overflow.
# IO_THREADS scales with readers: ~1.5x is a reasonable write-side ratio.
REQUESTED_READERS=8
IO_THREADS=8
# SPAWNERS=6 # Not using spawners explicitly. It fails to configure the links in a heavily pipelined workflow, as it seems
             # It would appear as " --spawners "$SPAWNERS" " in the pipeline, though

# --- SAFEGUARD: CAP READERS TO MAX FILES PER BATCH ---
# If we don't do this, program crashes!!!
if [ "$REQUESTED_READERS" -gt "$FILES_PER_BATCH" ]; then
    echo "  [I/O Tuning] WARNING: Requested readers ($REQUESTED_READERS) exceeds files per batch ($FILES_PER_BATCH). Capping READERS to $FILES_PER_BATCH."
    READERS=$FILES_PER_BATCH
else
    READERS=$REQUESTED_READERS
fi

echo "  Resource configuration:"
echo "    Total cores        : $TOTAL_CORES"
echo "    Total threads      : $TOTAL_THREADS"
echo "    Max usable cores   : $MAX_CORES_ALLOWED"
echo "    Max usable threads : $MAX_THREADS_ALLOWED"
echo "    Used threads       : $USED_THREADS"
echo "    Readers            : $READERS"
echo "    I/O Threads        : $IO_THREADS"
# echo "    Spawners           : $SPAWNERS"

# ==============================================================================
# 2. SIGNAL HANDLING & CLEANUP
# ==============================================================================
cleanup() {
    if [ -d "$TEMP_BASE" ]; then
        echo "  [Cleanup] Removing temp staging area..."
        rm -rf "$TEMP_BASE" 2>/dev/null
    fi
}

handle_interrupt() {
    echo ""
    echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
    echo "    Stopping current jobs..."
    pkill -P $$ 2>/dev/null
    exit 130
}

trap cleanup EXIT 
trap handle_interrupt INT TERM 

# ==============================================================================
# 3. PREPARE BATCHES
# ==============================================================================
echo "  Splitting input list..."
grep -vE '^\s*$' "$INPUT_LIST" | split -d -l "$FILES_PER_BATCH" - "$TEMP_BASE/batches/batch_"
NUM_BATCHES=$(find "$TEMP_BASE/batches" -name "batch_*" | wc -l)
echo "  Total Batches: $NUM_BATCHES"

# File to keep track of outputs for merging
MERGE_LIST_FILE="${WAGON_DIR}/results_producer/merge_list.txt"
> "$MERGE_LIST_FILE"

# ==============================================================================
# 4. MAIN LOOP
# ==============================================================================
for BATCH_FILE in $(ls "$TEMP_BASE/batches/batch_"* | sort); do
    
    # Extract Batch ID (e.g., batch_01 -> 01)
    BATCH_ID="${BATCH_FILE##*_}"
    
    echo "--------------------------------------------------------"
    echo " Processing Batch $BATCH_ID / $(printf "%02d" $((NUM_BATCHES-1)))"
    echo "--------------------------------------------------------"

    BATCH_WORK_DIR="$TEMP_BASE/work_$BATCH_ID"
    mkdir -p "$BATCH_WORK_DIR"
    LOCAL_INPUT_LIST="$BATCH_WORK_DIR/local_list.txt"
    BATCH_LOG="$WAGON_DIR/results_producer/logs/batch_$BATCH_ID.log"

    # ---------------------------------------------------------
    # A. STAGE-IN (1 Gbps LAN -> Local SSD)
    # ---------------------------------------------------------
    echo "  [I/O] Staging files from LAN to local SSD..."

    TOTAL_FILES=$(wc -l < "$BATCH_FILE")
    CURRENT_FILE=0

    while read -r LINE; do
        CLEAN_PATH="${LINE#file:}"
        FILENAME=$(basename "$CLEAN_PATH")

        ((CURRENT_FILE++))

        if [ -f "$CLEAN_PATH" ]; then
            cp "$CLEAN_PATH" "$BATCH_WORK_DIR/$FILENAME"
            echo "file:$BATCH_WORK_DIR/$FILENAME" >> "$LOCAL_INPUT_LIST"
        else
            echo -e "\n    WARNING: Source file missing: $CLEAN_PATH"
        fi

        # --- Progress bar ---
        PERCENT=$((100 * CURRENT_FILE / TOTAL_FILES))
        BAR_WIDTH=30
        FILLED=$((PERCENT * BAR_WIDTH / 100))
        EMPTY=$((BAR_WIDTH - FILLED))

        printf "\r  Progress: [%.*s%*s] %3d%% (%d/%d)" \
            "$FILLED" "##############################" \
            "$EMPTY" "" \
            "$PERCENT" "$CURRENT_FILE" "$TOTAL_FILES"

    done < "$BATCH_FILE"

    echo ""  # newline after progress bar

    # ---------------------------------------------------------
    # B. EXECUTE PIPELINE (Scaled safely for shared system)
    # ---------------------------------------------------------
    echo "  [O2] Launching controlled parallel pipeline..."
    if [ "$LOG_OUTPUT" = true ]; then
        echo "       (Check logs at: results_producer/logs/batch_$BATCH_ID.log)"
    else
        echo "       (Logging disabled. Killing output)"
    fi

    cd "$BATCH_WORK_DIR" || exit 1

    OPTION="-b --configuration json://$JSON_CONFIG_PATH --aod-writer-json ${AOD_WRITER_JSON}"

    if [ "$USE_TOF" = true ]; then
        if [ "$LOG_OUTPUT" = true ]; then
            time o2-analysis-event-selection-service ${OPTION} \
                --pipeline eventselection-run3:${PIPE_EVENTSEL} | \
            o2-analysis-multcenttable ${OPTION} \
                --pipeline mult-cent-table:${PIPE_MULTCENT} | \
            o2-analysis-propagationservice ${OPTION} \
                --pipeline propagation-service:${PIPE_PROPAGATION} | \
            o2-analysis-pid-tpc-service ${OPTION} \
                --pipeline pid-tpc-service:${PIPE_PIDTPC} | \
            o2-analysis-ft0-corrected-table ${OPTION} | \
            o2-analysis-pid-tof-base ${OPTION} \
                --pipeline tof-event-time:${PIPE_TOF} | \
            o2-analysis-lf-strangenesstofpid ${OPTION} \
                --pipeline strangenesstofpid:${PIPE_STRANGE} | \
            o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
                --pipeline lambdajetpolarizationions:${PIPE_LAMBDA} \
                --readers "$READERS" \
                --io-threads "$IO_THREADS" \
                --aod-file "@${LOCAL_INPUT_LIST}" \
                --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
                --shm-segment-size "$SHM_SIZE" \
                > "$BATCH_LOG" 2>&1
        else
            time o2-analysis-event-selection-service ${OPTION} \
                --pipeline eventselection-run3:${PIPE_EVENTSEL} | \
            o2-analysis-multcenttable ${OPTION} \
                --pipeline mult-cent-table:${PIPE_MULTCENT} | \
            o2-analysis-propagationservice ${OPTION} \
                --pipeline propagation-service:${PIPE_PROPAGATION} | \
            o2-analysis-pid-tpc-service ${OPTION} \
                --pipeline pid-tpc-service:${PIPE_PIDTPC} | \
            o2-analysis-ft0-corrected-table ${OPTION} | \
            o2-analysis-pid-tof-base ${OPTION} \
                --pipeline tof-event-time:${PIPE_TOF} | \
            o2-analysis-lf-strangenesstofpid ${OPTION} \
                --pipeline strangenesstofpid:${PIPE_STRANGE} | \
            o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
                --pipeline lambdajetpolarizationions:${PIPE_LAMBDA} \
                --readers "$READERS" \
                --io-threads "$IO_THREADS" \
                --aod-file "@${LOCAL_INPUT_LIST}" \
                --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
                --shm-segment-size "$SHM_SIZE" \
                > /dev/null 2>&1
        fi
    else
        # If debugging code, insert the flag "--fairmq-rate-logging 1" in the end of the workflow.
        # This can help diagnose bottlenecks
        if [ "$LOG_OUTPUT" = true ]; then
            time o2-analysis-event-selection-service ${OPTION} \
                --pipeline eventselection-run3:${PIPE_EVENTSEL} | \
            o2-analysis-multcenttable ${OPTION} \
                --pipeline mult-cent-table:${PIPE_MULTCENT} | \
            o2-analysis-propagationservice ${OPTION} \
                --pipeline propagation-service:${PIPE_PROPAGATION} | \
            o2-analysis-pid-tpc-service ${OPTION} \
                --pipeline pid-tpc-service:${PIPE_PIDTPC} | \
            o2-analysis-ft0-corrected-table ${OPTION} | \
            o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
                --pipeline lambdajetpolarizationions:${PIPE_LAMBDA} \
                --readers "$READERS" \
                --io-threads "$IO_THREADS" \
                --aod-file "@${LOCAL_INPUT_LIST}" \
                --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
                --shm-segment-size "$SHM_SIZE" \
                > "$BATCH_LOG" 2>&1
        else
            time o2-analysis-event-selection-service ${OPTION} \
                --pipeline eventselection-run3:${PIPE_EVENTSEL} | \
            o2-analysis-multcenttable ${OPTION} \
                --pipeline mult-cent-table:${PIPE_MULTCENT} | \
            o2-analysis-propagationservice ${OPTION} \
                --pipeline propagation-service:${PIPE_PROPAGATION} | \
            o2-analysis-pid-tpc-service ${OPTION} \
                --pipeline pid-tpc-service:${PIPE_PIDTPC} | \
            o2-analysis-ft0-corrected-table ${OPTION} | \
            o2-analysis-lf-lambdajetpolarizationions ${OPTION} \
                --pipeline lambdajetpolarizationions:${PIPE_LAMBDA} \
                --readers "$READERS" \
                --io-threads "$IO_THREADS" \
                --aod-file "@${LOCAL_INPUT_LIST}" \
                --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
                --shm-segment-size "$SHM_SIZE" \
                > /dev/null 2>&1
        fi
    fi

    EXIT_CODE=$?

    # ---------------------------------------------------------
    # C. STAGE-OUT (Results & Derived Data)
    # ---------------------------------------------------------
    if [ $EXIT_CODE -eq 0 ]; then
        echo "  [OK] O2 finished successfully."
        
        # 1. Handle QA AnalysisResults
        if [ -f "AnalysisResults.root" ]; then
            TARGET_RES="$WAGON_DIR/results_producer/AnalysisResults_${BATCH_ID}.root"
            mv "AnalysisResults.root" "$TARGET_RES"
            echo "$TARGET_RES" >> "$MERGE_LIST_FILE"
        fi

        # 2. Handle Derived AODs -> Placed in WAGON_DIR/AO2Ds/
        DERIVED_FILE=$(find . -maxdepth 1 -name "AO2D_*.root" | head -n 1)
        if [ -n "$DERIVED_FILE" ]; then
            TARGET_AOD="$WAGON_DIR/AO2Ds/AO2D_${BATCH_ID}.root"
            mv "$DERIVED_FILE" "$TARGET_AOD"
            echo "    Saved derived AOD to: AO2Ds/AO2D_${BATCH_ID}.root"
        else
            echo "    Warning: No derived AOD found in this batch."
        fi
    else
        echo "  [FAIL] Pipeline failed with code $EXIT_CODE"
        if [ "$LOG_OUTPUT" = true ]; then
            echo "  Check log: $BATCH_LOG"
        fi
    fi

    # ---------------------------------------------------------
    # D. LOCAL CLEANUP
    # ---------------------------------------------------------
    cd "$WAGON_DIR" || exit 1
    rm -rf "$BATCH_WORK_DIR"
    echo "  [Clean] Batch workspace deleted to free SSD space."
    echo ""

done

# ==============================================================================
# 5. FINAL MERGE & HYPERLOOP INTEGRATION
# ==============================================================================
echo "========================================================"
echo " Final Steps: Merging & Integration"
echo "========================================================"

FINAL_OUTPUT="$WAGON_DIR/AnalysisResults_merged.root"

if [ -s "$MERGE_LIST_FILE" ]; then
    echo "  Merging QA files..."
    hadd -f -j 16 "$FINAL_OUTPUT" @"$MERGE_LIST_FILE"
    
    if [ $? -eq 0 ]; then
        echo "  [OK] Merged QA saved to: AnalysisResults_merged.root"
        xargs rm < "$MERGE_LIST_FILE"
        rm "$MERGE_LIST_FILE"
    else
        echo "  [ERROR] hadd merge failed."
        exit 1
    fi
else
    echo "  [WARNING] No QA results generated to merge."
fi

# Generate the input_data_storage.txt for the DerivedDataHY consumers
echo "  Generating O2 input list for downstream consumers..."
> "$WAGON_DIR/input_data_storage.txt"
for AOD in "$WAGON_DIR"/AO2Ds/AO2D_*.root; do
    # Ensure we actually found files and not an unexpanded glob
    if [ -f "$AOD" ]; then
        echo "file:${AOD}" >> "$WAGON_DIR/input_data_storage.txt"
    fi
done

echo ""
echo "========================================================"
echo " SUCCESS! Wagon fully synthesized locally."
echo " Output Dir: $WAGON_DIR"
echo " You can now run the consumer on this wagon using runDerivedDataConsumer_HY.sh"
echo "========================================================"