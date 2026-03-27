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
# Hardware Target: Jarvis15 (AMD EPYC 9634 84-Core/336T) | 512GB RAM | 4TB SSD
# Strategy:
#   1. Stage LAN data sequentially to SSD (overcoming the 1 Gbps bottleneck).
#   2. Blast through the batch using O2 with scaled-up threads & 256GB SHM.
#   3. Save derived AO2Ds and QA directly into the HY-style folder structure.
#   4. Auto-generate input_data_storage.txt for immediate consumer use.
#
# Usage:
#   ./runProducerAsHY.sh <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG_JSON> <AOD_WRITER_JSON>
#
# Example:
#   ./runProducerAsHY.sh /home/users/cicerodm/RingPol/LHC25ae_pass2_local \
#                        /path/to/input_data_storage.txt \
#                        /home/users/cicerodm/RingPol/producer_configs/dpl-config-ITSandTPC.json \
#                        /home/users/cicerodm/PhD_codes/.../aod-writer.json
###############################################################################

# --- INPUT ARGUMENTS ---
OUTPUT_DATASET_DIR="$1"
INPUT_LIST="$2"
JSON_CONFIG_PATH="$3"
AOD_WRITER_JSON="$4"

if [[ -z "$OUTPUT_DATASET_DIR" ]] || [[ -z "$INPUT_LIST" ]] || [[ -z "$JSON_CONFIG_PATH" ]] || [[ -z "$AOD_WRITER_JSON" ]]; then
    echo "Usage: $0 <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG_JSON> <AOD_WRITER_JSON>"
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
rm -rf "$TEMP_BASE"
mkdir -p "$TEMP_BASE/batches"

# Copy the config for reproducibility
cp "$JSON_CONFIG_PATH" "$WAGON_DIR/results_producer/used_configs/"

# --- TUNING KNOBS (Supercharged for 336T / 512GB RAM) ---
FILES_PER_BATCH=50             # ~150GB per batch (safe for 4TB SSD even with overlap)
SHM_SIZE="200000000000"        # 200 GB (leave RAM headroom for system + cache)

# --- DETECT HARDWARE ---
TOTAL_CORES=$(lscpu | awk '/^Core\(s\) per socket:/ {cores=$4} /^Socket\(s\):/ {sockets=$2} END {print cores*sockets}')
TOTAL_THREADS=$(nproc)

# --- LIMIT: USE ONLY HALF OF CORES (IMPORTANT) ---
MAX_CORES_ALLOWED=$((TOTAL_CORES / 2))
MAX_THREADS_ALLOWED=$((MAX_CORES_ALLOWED * 2))   # assuming SMT (2 threads/core)

# --- PIPELINE PARALLELISM ---
PIPE_EVENTSEL=8
PIPE_PROPAGATION=24
PIPE_PIDTPC=8
PIPE_TOF=4
PIPE_STRANGE=6
PIPE_LAMBDA=16

# --- COMPUTE USED THREADS AUTOMATICALLY ---
USED_THREADS=$((PIPE_EVENTSEL + PIPE_PROPAGATION + PIPE_PIDTPC + PIPE_TOF + PIPE_STRANGE + PIPE_LAMBDA))

# --- I/O TUNING ---
READERS=6
MEM_RATE_LIMIT="3500000000" # 3.5 GB/s

echo "  Resource configuration:"
echo "    Total cores        : $TOTAL_CORES"
echo "    Total threads      : $TOTAL_THREADS"
echo "    Max usable cores   : $MAX_CORES_ALLOWED"
echo "    Max usable threads : $MAX_THREADS_ALLOWED"
echo "    Used threads       : $USED_THREADS"
echo "    Readers            : $READERS"

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
    while read -r LINE; do
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
    # B. EXECUTE PIPELINE (Scaled safely for shared system)
    # ---------------------------------------------------------
    echo "  [O2] Launching controlled parallel pipeline..."
    echo "       (Check logs at: results_producer/logs/batch_$BATCH_ID.log)"

    cd "$BATCH_WORK_DIR" || exit 1

    OPTION="-b --configuration json://$JSON_CONFIG_PATH --aod-writer-json ${AOD_WRITER_JSON}"

    # Use nice + ionice to prevent system lock-up
    time \
    nice -n 10 ionice -c2 -n7 \
    o2-analysis-event-selection-service ${OPTION} \
        --pipeline eventselection-run3:${PIPE_EVENTSEL} | \
    o2-analysis-multcenttable ${OPTION} | \
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
        --aod-file "@${LOCAL_INPUT_LIST}" \
        --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
        --shm-segment-size "$SHM_SIZE" \
        > "$BATCH_LOG" 2>&1

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
        echo "  Check log: $BATCH_LOG"
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