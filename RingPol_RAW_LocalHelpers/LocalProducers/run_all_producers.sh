#!/usr/bin/env bash
###############################################################################
# Script name:
#   run_all_producers.sh
#
# Purpose:
#   Orchestrator script that sequentially runs multiple O2 DPL producer 
#   configurations over a specified dataset. It handles dataset sampling, 
#   progress tracking, and error reporting. 
#
# Hardware Strategy & NUMA Isolation:
#   Targeted for dual-socket AMD EPYC systems (e.g., jarvis15 with 2x 9634).
#   To prevent catastrophic OS context-switching and cross-socket memory 
#   latency (Infinity Fabric penalty), this coordinator strictly enforces 
#   NUMA isolation. It wraps every worker process in:
#       numactl --cpunodebind=0 --preferred=0
#   This locks the ~56-process pipeline and its 128 GB shared memory pool 
#   entirely onto Node 0, maximizing CPU efficiency and utilizing the Linux 
#   Page Cache safely without triggering OOM kills.
#   CRITICAL: Do not run multiple instances of this script concurrently!
#
# Usage:
#   ./run_all_producers.sh [--do-log] [PERCENTAGE] [OUTPUT_DIR] [INPUT_LIST] [CONFIGS_DIR]
#
# Examples:
#   1. Run everything with default paths (100% of data):
#      ./run_all_producers.sh
#
#   2. Run a quick 10% test with terminal logging enabled:
#      ./run_all_producers.sh --do-log 10
#
#   3. Run 100% of data with custom paths, no logging:
#      ./run_all_producers.sh 100 /data/output /data/input.txt /data/configs
###############################################################################

# ==============================================================================
# DEFAULTS
# ==============================================================================
DEFAULT_DATASET_DIR="/home/users/cicerodm/RingPol/LHC25ae_pass2_local"
DEFAULT_INPUT_LIST="${DEFAULT_DATASET_DIR}/input_data_storage.txt"
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/producer_configs"

# --- ARGUMENT PARSING ---
DO_LOG_FLAG=""
POS_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    --do-log)
      DO_LOG_FLAG="--do-log"
      shift
      ;;
    *)
      POS_ARGS+=("$1")
      shift
      ;;
  esac
done

PERCENTAGE="${POS_ARGS[0]:-100}"
OUTPUT_DATASET_DIR="${POS_ARGS[1]:-$DEFAULT_DATASET_DIR}"
INPUT_LIST="${POS_ARGS[2]:-$DEFAULT_INPUT_LIST}"
PRODUCER_CONFIGS_DIR="${POS_ARGS[3]:-$DEFAULT_CONFIGS_DIR}"

# ==============================================================================
# PATHS
# ==============================================================================
REPO_DIR="/home/users/cicerodm/PhD_codes"
FRAMEWORK_DIR="${REPO_DIR}/RingPol_RAW_LocalHelpers/LocalProducers"
RUN_SCRIPT="${FRAMEWORK_DIR}/runProducerAsHY.sh"
AOD_WRITER="${REPO_DIR}/RingPol_RAW_LocalHelpers/Local_framework/aod-writer.json"

# ==============================================================================
# GLOBAL TIMERS
# ==============================================================================
SCRIPT_START_TIME=$(date +%s)
TOTAL_RUNTIME=0
PROCESSED=0

# ==============================================================================
# SIGNAL HANDLING
# ==============================================================================
FAILURES=()

handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  echo "    Stopping. Partial failure report:"
  echo ""

  if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "  No failures recorded before interrupt."
  else
    print_failure_table
  fi

  echo "========================================================"
  exit 130
}
trap handle_interrupt INT TERM

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
print_failure_table() {
  printf "  %-30s  %-40s  %-20s\n" "TYPE" "CONFIG" "FAILED STEP"
  printf "  %-30s  %-40s  %-20s\n" \
    "------------------------------" \
    "----------------------------------------" \
    "--------------------"

  for ENTRY in "${FAILURES[@]}"; do
    TYPE_COL=$( echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
    CONF_COL=$( echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
    STEP_COL=$( echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$3); print $3}')

    printf "  %-30s  %-40s  %-20s\n" "$TYPE_COL" "$CONF_COL" "$STEP_COL"
  done
}

format_time() {
  local T=$1
  printf "%02dh:%02dm:%02ds" $((T/3600)) $((T%3600/60)) $((T%60))
}

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================
for REQUIRED in "$RUN_SCRIPT" "$AOD_WRITER" "$INPUT_LIST"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: required file not found: ${REQUIRED}"
    exit 1
  fi
done

if [ ! -d "$PRODUCER_CONFIGS_DIR" ]; then
  echo "Error: PRODUCER_CONFIGS_DIR not found: ${PRODUCER_CONFIGS_DIR}"
  exit 1
fi

if [ ! -d "$OUTPUT_DATASET_DIR" ]; then
  echo "Warning: OUTPUT_DATASET_DIR does not exist, creating it."
  mkdir -p "$OUTPUT_DATASET_DIR" || exit 1
fi

# Collect configs
shopt -s nullglob
CONFIG_FILES=("${PRODUCER_CONFIGS_DIR}"/dpl-config-*.json)
shopt -u nullglob

TOTAL_CONFIGS=${#CONFIG_FILES[@]}

if [ $TOTAL_CONFIGS -eq 0 ]; then
  echo "Error: no dpl-config-*.json files found."
  exit 1
fi

# ==============================================================================
# INPUT SAMPLING (PERCENTAGE CONTROL)
# ==============================================================================

if ! awk "BEGIN {exit !($PERCENTAGE+0 == $PERCENTAGE)}"; then
  echo "Error: PERCENTAGE must be numeric (got '$PERCENTAGE')"
  exit 1
fi

if (( $(awk "BEGIN {print ($PERCENTAGE <= 0 || $PERCENTAGE > 100)}") )); then
  echo "Error: PERCENTAGE must be in (0, 100]"
  exit 1
fi

TOTAL_LINES=$(grep -cvE '^\s*$' "$INPUT_LIST")

if [ "$TOTAL_LINES" -eq 0 ]; then
  echo "Error: INPUT_LIST is empty."
  exit 1
fi

# Compute number of lines to keep (ceil)
LINES_TO_KEEP=$(awk -v total="$TOTAL_LINES" -v perc="$PERCENTAGE" \
  'BEGIN {printf "%d", (total * perc / 100.0) + 0.999}')

# Ensure at least 1 line
if [ "$LINES_TO_KEEP" -lt 1 ]; then
  LINES_TO_KEEP=1
fi

# Create temporary input list if not 100%
if (( $(awk "BEGIN {print ($PERCENTAGE < 100)}") )); then
  TEMP_INPUT_LIST=$(mktemp /tmp/input_subset_XXXX.txt)

  echo "  Sampling input list: ${PERCENTAGE}% (${LINES_TO_KEEP}/${TOTAL_LINES} files)"

  grep -vE '^\s*$' "$INPUT_LIST" | head -n "$LINES_TO_KEEP" > "$TEMP_INPUT_LIST"

  INPUT_LIST="$TEMP_INPUT_LIST"

  # Ensure cleanup
  trap 'rm -f "$TEMP_INPUT_LIST"' EXIT
else
  echo "  Using full dataset (100%)"
fi

echo ""

# ==============================================================================
# HEADER
# ==============================================================================
echo "========================================================"
echo "  run_all_producers.sh"
echo "  Configs found : ${TOTAL_CONFIGS}"
echo "  Input List    : ${INPUT_LIST}"
echo "  Data fraction : ${PERCENTAGE}%"
echo "  Used files    : ${LINES_TO_KEEP}"
if [ -n "$DO_LOG_FLAG" ]; then
  echo "  Logging       : ENABLED (--do-log)"
else
  echo "  Logging       : DISABLED"
fi
echo "========================================================"
echo ""

# ==============================================================================
# MAIN LOOP
# ==============================================================================
for CONFIG_FILE in "${CONFIG_FILES[@]}"; do
  CONFIG_BASENAME=$(basename "$CONFIG_FILE" .json)

  JOB_START=$(date +%s)

  echo "--------------------------------------------------------"
  echo "  Running Producer: ${CONFIG_BASENAME}"
  echo "--------------------------------------------------------"

  # Notice the unquoted $DO_LOG_FLAG here. If it's empty, bash just ignores it.
  numactl --cpunodebind=0 --preferred=0 "$RUN_SCRIPT" $DO_LOG_FLAG "$OUTPUT_DATASET_DIR" "$INPUT_LIST" "$CONFIG_FILE" "$AOD_WRITER"
  EXIT_CODE=$?

  JOB_END=$(date +%s)
  JOB_TIME=$((JOB_END - JOB_START))

  TOTAL_RUNTIME=$((TOTAL_RUNTIME + JOB_TIME))
  PROCESSED=$((PROCESSED + 1))

  # Progress %
  PERCENT=$(awk "BEGIN {printf \"%.2f\", (${PROCESSED}/${TOTAL_CONFIGS})*100}")

  # Average time per job
  AVG_TIME=$((TOTAL_RUNTIME / PROCESSED))

  # Remaining jobs
  REMAINING=$((TOTAL_CONFIGS - PROCESSED))

  # ETA
  ETA=$((AVG_TIME * REMAINING))

  if [ $EXIT_CODE -ne 0 ]; then
    echo "  -> FAILED: ${CONFIG_BASENAME}"
    FAILURES+=("TableProducer | ${CONFIG_BASENAME} | runProducerAsHY")
  else
    echo "  -> OK: ${CONFIG_BASENAME}"
  fi

  echo "  Time (this job) : $(format_time $JOB_TIME)"
  echo "  Total elapsed   : $(format_time $TOTAL_RUNTIME)"
  echo "  Progress        : ${PROCESSED}/${TOTAL_CONFIGS} (${PERCENT}%)"
  echo "  ETA remaining   : $(format_time $ETA)"
  echo ""
done

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
TOTAL_ELAPSED=$(( $(date +%s) - SCRIPT_START_TIME ))

echo "========================================================"
echo "  TOTAL ELAPSED TIME: $(format_time $TOTAL_ELAPSED)"
echo ""

if [ ${#FAILURES[@]} -eq 0 ]; then
  echo "  All producer configurations completed successfully."
else
  echo "  FAILURES (${#FAILURES[@]} total)"
  echo ""
  print_failure_table
fi

echo "========================================================"