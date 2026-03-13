#!/usr/bin/env bash

###############################################################################
# Script name:
#   runDerivedDataConsumer_HY.sh
#
# Purpose:
#   Runs the ALICE O2 derived-data consumer task on AO2D files downloaded
#   from Hyperloop. Unlike the raw-data workflow, there is no intermediate
#   "derived AOD producer" step: the files in AO2Ds/ are already the derived
#   input and are passed directly to the consumer task.
#
#   Input  : <WORK_DIR>/AO2Ds/AO2D_*.root
#   Output : <WORK_DIR>/results_consumer/ConsumerResults_<CONS_SUFFIX>.root
#   Log    : <WORK_DIR>/results_consumer/logs/consumer_run_<CONS_SUFFIX>.log
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
SHM_SIZE="32000000000"      # Shared memory: 32 GB (increased for 3D histograms)
MEM_RATE_LIMIT="2000000000" # AOD memory rate limit: 2 GB/s
THREADS=32                  # Pipeline threads for the derived consumer task

# ==============================================================================
# 1. ARGUMENT PARSING AND VALIDATION
# ==============================================================================
WORK_DIR="$1"
JSON_CONSUMER_CONFIG="$2"

if [ -z "$WORK_DIR" ] || [ -z "$JSON_CONSUMER_CONFIG" ]; then
  echo "Usage: $0 <WORK_DIR> <consumer_config.json>"
  echo ""
  echo "  WORK_DIR          : per-wagon working directory containing AO2Ds/"
  echo "  consumer_config   : dpl-config JSON for the consumer task"
  exit 1
fi

# Resolve to absolute paths so cd later does not break them
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
# The output file and log are named after the consumer config suffix only.
# For HY-derived data there is no data-producer suffix (unlike the raw
# workflow where a second JSON tracked the aod-writer configuration).
#
# Example:
#   dpl-config-DerivedConsumer-JustLambda.json  ->  CONS_SUFFIX = JustLambda
#   Output: ConsumerResults_JustLambda.root
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

# Remove any leftover staging directory from a previous run
rm -rf "$CONSUMER_TEMP"
mkdir -p "$CONSUMER_TEMP"

echo "========================================================"
echo "  HY Derived Data Consumer"
echo "  WORK_DIR         : ${WORK_DIR}"
echo "  Consumer suffix  : ${CONS_SUFFIX}"
echo "  Input            : ${AOD_DIR}/AO2D_*.root"
echo "  Output           : ${CONSUMER_RESULTS}/ConsumerResults_${CONS_SUFFIX}.root"
echo "========================================================"

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
# 5. BUILD INPUT FILE LIST
#
# Finds all AO2D_*.root files in AO2Ds/, sorts them numerically by index,
# and prepends "file:" so O2 treats them as local filesystem paths.
# ==============================================================================
DERIVED_LIST="${CONSUMER_TEMP}/derived_inputs.txt"

find "$AOD_DIR" -maxdepth 1 -name "AO2D_*.root" \
  | sort -t_ -k2 -V \
  | awk '{print "file:"$0}' \
  > "$DERIVED_LIST"

TOTAL_FILES=$(wc -l < "$DERIVED_LIST")

if [ "$TOTAL_FILES" -eq 0 ]; then
  echo "Error: No AO2D_*.root files found in: ${AOD_DIR}"
  exit 1
fi

echo "  Found ${TOTAL_FILES} input files. Processing all at once."
echo "--------------------------------------------------------"

# ==============================================================================
# 6. RUN THE CONSUMER PIPELINE
# ==============================================================================
MAIN_LOG="${CONSUMER_LOGS}/consumer_run_${CONS_SUFFIX}.log"
FINAL_OUTPUT="${CONSUMER_RESULTS}/ConsumerResults_${CONS_SUFFIX}.root"

# Run inside the staging directory so O2 scratch files stay contained
cd "$CONSUMER_TEMP" || exit 1

echo "  [O2] Launching consumer pipeline..."

time \
o2-analysis-lf-lambdajetpolarizationionsderived \
    -b \
    --configuration "json://${JSON_CONSUMER_CONFIG}" \
    --pipeline "lambdajetpolarizationionsderived:${THREADS}" \
    --aod-file "@${DERIVED_LIST}" \
    --aod-memory-rate-limit "$MEM_RATE_LIMIT" \
    --shm-segment-size "$SHM_SIZE" \
    > "$MAIN_LOG" 2>&1

EXIT_CODE=$?

# ==============================================================================
# 7. STAGE-OUT
# ==============================================================================
cd "$WORK_DIR" || exit 1

if [ $EXIT_CODE -eq 0 ]; then
  echo "  [OK] Consumer finished successfully."

  if [ -f "${CONSUMER_TEMP}/AnalysisResults.root" ]; then
    mv "${CONSUMER_TEMP}/AnalysisResults.root" "$FINAL_OUTPUT"
    # Save a copy of the config used, for reproducibility
    cp "$JSON_CONSUMER_CONFIG" "${CONSUMER_CONFIGS}/consumer-config_${CONS_SUFFIX}.json"
    echo ""
    echo "SUCCESS."
    echo "  Output : ${FINAL_OUTPUT}"
    echo "  Log    : ${MAIN_LOG}"
    echo "  Config : ${CONSUMER_CONFIGS}/consumer-config_${CONS_SUFFIX}.json"
  else
    echo "  Error: AnalysisResults.root not found in staging area."
    echo "  Check log: ${MAIN_LOG}"
    exit 1
  fi

else
  echo "  [FAIL] Consumer exited with code ${EXIT_CODE}."
  echo "  Check log: ${MAIN_LOG}"
  exit $EXIT_CODE
fi
