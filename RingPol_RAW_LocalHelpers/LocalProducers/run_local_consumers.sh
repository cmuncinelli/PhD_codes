#!/usr/bin/env bash
###############################################################################
# Script name:
#   run_local_consumers.sh
#
# Purpose:
#   Coordinator script for running the derived data consumer on LOCALLY
#   produced wagons. It bridges the gap between your local producer outputs
#   and the existing DerivedDataHY framework.
#
#   How it tricks the old framework:
#     1. It auto-discovers valid wagons inside your local dataset folder.
#     2. It generates a "fake" 'local_registry.conf' to maintain a paper trail.
#     3. It calls the EXACT SAME runDerivedDataConsumer_HY.sh script from the 
#        old framework, passing it the expected WORK_DIR. The old script is 
#        none the wiser that the data didn't come from Hyperloop!
#     4. Runs the ROOT post-processing macros (delta errors & sig extract).
#
# Usage:
#   ./run_local_consumers.sh [LOCAL_DATASET_DIR] [CONSUMER_CONFIGS_DIR]
#
# Example:
#   ./run_local_consumers.sh /home/users/cicerodm/RingPol/LHC25ae_pass2_local
###############################################################################

# ==============================================================================
# CONFIGURATION & PATHS
# ==============================================================================
# Root of the PhD_codes git repository.
REPO_DIR="/home/users/cicerodm/PhD_codes"

# We strictly point to the OLD framework for the heavy lifting scripts
HY_FRAMEWORK_DIR="${REPO_DIR}/RingPol_RAW_LocalHelpers/DerivedDataHY"
CONSUMER_SCRIPT="${HY_FRAMEWORK_DIR}/runDerivedDataConsumer_HY.sh"

# Post-processing ROOT macros
EXTRACT_DELTA_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/extractDeltaErrors.cxx"
SIGNAL_EXTRACT_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/signalExtractionRing.cxx"

# Default Arguments
DEFAULT_DATASET_DIR="/home/users/cicerodm/RingPol/LHC25ae_pass2_local"
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/consumer_configs"

LOCAL_DATASET_DIR="${1:-$DEFAULT_DATASET_DIR}"
CONSUMER_CONFIGS_DIR="${2:-$DEFAULT_CONFIGS_DIR}"

# ==============================================================================
# SIGNAL HANDLING
# ==============================================================================
handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  echo "    Stopping. Partial failure report:"
  echo ""
  if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "  No failures recorded before interrupt."
  else
    printf "  %-45s  %-40s  %-20s\n" "WAGON" "CONFIG" "FAILED STEP"
    printf "  %-45s  %-40s  %-20s\n" \
      "---------------------------------------------" \
      "----------------------------------------" \
      "--------------------"
    for ENTRY in "${FAILURES[@]}"; do
      WAGON_COL=$(echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
      CONF_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
      STEP_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$3); print $3}')
      printf "  %-45s  %-40s  %-20s\n" "$WAGON_COL" "$CONF_COL" "$STEP_COL"
    done
  fi
  echo "========================================================"
  exit 130
}
trap handle_interrupt INT TERM

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================
for REQUIRED in "$CONSUMER_SCRIPT" "$EXTRACT_DELTA_MACRO" "$SIGNAL_EXTRACT_MACRO"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: required legacy script not found: ${REQUIRED}"
    exit 1
  fi
done

if [ ! -d "$CONSUMER_CONFIGS_DIR" ]; then
  echo "Error: CONSUMER_CONFIGS_DIR not found: ${CONSUMER_CONFIGS_DIR}"
  exit 1
fi

shopt -s nullglob
CONFIG_FILES=("${CONSUMER_CONFIGS_DIR}"/dpl-config-DerivedConsumer-*.json)
shopt -u nullglob

if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
  echo "Error: no dpl-config-DerivedConsumer-*.json files found in: ${CONSUMER_CONFIGS_DIR}"
  exit 1
fi

# ==============================================================================
# AUTO-GENERATE LOCAL REGISTRY (TRICKING THE FRAMEWORK)
# ==============================================================================
# Instead of manually editing a .conf file, we scan the local dataset and 
# dynamically build a registry in memory/on-disk for this run.

DATASET_NAME=$(basename "$LOCAL_DATASET_DIR")
STORAGE_BASE=$(dirname "$LOCAL_DATASET_DIR")
LOCAL_REGISTRY="${LOCAL_DATASET_DIR}/local_registry_auto.conf"

echo "# Auto-generated local registry for testing" > "$LOCAL_REGISTRY"
echo "# TRAIN_RUN_ID  DATASET_NAME  WAGON_SHORTNAME  STORAGE_BASE" >> "$LOCAL_REGISTRY"

WAGON_LINES=()
# Find every folder that has an AO2Ds subfolder
for WAGON_DIR in "$LOCAL_DATASET_DIR"/*; do
  if [ -d "$WAGON_DIR/AO2Ds" ]; then
    WAGON_SHORTNAME=$(basename "$WAGON_DIR")
    # We assign a dummy train ID '999999' to satisfy the formatting
    LINE="999999  ${DATASET_NAME}  ${WAGON_SHORTNAME}  ${STORAGE_BASE}"
    echo "$LINE" >> "$LOCAL_REGISTRY"
    WAGON_LINES+=("$LINE")
  fi
done

if [ ${#WAGON_LINES[@]} -eq 0 ]; then
  echo "Error: no valid wagons (folders containing AO2Ds/) found in: ${LOCAL_DATASET_DIR}"
  exit 1
fi

# ==============================================================================
# SUMMARY HEADER
# ==============================================================================
echo "========================================================"
echo "  run_local_consumers.sh (Tricking HY Framework)"
echo "  Dataset   : ${DATASET_NAME}"
echo "  Wagons    : ${#WAGON_LINES[@]} (Auto-discovered)"
echo "  Configs   : ${#CONFIG_FILES[@]}"
echo "========================================================"
echo ""

FAILURES=()

# ==============================================================================
# MAIN LOOP (Exact equivalent of the old HY loop)
# ==============================================================================
for LINE in "${WAGON_LINES[@]}"; do
  
  DATASET_NAME=$(   echo "$LINE" | awk '{print $2}')
  WAGON_SHORTNAME=$(echo "$LINE" | awk '{print $3}')
  STORAGE_BASE=$(   echo "$LINE" | awk '{print $4}')
  WORK_DIR="${STORAGE_BASE}/${DATASET_NAME}/${WAGON_SHORTNAME}"

  echo "--------------------------------------------------------"
  echo "  Wagon: ${DATASET_NAME}/${WAGON_SHORTNAME}"
  echo "--------------------------------------------------------"

  for CONFIG_FILE in "${CONFIG_FILES[@]}"; do

    CONFIG_BASENAME=$(basename "$CONFIG_FILE" .json)
    CONS_SUFFIX="${CONFIG_BASENAME#dpl-config-DerivedConsumer-}"

    CONSUMER_RESULT="${WORK_DIR}/results_consumer/ConsumerResults_${CONS_SUFFIX}.root"
    SIGNAL_EXTRACT_DIR="${WORK_DIR}/results_SigExtract"
    WRAPPER_LOG="${WORK_DIR}/results_consumer/logs/wrapper_${CONS_SUFFIX}.log"
    mkdir -p "${WORK_DIR}/results_consumer/logs"

    # ------------------------------------------------------------------
    # Step 1: run consumer (Calls the old HY script!)
    # ------------------------------------------------------------------
    echo -n "  [1/3] consumer        : ${CONS_SUFFIX}"
    
    # This invokes the legacy script. It looks at $WORK_DIR/AO2Ds and processes them.
    "$CONSUMER_SCRIPT" "$WORK_DIR" "$CONFIG_FILE" > "$WRAPPER_LOG" 2>&1
    CONSUMER_EXIT=$?

    if [ $CONSUMER_EXIT -ne 0 ] || [ ! -f "$CONSUMER_RESULT" ]; then
      echo "  -> FAILED  (log: ${WRAPPER_LOG})"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | consumer")
      continue
    fi
    echo "  -> OK"

    # ------------------------------------------------------------------
    # Step 2: extractDeltaErrors
    # ------------------------------------------------------------------
    echo -n "  [2/3] extractDeltaErr : ${CONS_SUFFIX}"
    root -l -b -q "${EXTRACT_DELTA_MACRO}(\"${CONSUMER_RESULT}\")" > /dev/null 2>&1
    DELTA_EXIT=$?

    if [ $DELTA_EXIT -ne 0 ]; then
      echo "  -> FAILED"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | extractDeltaErrors")
    else
      echo "  -> OK"
    fi

    # ------------------------------------------------------------------
    # Step 3: signalExtractionRing
    # ------------------------------------------------------------------
    echo -n "  [3/3] sigExtract      : ${CONS_SUFFIX}"
    root -l -b -q "${SIGNAL_EXTRACT_MACRO}(\"${CONSUMER_RESULT}\", \"${SIGNAL_EXTRACT_DIR}/\")" > /dev/null 2>&1
    SIG_EXIT=$?

    if [ $SIG_EXIT -ne 0 ]; then
      echo "  -> FAILED"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | signalExtractionRing")
    else
      echo "  -> OK"
    fi

    echo ""

  done  # configs
  echo ""
done  # wagons

# ==============================================================================
# FAILURE SUMMARY TABLE
# ==============================================================================
echo "========================================================"
if [ ${#FAILURES[@]} -eq 0 ]; then
  echo "  All local wagons completed successfully!"
else
  echo "  FAILURES (${#FAILURES[@]} total)"
  echo ""
  printf "  %-45s  %-40s  %-20s\n" "WAGON" "CONFIG" "FAILED STEP"
  printf "  %-45s  %-40s  %-20s\n" \
    "---------------------------------------------" \
    "----------------------------------------" \
    "--------------------"
  for ENTRY in "${FAILURES[@]}"; do
    WAGON_COL=$(echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
    CONF_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
    STEP_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$3); print $3}')
    printf "  %-45s  %-40s  %-20s\n" "$WAGON_COL" "$CONF_COL" "$STEP_COL"
  done
fi
echo "========================================================"