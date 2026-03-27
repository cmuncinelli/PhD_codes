#!/usr/bin/env bash

###############################################################################
# Script name:
#   run_all_wagons.sh
#
# Purpose:
#   Optional end-to-end wrapper that runs the full analysis chain for every
#   wagon registered in a train registry file, against every consumer config
#   found in a given directory.
#
#   For each (wagon, config) pair it will:
#
#     1. Run runDerivedDataConsumer_HY.sh   -> ConsumerResults_<SUFFIX>.root
#     2. Run extractDeltaErrors.cxx (ROOT)  -> delta/error plots
#     3. Run signalExtractionRing.cxx (ROOT)-> results_SigExtract/
#
#   Steps 2 and 3 are skipped for a given pair if step 1 fails.
#   All failures are collected and printed as a summary table at the end.
#
# Usage:
#   ./run_all_wagons.sh [REGISTRY] [CONSUMER_CONFIGS_DIR]
#
# Arguments:
#   REGISTRY (optional):
#     Path to a train_registry.conf file.
#     Defaults to:
#       ${FRAMEWORK_DIR}/train_registry.conf
#
#   CONSUMER_CONFIGS_DIR (optional):
#     Directory containing dpl-config-DerivedConsumer-*.json files.
#     Defaults to:
#       /home/users/cicerodm/RingPol/consumer_configs
#
# Examples:
#   ./run_all_wagons.sh
#   ./run_all_wagons.sh /path/to/my_registry.conf
#   ./run_all_wagons.sh /path/to/my_registry.conf /path/to/my_configs
#
###############################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Root of the PhD_codes git repository. All macro paths are derived from here.
REPO_DIR="/home/users/cicerodm/PhD_codes"

# Absolute paths to the ROOT macros used in the post-processing steps.
EXTRACT_DELTA_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/extractDeltaErrors.cxx"
SIGNAL_EXTRACT_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/signalExtractionRing.cxx"

# Absolute path to the consumer launcher script.
FRAMEWORK_DIR="${REPO_DIR}/RingPol_RAW_LocalHelpers/DerivedDataHY"
CONSUMER_SCRIPT="${FRAMEWORK_DIR}/runDerivedDataConsumer_HY.sh"
DEFAULT_REGISTRY="${FRAMEWORK_DIR}/train_registry.conf"

# Use the first argument if provided, otherwise fall back to the legacy default
# (this is meant essentially to be able to choose a .conf file for the locally
#  produced derived data instead of the HY produced datasets)
REGISTRY="${1:-$DEFAULT_REGISTRY}"

# Directory containing dpl-config-DerivedConsumer-*.json files.
# Can be overridden by passing a path as $2.
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/consumer_configs"
CONSUMER_CONFIGS_DIR="${2:-$DEFAULT_CONFIGS_DIR}"

# ==============================================================================
# SIGNAL HANDLING
# ==============================================================================
# Note on trap interaction with child scripts:
#   Ctrl+C sends SIGINT to the entire foreground process group, so both this
#   script and whichever child (.sh or ROOT macro) is currently running will
#   each receive the signal and fire their own traps independently. There is
#   no need to forward signals manually for the normal interactive use case.
#   Killing this script by PID from another terminal would only hit this
#   process; in that case child traps do not fire automatically.

handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  echo "    Stopping. Partial failure report:"
  echo ""
  if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "  No failures recorded before interrupt."
  else
    printf "  %-45s  %-40s  %-20s
" "WAGON" "CONFIG" "FAILED STEP"
    printf "  %-45s  %-40s  %-20s
"       "---------------------------------------------"       "----------------------------------------"       "--------------------"
    for ENTRY in "${FAILURES[@]}"; do
      WAGON_COL=$(echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
      CONF_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
      STEP_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$3); print $3}')
      printf "  %-45s  %-40s  %-20s
" "$WAGON_COL" "$CONF_COL" "$STEP_COL"
    done
  fi
  echo "========================================================"
  exit 130
}

trap handle_interrupt INT TERM

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================
for REQUIRED in "$REGISTRY" "$CONSUMER_SCRIPT" \
                "$EXTRACT_DELTA_MACRO" "$SIGNAL_EXTRACT_MACRO"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: required file not found: ${REQUIRED}"
    exit 1
  fi
done

if [ ! -d "$CONSUMER_CONFIGS_DIR" ]; then
  echo "Error: CONSUMER_CONFIGS_DIR not found: ${CONSUMER_CONFIGS_DIR}"
  exit 1
fi

# Collect consumer config files
shopt -s nullglob
CONFIG_FILES=("${CONSUMER_CONFIGS_DIR}"/dpl-config-DerivedConsumer-*.json)
shopt -u nullglob

if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
  echo "Error: no dpl-config-DerivedConsumer-*.json files found in: ${CONSUMER_CONFIGS_DIR}"
  exit 1
fi

# Collect wagon entries from registry (skip comments and blanks)
WAGON_LINES=()
while IFS= read -r LINE || [ -n "$LINE" ]; do
  [[ "$LINE" =~ ^\s*# ]] && continue
  [[ -z "${LINE// }"   ]] && continue
  WAGON_LINES+=("$LINE")
done < "$REGISTRY"

if [ ${#WAGON_LINES[@]} -eq 0 ]; then
  echo "Error: no wagons found in registry: ${REGISTRY}"
  exit 1
fi

# ==============================================================================
# SUMMARY HEADER
# ==============================================================================
echo "========================================================"
echo "  run_all_wagons.sh"
echo "  Wagons    : ${#WAGON_LINES[@]}"
echo "  Configs   : ${#CONFIG_FILES[@]}"
echo "  Configs   : ${CONSUMER_CONFIGS_DIR}"
echo "========================================================"
echo ""

# Failure log: each entry is "WAGON_SHORTNAME | CONFIG_SUFFIX | STAGE"
FAILURES=()

# ==============================================================================
# MAIN LOOP
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

    # Derive the output suffix the consumer will use, e.g. "JustLambda"
    CONFIG_BASENAME=$(basename "$CONFIG_FILE" .json)
    CONS_SUFFIX="${CONFIG_BASENAME#dpl-config-DerivedConsumer-}"

    CONSUMER_RESULT="${WORK_DIR}/results_consumer/ConsumerResults_${CONS_SUFFIX}.root"
    SIGNAL_EXTRACT_DIR="${WORK_DIR}/results_SigExtract"
    WRAPPER_LOG="${WORK_DIR}/results_consumer/logs/wrapper_${CONS_SUFFIX}.log"
    mkdir -p "${WORK_DIR}/results_consumer/logs"

    # ------------------------------------------------------------------
    # Step 1: run consumer
    # ------------------------------------------------------------------
    echo -n "  [1/3] consumer        : ${CONS_SUFFIX}"
    # Consumer output goes to a per-config log file so failures are inspectable.
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
    root -l -b -q \
      "${EXTRACT_DELTA_MACRO}(\"${CONSUMER_RESULT}\")" \
      > /dev/null 2>&1
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
    root -l -b -q \
      "${SIGNAL_EXTRACT_MACRO}(\"${CONSUMER_RESULT}\", \"${SIGNAL_EXTRACT_DIR}/\")" \
      > /dev/null 2>&1
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
  echo "  All steps completed successfully."
else
  echo "  FAILURES (${#FAILURES[@]} total)"
  echo ""
  # Column headers
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