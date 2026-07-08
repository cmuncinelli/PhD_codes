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
#     1. Run runDerivedDataConsumer_HY.sh -> ConsumerResults_<SUFFIX>.root
#     2. Run extractDeltaErrors.cxx (ROOT) -> delta/error plots
#     3. Run signalExtractionRing.cxx (ROOT) -> results_SigExtract/
#     4. Run makeCumulativeDCAdauProfile.cxx (ROOT)-> CumulativeProfiles_<SUFFIX>.root
#     5. Run auxiliaryPlots.cxx (ROOT) -> AuxiliaryPlots.root (Cross-config)
#
#   Steps 2, 3 and 4 are skipped for a given pair if step 1 fails or if the 
#   ConsumerResults file is missing.
#
# Usage:
#   ./run_all_wagons.sh [OPTIONS] [REGISTRY] [CONSUMER_CONFIGS_DIR]
#
# Options:
#   -h, --help                 Show this help message.
#   -p, --post-process-only    Skip Step 1 (consumer) and run only the ROOT
#                              macros (Steps 2, 3, 4) on existing output.
#   -s, --skip-sig-extract     Skip Step 3 (signal extraction) during execution.
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
CUMUL_DCA_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/makeCumulativeDCAdauProfile.cxx"
AUXILIARY_PLOTS_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliaryPlots.cxx"

# Absolute path to the consumer launcher script.
FRAMEWORK_DIR="${REPO_DIR}/RingPol_RAW_LocalHelpers/DerivedDataHY"
CONSUMER_SCRIPT="${FRAMEWORK_DIR}/runDerivedDataConsumer_HY.sh"
DEFAULT_REGISTRY="${FRAMEWORK_DIR}/train_registry.conf"
# Directory containing dpl-config-DerivedConsumer-*.json files.
# Can be overridden by passing a path as $2.
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/consumer_configs"

# Directory containing the MC reference ConsumerResults files for Auxiliary Plots
# Leave empty ("") if you do not want to overlay MC.
MC_REF_DIR="/home/users/cicerodm/RingPol/LHC25h3c/ITSandTPC_min3ITS/results_consumer" # TODO: add a remark in the README about this hardcoded path!!!
                                                                                      # Same for the DEFAULT_CONFIGS_DIR path could probably come in hand

# Absolute path to the Toy Model ROOT file for Auxiliary Plots (chose the Toy Model representative whose configurations are 
# closest to what an actual data selection would do in data processing, as an attempt to keep everything consistent)
# Leave empty ("") if you do not want to overlay the Toy Model.
TOY_MODEL_PATH="/home/users/cicerodm/RingPol/HelicityToyModel/9_RealisticAlice/alice_an_std/helicity_alice_an_std.root" # TODO: Document this path in the README too!
# For testing with a model configuration that is already on disk:
# TOY_MODEL_PATH="/home/users/cicerodm/RingPol/HelicityToyModel/9_RealisticAlice/alice_std/helicity_alice_std.root"

# Set executable paths
EXTRACT_DELTA_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/extractDeltaErrors.exe"
SIGNAL_EXTRACT_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/signalExtractionRing.exe"
CUMUL_DCA_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/makeCumulativeDCAdauProfile.exe"
AUXILIARY_PLOTS_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliaryPlots.exe"

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================
POST_PROCESS_ONLY=0
SKIP_SIG_EXTRACT=0
REGISTRY_ARG=""
CONFIGS_DIR_ARG=""

print_help() {
    awk '/^# Usage:/,/^#+$/' "$0" | head -n -1 | sed 's/^#//'
    exit 0
}

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            print_help
            ;;
        --post-process-only|-p)
            POST_PROCESS_ONLY=1
            shift
            ;;
        --skip-sig-extract|-s)
            SKIP_SIG_EXTRACT=1
            shift
            ;;
        -*)
            echo "Unknown option: $1. Use --help for usage." >&2
            exit 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Use the first argument if provided, otherwise fall back to the legacy default
# (this is meant essentially to be able to choose a .conf file for the locally
#  produced derived data instead of the HY produced datasets)
REGISTRY="${1:-$DEFAULT_REGISTRY}"
REGISTRY="${POSITIONAL_ARGS[0]:-$DEFAULT_REGISTRY}"
CONSUMER_CONFIGS_DIR="${POSITIONAL_ARGS[1]:-$DEFAULT_CONFIGS_DIR}"

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

cleanup() {
  # Removing the compiled binaries after usage:
  rm -f "$EXTRACT_DELTA_EXE" "$SIGNAL_EXTRACT_EXE" "$CUMUL_DCA_EXE" "$AUXILIARY_PLOTS_EXE"
}

trap cleanup EXIT
trap handle_interrupt INT TERM

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================
for REQUIRED in "$REGISTRY" "$CONSUMER_SCRIPT" \
                "$EXTRACT_DELTA_MACRO" "$SIGNAL_EXTRACT_MACRO" "$CUMUL_DCA_MACRO" "$AUXILIARY_PLOTS_MACRO"; do
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
# Compile Analysis Executables AoT
# ==============================================================================
echo "Compiling ROOT macros into native executables..."

COMPILE_WARN_FLAGS="-Wall -Wextra -Wpedantic -Wshadow" # For bug catching (TODO: check some unused variables that were revealed! May be buggy code!)
# COMPILE_WARN_FLAGS="" # For quick compiling without warnings (do NOT use me without care!)
ROOT_AWARE_FLAGS="-DNDEBUG $(root-config --cflags)" # Strips out all assert() checks. ROOT headers respect this flag!
OPTIMIZATION_FLAGS="-O3 -march=native -fno-math-errno -flto -fno-trapping-math" # Using the fno-math-errno to disable the errno verification of passing negative numbers to sqrts and the such.
                                                                                # My code should be safe, so this will never be used!
                                                                                # Also introduced -fno-trapping-math as there is no explicit usage of NaNs
                                                                                # nor any expected unsafe divisions by zero appearing anywhere in the code!
ROOT_LIBS="$(root-config --glibs)"

# Compile Step (Stops script if compilation fails)
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$EXTRACT_DELTA_EXE" "$EXTRACT_DELTA_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$SIGNAL_EXTRACT_EXE" "$SIGNAL_EXTRACT_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$CUMUL_DCA_EXE" "$CUMUL_DCA_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$AUXILIARY_PLOTS_EXE" "$AUXILIARY_PLOTS_MACRO" $ROOT_LIBS || exit 1

echo "Compilation successful!"
echo ""

# ==============================================================================
# SUMMARY HEADER
# ==============================================================================
echo "========================================================"
echo "  run_all_wagons.sh"
echo "  Wagons       : ${#WAGON_LINES[@]}"
echo "  Configs      : ${#CONFIG_FILES[@]}"
echo "  Configs Dir  : ${CONSUMER_CONFIGS_DIR}"
echo "  Toy Model    : ${TOY_MODEL_PATH:-None}"
echo "  Mode         : $( [ $POST_PROCESS_ONLY -eq 1 ] && echo 'POST-PROCESS ONLY' || echo 'FULL CHAIN' )$( [ $SKIP_SIG_EXTRACT -eq 1 ] && echo ' [SKIP SIG EXTRACT]' )"
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
    
    # Logging Setup
    LOG_DIR="${WORK_DIR}/results_consumer/logs"
    mkdir -p "${LOG_DIR}"
      # Creating smaller log folders to keep everything tidy (it was getting really messy!)
    mkdir -p "${LOG_DIR}/wrappers/"
    mkdir -p "${LOG_DIR}/batches/"
    mkdir -p "${LOG_DIR}/deltaErr/"
    mkdir -p "${LOG_DIR}/sigExtract/"
    mkdir -p "${LOG_DIR}/cumulHist/"
    WRAPPER_LOG="${LOG_DIR}/wrappers/wrapper_${CONS_SUFFIX}.log"
    DELTA_LOG="${LOG_DIR}/deltaErr/extractDeltaErr_${CONS_SUFFIX}.log"
    SIG_LOG="${LOG_DIR}/sigExtract/sigExtract_${CONS_SUFFIX}.log"
    CUMUL_LOG="${LOG_DIR}/cumulHist/cumulDCA_${CONS_SUFFIX}.log"

    # ------------------------------------------------------------------
    # Step 1: consumer
    # ------------------------------------------------------------------
    if [ $POST_PROCESS_ONLY -eq 0 ]; then
      echo -n "  [1/5] consumer        : ${CONS_SUFFIX}"
      # Consumer output goes to a per-config log file so failures are inspectable.
      if [ -d /sys/devices/system/node/node1 ]; then
        # Binding consumer to the NUMA node1 (just convenience: producers are running in node 0 on jarvis15 right now)
          numactl --cpunodebind=1 --preferred=1  "$CONSUMER_SCRIPT" "$WORK_DIR" "$CONFIG_FILE" > "$WRAPPER_LOG" 2>&1
      else # If does not have more than one node, just revert to usual behavior!
          "$CONSUMER_SCRIPT" "$WORK_DIR" "$CONFIG_FILE" > "$WRAPPER_LOG" 2>&1
      fi
      CONSUMER_EXIT=$?

      if [ $CONSUMER_EXIT -ne 0 ] || [ ! -f "$CONSUMER_RESULT" ]; then
        echo "  -> FAILED  (log: ${WRAPPER_LOG})"
        FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | consumer")
        continue
      fi
      echo "  -> OK"
    else
      # Post-process only mode: verify file exists
      echo -n "  [1/5] consumer        : ${CONS_SUFFIX} (SKIPPED)"
      if [ ! -f "$CONSUMER_RESULT" ]; then
        echo "  -> FAILED  (File missing)"
        FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | missing_consumer_result")
        continue
      fi
      echo "  -> OK (Found file)"
    fi

    # ------------------------------------------------------------------
    # Step 2: extractDeltaErrors
    # ------------------------------------------------------------------
    echo -n "  [2/5] extractDeltaErr : ${CONS_SUFFIX}"
    "$EXTRACT_DELTA_EXE" "${CONSUMER_RESULT}" > "$DELTA_LOG" 2>&1
    DELTA_EXIT=$?

    if [ $DELTA_EXIT -ne 0 ]; then
      echo "  -> FAILED  (log: ${DELTA_LOG})"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | extractDeltaErrors")
    else
      echo "  -> OK"
    fi

    # ------------------------------------------------------------------
    # Step 3: signalExtractionRing
    # ------------------------------------------------------------------
    if [ $SKIP_SIG_EXTRACT -eq 0 ]; then
      echo -n "  [3/5] sigExtract      : ${CONS_SUFFIX}"
      "$SIGNAL_EXTRACT_EXE" "${CONSUMER_RESULT}" "${SIGNAL_EXTRACT_DIR}/" > "$SIG_LOG" 2>&1
      SIG_EXIT=$?

      if [ $SIG_EXIT -ne 0 ]; then
        echo "  -> FAILED  (log: ${SIG_LOG})"
        FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | signalExtractionRing")
      else
        echo "  -> OK"
      fi
    else
      echo "  [3/5] sigExtract      : ${CONS_SUFFIX} (SKIPPED)"
    fi

    # ------------------------------------------------------------------
    # Step 4: makeCumulativeDCAdauProfile
    # ------------------------------------------------------------------
    echo -n "  [4/5] cumulDCA        : ${CONS_SUFFIX}"
    "$CUMUL_DCA_EXE" "${CONSUMER_RESULT}" > "$CUMUL_LOG" 2>&1
    CUMUL_EXIT=$?

    if [ $CUMUL_EXIT -ne 0 ]; then
      echo "  -> FAILED  (log: ${CUMUL_LOG})"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | makeCumulDCA")
    else
      echo "  -> OK"
    fi

    echo ""

  done  # configs

  # ------------------------------------------------------------------
  # Step 5: auxiliaryPlots (Cross-configuration aggregation)
  # ------------------------------------------------------------------
  # We run this ONCE per wagon, passing the wagon's base working dir
  AUX_LOG="${WORK_DIR}/results_consumer/logs/auxPlots.log" # Saves a single log, to the root of the logs folder, because there will be only one single log for a given wagon
  echo -n "  [5/5] auxiliaryPlots  : (Cross-config summary)"
  
  # Forwarding the consumer results directory, the MC reference, and now the Toy Model path:
  "$AUXILIARY_PLOTS_EXE" "${WORK_DIR}/results_consumer" "${MC_REF_DIR}" "${TOY_MODEL_PATH}" > "$AUX_LOG" 2>&1
  AUX_EXIT=$?

  if [ $AUX_EXIT -ne 0 ]; then
    echo "  -> FAILED  (log: ${AUX_LOG})"
    FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ALL_CONFIGS | auxiliaryPlots")
  else
    echo "  -> OK"
  fi

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
  echo "  (Check the 'results_consumer/logs' folder for the respective wagon)"
fi
echo "========================================================"