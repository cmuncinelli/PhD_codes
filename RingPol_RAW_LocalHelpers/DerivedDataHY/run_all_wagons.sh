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
#     1. Run runDerivedDataConsumer_HY.sh -> results_consumer/ConsumerResults_<SUFFIX>.root
#     2. Run extractDeltaErrors.cxx (ROOT) -> results_DeltaErr/
#     3. Run signalExtractionRing.cxx (ROOT) -> results_SigExtract/
#     4. Run makeCumulativeDCAdauProfile.cxx (ROOT)-> results_CumulativePlots/
#     5. Run auxiliaryPerConfigPlots.cxx (ROOT) -> results_AuxPerConfig/
#     6. Run auxiliarySummaryPlots.cxx (ROOT) -> results_consumer/auxiliarySummaryPlots.root (Cross-config)
#     7. Run zvtxBitForensics.cxx (ROOT) -> results_consumer/zvtxBitForensics.root (Per-wagon, parallel)
#
#   Steps 2, 4, and 5 each write into their own dedicated folder, mirroring how step 3 already uses results_SigExtract/,
#   instead of writing back into results_consumer/ itself.
#
#   Steps 2, 3, 4 and 5 are skipped for a given pair if step 1 fails or if the 
#   ConsumerResults file is missing.
#
#   Steps 6 and 7 run once per wagon rather than once per config. Step 7 reads the
#   raw AO2Ds instead of the consumer output, so it runs even in post-process-only
#   mode; it is skipped (not failed) when the AO2Ds/ folder is absent or empty.
#
# Usage:
#   ./run_all_wagons.sh [OPTIONS] [REGISTRY] [CONSUMER_CONFIGS_DIR]
#
# Options:
#   -h, --help                 Show this help message.
#   -p, --post-process-only    Skip Step 1 (consumer) and run only the ROOT
#                              macros (Steps 2, 3, 4, 5) on existing output.
#   -s, --skip-sig-extract     Skip Step 3 (signal extraction) during execution.
#   -f, --skip-forensics       Skip Step 7 (AO2D bit forensics). Useful when the
#                              AO2Ds have already been checked or deleted.
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
AUX_PERCONFIG_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliaryPerConfigPlots.cxx"
AUXILIARY_PLOTS_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliarySummaryPlots.cxx"
FORENSICS_MACRO="${REPO_DIR}/RingPol_RAW_LocalHelpers/zvtxBitForensics.cxx"

# Absolute path to the consumer launcher script.
FRAMEWORK_DIR="${REPO_DIR}/RingPol_RAW_LocalHelpers/DerivedDataHY"
CONSUMER_SCRIPT="${FRAMEWORK_DIR}/runDerivedDataConsumer_HY.sh"
DEFAULT_REGISTRY="${FRAMEWORK_DIR}/train_registry.conf"
# Directory containing dpl-config-DerivedConsumer-*.json files.
# Can be overridden by passing a path as $2.
DEFAULT_CONFIGS_DIR="/home/users/cicerodm/RingPol/consumer_configs"

# Directory containing the MC reference ConsumerResults files for Auxiliary Plots
# Leave empty ("") if you do not want to overlay MC.
MC_REF_DIR="/home/users/cicerodm/RingPol/LHC25h3c/ITSandTPC_min3ITS/results_consumer"

# Same thing, but for the pp baseline dataset:
# PP_REF_DIR="/home/users/cicerodm/RingPol/LHC23_pass4_Thin_small/ITSandTPC_min3ITS/results_consumer"
PP_REF_DIR="/home/users/cicerodm/RingPol/LHC23_pass4_Thin_medium/ITSandTPC-pp/results_consumer"


# Absolute path to the Toy Model ROOT file for Auxiliary Plots (chose the Toy Model representative whose configurations are 
# closest to what an actual data selection would do in data processing, as an attempt to keep everything consistent)
# Leave empty ("") if you do not want to overlay the Toy Model.
TOY_MODEL_PATH="/home/users/cicerodm/RingPol/HelicityToyModel/9_RealisticAlice/alice_an_std/helicity_alice_an_std.root"
# Testing with other model configurations:
# TOY_MODEL_PATH="/home/users/cicerodm/RingPol/HelicityToyModel/9_RealisticAlice/alice_std/helicity_alice_std.root"

# Set executable paths
EXTRACT_DELTA_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/extractDeltaErrors.exe"
SIGNAL_EXTRACT_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/signalExtractionRing.exe"
CUMUL_DCA_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/makeCumulativeDCAdauProfile.exe"
AUX_PERCONFIG_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliaryPerConfigPlots.exe"
AUXILIARY_PLOTS_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/auxiliarySummaryPlots.exe"
FORENSICS_EXE="${REPO_DIR}/RingPol_RAW_LocalHelpers/zvtxBitForensics.exe"

# Parallel workers for Step 7. Only a handful of branches are read per file
# (SetBranchStatus("*", 0)), so the per-file cost is basket decompression rather
# than streaming whole AO2Ds. 24 is about enough for a SATA SSD and leaves ample headroom
# on NVMe; raising it further only helps if the step proves CPU-bound.
FORENSICS_JOBS=24

# Set by Step 7 so the EXIT trap can remove the batch folder after an interrupt.
FORENSICS_TMP_DIR=""

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================
POST_PROCESS_ONLY=0
SKIP_SIG_EXTRACT=0
SKIP_FORENSICS=0
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
        --skip-forensics|-f)
            SKIP_FORENSICS=1
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
  rm -f "$EXTRACT_DELTA_EXE" "$SIGNAL_EXTRACT_EXE" "$CUMUL_DCA_EXE" "$AUX_PERCONFIG_EXE" "$AUXILIARY_PLOTS_EXE" "$FORENSICS_EXE"
  # Step 7's batch folder is temporary by design; remove it even on interrupt.
  if [ -n "$FORENSICS_TMP_DIR" ] && [ -d "$FORENSICS_TMP_DIR" ]; then
    rm -rf "$FORENSICS_TMP_DIR"
  fi
}

trap cleanup EXIT
trap handle_interrupt INT TERM

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================
for REQUIRED in "$REGISTRY" "$CONSUMER_SCRIPT" \
                "$EXTRACT_DELTA_MACRO" "$SIGNAL_EXTRACT_MACRO" "$CUMUL_DCA_MACRO" "$AUX_PERCONFIG_MACRO" "$AUXILIARY_PLOTS_MACRO" "$FORENSICS_MACRO"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: required file not found: ${REQUIRED}"
    exit 1
  fi
done

if [ ! -d "$CONSUMER_CONFIGS_DIR" ]; then
  echo "Error: CONSUMER_CONFIGS_DIR not found: ${CONSUMER_CONFIGS_DIR}"
  exit 1
fi

# GNU parallel is preferred for Step 7: it runs in the foreground, so Ctrl+C
# propagates to the workers normally, and it provides an ETA.
USE_PARALLEL=0
if command -v parallel >/dev/null 2>&1; then
  USE_PARALLEL=1
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
OPTIMIZATION_FLAGS="-O3 -march=native -fno-math-errno -flto=auto -fno-trapping-math" # Using the fno-math-errno to disable the errno verification of passing negative numbers to sqrts and the such.
                                                                                     # My code should be safe, so this will never be used!
                                                                                     # Also introduced -fno-trapping-math as there is no explicit usage of NaNs
                                                                                     # nor any expected unsafe divisions by zero appearing anywhere in the code!
ROOT_LIBS="$(root-config --glibs)"

# Compile Step (Stops script if compilation fails)
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$EXTRACT_DELTA_EXE" "$EXTRACT_DELTA_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$SIGNAL_EXTRACT_EXE" "$SIGNAL_EXTRACT_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$CUMUL_DCA_EXE" "$CUMUL_DCA_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$AUX_PERCONFIG_EXE" "$AUX_PERCONFIG_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$AUXILIARY_PLOTS_EXE" "$AUXILIARY_PLOTS_MACRO" $ROOT_LIBS || exit 1
g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$FORENSICS_EXE" "$FORENSICS_MACRO" $ROOT_LIBS || exit 1

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
echo "  Mode         : $( [ $POST_PROCESS_ONLY -eq 1 ] && echo 'POST-PROCESS ONLY' || echo 'FULL CHAIN' )$( [ $SKIP_SIG_EXTRACT -eq 1 ] && echo ' [SKIP SIG EXTRACT]' )$( [ $SKIP_FORENSICS -eq 1 ] && echo ' [SKIP FORENSICS]' )"
echo "  Forensics    : $( [ $SKIP_FORENSICS -eq 1 ] && echo 'skipped' || echo "${FORENSICS_JOBS} workers$( [ $USE_PARALLEL -eq 1 ] && echo ' (GNU parallel, with ETA)' || echo ' (bash fallback)' )" )"
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

    # Defining folders for each post-processing step:
    SIGNAL_EXTRACT_DIR="${WORK_DIR}/results_SigExtract"
    DELTA_ERR_DIR="${WORK_DIR}/results_DeltaErr"
    CUMUL_DIR="${WORK_DIR}/results_CumulativePlots"
    AUX_PERCONFIG_DIR="${WORK_DIR}/results_AuxPerConfig"

    # Logging Setup
    # Each step's log lives under its own "results_*/logs/" folder
    # Only the consumer's own wrapper/batch logs (and the wagon-level auxiliarySummaryPlots
    # log, see Step 6 below) remain under results_consumer/logs/.
    LOG_DIR="${WORK_DIR}/results_consumer/logs"
    DELTA_LOG_DIR="${DELTA_ERR_DIR}/logs"
    SIG_LOG_DIR="${SIGNAL_EXTRACT_DIR}/logs"
    CUMUL_LOG_DIR="${CUMUL_DIR}/logs"
    AUX_PERCONFIG_LOG_DIR="${AUX_PERCONFIG_DIR}/logs"

    mkdir -p "${LOG_DIR}"
      # Creating smaller log folders to keep everything tidy (it was getting really messy!)
    mkdir -p "${LOG_DIR}/wrappers/"
    mkdir -p "${LOG_DIR}/batches/"
    mkdir -p "${DELTA_LOG_DIR}"
    mkdir -p "${SIG_LOG_DIR}"
    mkdir -p "${CUMUL_LOG_DIR}"
    mkdir -p "${AUX_PERCONFIG_LOG_DIR}"
    WRAPPER_LOG="${LOG_DIR}/wrappers/wrapper_${CONS_SUFFIX}.log"
    DELTA_LOG="${DELTA_LOG_DIR}/extractDeltaErr_${CONS_SUFFIX}.log"
    SIG_LOG="${SIG_LOG_DIR}/sigExtract_${CONS_SUFFIX}.log"
    CUMUL_LOG="${CUMUL_LOG_DIR}/cumulDCA_${CONS_SUFFIX}.log"
    AUX_PERCONFIG_LOG="${AUX_PERCONFIG_LOG_DIR}/auxPerConfig_${CONS_SUFFIX}.log"

    # ------------------------------------------------------------------
    # Step 1: consumer
    # ------------------------------------------------------------------
    if [ $POST_PROCESS_ONLY -eq 0 ]; then
      echo -n "  [1/7] consumer        : ${CONS_SUFFIX}"
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
      echo -n "  [1/7] consumer        : ${CONS_SUFFIX} (SKIPPED)"
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
    echo -n "  [2/7] extractDeltaErr : ${CONS_SUFFIX}"
    "$EXTRACT_DELTA_EXE" "${CONSUMER_RESULT}" "${DELTA_ERR_DIR}/" > "$DELTA_LOG" 2>&1
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
      echo -n "  [3/7] sigExtract      : ${CONS_SUFFIX}"
      "$SIGNAL_EXTRACT_EXE" "${CONSUMER_RESULT}" "${SIGNAL_EXTRACT_DIR}/" > "$SIG_LOG" 2>&1
      SIG_EXIT=$?

      if [ $SIG_EXIT -ne 0 ]; then
        echo "  -> FAILED  (log: ${SIG_LOG})"
        FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | signalExtractionRing")
      else
        echo "  -> OK"
      fi
    else
      echo "  [3/7] sigExtract      : ${CONS_SUFFIX} (SKIPPED)"
    fi

    # ------------------------------------------------------------------
    # Step 4: makeCumulativeDCAdauProfile
    # ------------------------------------------------------------------
    echo -n "  [4/7] cumulDCA        : ${CONS_SUFFIX}"
    "$CUMUL_DCA_EXE" "${CONSUMER_RESULT}" "${CUMUL_DIR}/" > "$CUMUL_LOG" 2>&1
    CUMUL_EXIT=$?

    if [ $CUMUL_EXIT -ne 0 ]; then
      echo "  -> FAILED  (log: ${CUMUL_LOG})"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | makeCumulDCA")
    else
      echo "  -> OK"
    fi

    # ------------------------------------------------------------------
    # Step 5: auxiliaryPerConfigPlots
    # ------------------------------------------------------------------
    echo -n "  [5/7] auxPerConfig    : ${CONS_SUFFIX}"
    "$AUX_PERCONFIG_EXE" "${CONSUMER_RESULT}" "${AUX_PERCONFIG_DIR}/" > "$AUX_PERCONFIG_LOG" 2>&1
    AUX_PERCONFIG_EXIT=$?

    if [ $AUX_PERCONFIG_EXIT -ne 0 ]; then
      echo "  -> FAILED  (log: ${AUX_PERCONFIG_LOG})"
      FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ${CONS_SUFFIX} | auxPerConfigPlots")
    else
      echo "  -> OK"
    fi

    echo ""

  done  # configs

  # ------------------------------------------------------------------
  # Step 6: auxiliarySummaryPlots (Cross-configuration aggregation)
  # ------------------------------------------------------------------
  # We run this ONCE per wagon, passing the wagon's base working dir
  AUX_LOG="${WORK_DIR}/results_consumer/logs/auxSummaryPlots.log" # Saves a single log, to the root of the logs folder, because there will be only one single log for a given wagon
  echo -n "  [6/7] auxiliarySummaryPlots  : (Cross-config summary)"
  
  # Forwarding the consumer results directory, the MC reference, and now the Toy Model path:
  "$AUXILIARY_PLOTS_EXE" "${WORK_DIR}/results_consumer" "${MC_REF_DIR}" "${PP_REF_DIR}" "${TOY_MODEL_PATH}" > "$AUX_LOG" 2>&1
  AUX_EXIT=$?

  if [ $AUX_EXIT -ne 0 ]; then
    echo "  -> FAILED  (log: ${AUX_LOG})"
    FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ALL_CONFIGS | auxiliarySummaryPlots")
  else
    echo "  -> OK"
  fi

  # ------------------------------------------------------------------
  # Step 7: zvtxBitForensics (AO2D bit-level and integrity QA)
  # ------------------------------------------------------------------
  # Runs once per wagon on the raw AO2Ds, independently of the consumer output.
  # Files are distributed round-robin across NBATCHES manifests, processed by
  # FORENSICS_JOBS concurrent workers, each writing its own partial .root; hadd
  # merges them and a --finalize pass turns the additive numerator/denominator
  # counters into ratio histograms. Ratios cannot be produced before the merge
  # because hadd can only sum.
  FORENSICS_DIR="${WORK_DIR}/AO2Ds"
  FORENSICS_OUT="${WORK_DIR}/results_consumer/zvtxBitForensics.root"
  FORENSICS_LOG="${WORK_DIR}/results_consumer/logs/zvtxBitForensics.log"
  FORENSICS_LOG_DIR="${WORK_DIR}/results_consumer/logs/forensics_logs"
  mkdir -p "${WORK_DIR}/results_consumer/logs" "$FORENSICS_LOG_DIR"

  echo -n "  [7/7] zvtxBitForensics       : (AO2D integrity QA)"

  if [ $SKIP_FORENSICS -eq 1 ]; then
    echo "  -> SKIPPED (--skip-forensics)"
  else
    shopt -s nullglob
    AOD_FILES=("${FORENSICS_DIR}"/AO2D_*.root)
    shopt -u nullglob

    # EDGE CASE 1: 0 Files
    if [ ${#AOD_FILES[@]} -eq 0 ]; then
      echo "  -> SKIPPED (no AO2Ds in ${FORENSICS_DIR})"
    else
      # More batches than workers: parallel then has something to schedule,
      # which gives both load balancing and a usable ETA. Capped so the number
      # of partial .root files stays modest for hadd.
      NBATCHES=$(( FORENSICS_JOBS * 4 ))

      # EDGE CASE 2: Fewer files than max batches (e.g., 1 file -> 1 batch)
      if [ ${#AOD_FILES[@]} -lt $NBATCHES ]; then
        NBATCHES=${#AOD_FILES[@]}
      fi

      FORENSICS_TMP_DIR="${WORK_DIR}/results_consumer/.forensics_batches_$$"
      mkdir -p "$FORENSICS_TMP_DIR"

      # The batch count varies between runs, so logs from a previous, larger run
      # would otherwise linger here and be misread as belonging to this one.
      rm -f "${FORENSICS_LOG_DIR}"/batch_*.log

      IDX=0
      for AOD_FILE in "${AOD_FILES[@]}"; do
        echo "$AOD_FILE" >> "${FORENSICS_TMP_DIR}/batch_$((IDX % NBATCHES)).txt"
        IDX=$((IDX + 1))
      done

      : > "$FORENSICS_LOG"
      FORENSICS_EXIT=0

      if [ $USE_PARALLEL -eq 1 ]; then
        # parallel runs in the foreground, so Ctrl+C reaches it and its children
        # through the normal foreground process group. Worker stdout goes to
        # per-batch logs; --eta writes progress to the terminal.
        #
        # The command string is single-quoted and expanded by a FRESH shell, so
        # every variable it names must be exported -- an unexported one silently
        # becomes the empty string.
        #
        # {/.} is the manifest basename without extension (batch_0.txt -> batch_0),
        # keeping manifest, partial .root and log names in agreement. {#} would not:
        # it is the 1-based job sequence number, so batch_0.txt would produce
        # batch_1.log.
        export FORENSICS_EXE FORENSICS_TMP_DIR FORENSICS_LOG_DIR
        parallel --will-cite --eta -j "$FORENSICS_JOBS" \
          '"$FORENSICS_EXE" {} "$FORENSICS_TMP_DIR/zvtxForensics_{/.}.root" > "$FORENSICS_LOG_DIR/{/.}.log" 2>&1' \
          ::: "${FORENSICS_TMP_DIR}"/batch_*.txt
        # Captured immediately: parallel exits with the number of failed jobs.
        PARALLEL_EXIT=$?
        [ $PARALLEL_EXIT -eq 0 ] || FORENSICS_EXIT=1
      else
        # Fallback: background jobs have SIGINT set to SIG_IGN by bash, so the
        # PIDs are recorded and the trap SIGTERMs them explicitly.
        #
        # The tag comes from the manifest name, not a counter: the glob expands
        # lexically (batch_0, batch_1, batch_10, batch_2, ...), so a counter would
        # pair each log with the wrong manifest.
        for BATCH_MANIFEST in "${FORENSICS_TMP_DIR}"/batch_*.txt; do
          BATCH_TAG=$(basename "$BATCH_MANIFEST" .txt)
          "$FORENSICS_EXE" "$BATCH_MANIFEST" "${FORENSICS_TMP_DIR}/zvtxForensics_${BATCH_TAG}.root" \
            > "${FORENSICS_LOG_DIR}/${BATCH_TAG}.log" 2>&1 &
          FORENSICS_PIDS+=($!)

          # Throttle to FORENSICS_JOBS concurrent workers.
          while [ "$(jobs -rp | wc -l)" -ge "$FORENSICS_JOBS" ]; do
            sleep 0.2
          done
        done
        for PID in "${FORENSICS_PIDS[@]}"; do
          wait "$PID" || FORENSICS_EXIT=1
        done
        FORENSICS_PIDS=()
      fi

      # Merge, then finalize. Both must succeed for the step to count as OK.
      if [ $FORENSICS_EXIT -eq 0 ]; then
        shopt -s nullglob
        BATCH_ROOTS=("${FORENSICS_TMP_DIR}"/zvtxForensics_batch_*.root)
        shopt -u nullglob

        if [ ${#BATCH_ROOTS[@]} -eq 0 ]; then
          FORENSICS_EXIT=1
        # EDGE CASE 3: Exactly 1 output file. Skip hadd and just move it.
        elif [ ${#BATCH_ROOTS[@]} -eq 1 ]; then
          mv "${BATCH_ROOTS[0]}" "$FORENSICS_OUT" || FORENSICS_EXIT=1
        else
          hadd -f "$FORENSICS_OUT" "${BATCH_ROOTS[@]}" >> "$FORENSICS_LOG" 2>&1 || FORENSICS_EXIT=1
        fi
      fi

      if [ $FORENSICS_EXIT -eq 0 ]; then
        "$FORENSICS_EXE" --finalize "$FORENSICS_OUT" >> "$FORENSICS_LOG" 2>&1 || FORENSICS_EXIT=1
      fi

      rm -rf "$FORENSICS_TMP_DIR"
      FORENSICS_TMP_DIR=""

      if [ $FORENSICS_EXIT -ne 0 ]; then
        echo "  -> FAILED  (log: ${FORENSICS_LOG})"
        FAILURES+=("${DATASET_NAME}/${WAGON_SHORTNAME} | ALL_CONFIGS | zvtxBitForensics")
      else
        echo "  -> OK"
      fi
    fi
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
  echo "  (Each step's log lives under its own results_*/logs/ folder inside the wagon"
  echo "   directory -- e.g. results_DeltaErr/logs/, results_SigExtract/logs/, etc. -- except"
  echo "   the consumer and auxiliarySummaryPlots logs, which stay under results_consumer/logs/)"
fi
echo "========================================================"