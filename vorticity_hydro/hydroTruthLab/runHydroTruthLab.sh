#!/usr/bin/env bash

###############################################################################
# Script name:
#   runHydroTruthLab.sh
#
# Purpose:
#   Runs the hydroTruthLab chain for every run of the runs registry, each on one dataset of the datasets registry (Vitor's hydro events):
#     1. hydroTruthLab.cxx (one process per sample)  -> <run>/hydroTruthLab_<run>_bullet.root, hydroTruthLab_<run>_noBullet.root
#     2. hydroTruthLabPerConfig.cxx (both samples)   -> <run>/hydroTruthLabPerConfig_<run>.root
#     3. hydroTruthLabSummaryPlots.cxx (every run)   -> hydroTruthLabSummaryPlots.root, the compare lines of the runs registry
#
#   with <run> = <dataset>_<variation>. The defaults of every configurable live in the CONFIGURATION block below, and a run of the registry
#   lists only what it changes. Runs go in parallel, MAX_PARALLEL_RUNS at a time, each with its two producers at once, sharing one seed.
#   Each sample runs to its last event with hydro output (and to the last entry of the jet file for the jet trigger), except with mixing,
#   where both samples take the range of the bullet sample. The producers skip the events of the range without hydro output, with a warning.
#   Design, definitions and the output reference: README.md in this folder.
#
# Usage:
#   ./runHydroTruthLab.sh [OPTIONS]
#
# Options:
#   -h, --help                 Show this help message.
#   -p, --post-process-only    Skip the producers and rerun PerConfig and SummaryPlots on the existing producer files.
#   -s, --summary-plots-only   Rerun SummaryPlots only, on the existing PerConfig files.
#   -t, --tag NAME             Name of the campaign folder under OUT_BASE (default: "default").
#   -o, --only RUN             Run only this run (<dataset>_<variation>). Repeatable.
#   -j, --jobs N               Runs at a time (default: MAX_PARALLEL_RUNS). 1 runs them one after the other.
#       --runs FILE            Runs registry (default: hydroTruthLab_runs.conf next to this script).
#       --datasets FILE        Datasets registry (default: hydroTruthLab_datasets.conf next to this script).
#       --cleanup              Remove the compiled executables and exit. g++ rebuilds them on every run anyway,
#                              so a stale binary can never be picked up.
#
# Output (inside ${OUT_BASE}/${RUN_TAG}/):
#   logHydroTruthLab.log                  this script's log (appended, one block per invocation)
#   backupOf_runHydroTruthLab.sh          copies of this script and of both registries as they were last run: the complete configuration.
#   backupOf_hydroTruthLab_runs.conf      Records only, not meant to be run from there. A -p run overwrites them, but the producer files
#   backupOf_hydroTruthLab_datasets.conf  keep their own options in Config/ (and PerConfig copies them into its own Config/)
#   <run>/                                one flat folder per run, every file named after it:
#                                           hydroTruthLab_<run>_bullet.root (+ .log), hydroTruthLab_<run>_noBullet.root (+ .log),
#                                           hydroTruthLabPerConfig_<run>.root (+ .log)
#   hydroTruthLabSummaryPlots.root (+ .log)  one folder per compare line of the runs registry
#   *.root.part                           never left behind: outputs are written under this name and renamed only on success
#
# Examples:
#   ./runHydroTruthLab.sh
#   ./runHydroTruthLab.sh --only lumpy4050_default --only lumpy6070_default
#   ./runHydroTruthLab.sh -p --tag smudgeStudy --jobs 4
#   ./runHydroTruthLab.sh -s
#
###############################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# The macros and the registries sit next to this script, so the repository path is never hardcoded
LAB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRODUCER_MACRO="${LAB_DIR}/hydroTruthLab.cxx"
PERCONFIG_MACRO="${LAB_DIR}/hydroTruthLabPerConfig.cxx"
SUMMARY_MACRO="${LAB_DIR}/hydroTruthLabSummaryPlots.cxx"
COMMON_HEADER="${LAB_DIR}/hydroTruthLabCommon.h" # included by PerConfig and SummaryPlots from this folder: no flag needed
PRODUCER_EXE="${LAB_DIR}/hydroTruthLab.exe"
PERCONFIG_EXE="${LAB_DIR}/hydroTruthLabPerConfig.exe"
SUMMARY_EXE="${LAB_DIR}/hydroTruthLabSummaryPlots.exe"
DATASETS_FILE="${LAB_DIR}/hydroTruthLab_datasets.conf"
RUNS_FILE="${LAB_DIR}/hydroTruthLab_runs.conf"

# Output base folder: the only place where it is set
OUT_BASE="${HOME}/RingPol/hydroTruthLab"

# Hydro inputs (Vitor's events in jarvis4). Sample folders and jet files of the datasets registry are relative to this folder
HYDRO_BASE="/storage1/vribeiro"
INPUT_NAME="Smu_dpTdphidy_Thermal_rapidity_3122.dat"

# Event range of a sample: from FIRST_EVENT to its last event that holds INPUT_NAME, capped by the last entry of the jet file when the jet
# trigger or mixing needs one per event. With mixing, both samples take the range of the bullet sample: the partners of an event are chosen
# among the events of the range, so the two samples must share it to stay paired
FIRST_EVENT=0
MAX_EVENTS=0 # 0 = no cap. A positive value caps the range (quick tests)

# Runs at a time: half the cores, since each run is two single-threaded producers at once (then PerConfig), about 30 s for 250 events.
# jarvis4 is shared, so lower it (or pass --jobs) when others need the machine
MAX_PARALLEL_RUNS=$(( $(nproc) / 2 > 0 ? $(nproc) / 2 : 1 ))

# Constants copied from the consumer (lambdaJetPolarizationIonsDerived.cxx): LambdaMass = o2::constants::physics::MassLambda,
# LambdaWeakDecayConstant (DPG 2025 update). Keep them in sync with it by hand: this chain cannot read the O2 headers.
LAMBDA_MASS=1.115683
ALPHA_LAMBDA=0.749

# Producer defaults (hydroTruthLab.exe --help lists them all). A run changes any of them by name. --jetFile, --inputName, --firstEvent,
# --nEvents and --seed are set per run, from the registries
PRODUCER_OPTS=(
  --trigger=jet # jet: phi_t = phi_J. rp: phi_t = 0, along the impact parameter (README, "Trigger modes"), which turns mixing off
  # Cuts of Vitor's PRC paper (PhysRevC.109.014905), repeated as a closure test of our code
  --ptMin=0.5
  --ptMax=1.5
  --absYMax=0.5
  --lambdaMass="${LAMBDA_MASS}"
  --deltaScanSteps=100 # one step per phi cell (3.6 deg)
  --mixDepth=10 # mixing partners per event. With the rp trigger every mixed trigger is phi_t = 0, so mixing is turned off there
  --smudgeSigma=0.1 # Higher than the consumer's forceJetDirectionSmudge width, which would be 5% * jetR ~ 0.02, with jetR = 0.4
  --smudgeMode=cone
  --smudgeDraws=10
  --mapSubdivisions=3 # PolarizationMaps/: each cell spread over 3x3x3 points of its momentum volume (1 = nodes only, with empty bins at high p)
)

# PerConfig defaults (hydroTruthLabPerConfig.exe --help lists them all). --centralityLabel, --initialConditions and --runLabel come from the
# registries
PERCONFIG_OPTS=(
  --paired=1 # key hypothesis: paired initial conditions (README)
  --alphaLambda="${ALPHA_LAMBDA}"
  # Lambdas left after every cut in the O-O analysis (topological, kinematic, PID, Armenteros): the significance of a signal as pure as hydro
  # allows. A ceiling, not a prediction: Pb-Pb hydro, ideal detector, one measurement (README, "Expected precision")
  --nLambdaPrecision=288574517
)

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

# Numbers are written and parsed in the C locale whatever the user's: a decimal comma would break awk and the macros' options
export LC_ALL=C
shopt -s nullglob # a glob that matches nothing expands to nothing, not to itself

POST_PROCESS_ONLY=0
SUMMARY_ONLY=0
CLEANUP_ONLY=0
RUN_TAG="default"
JOBS=$MAX_PARALLEL_RUNS
ONLY_RUNS=()

# Prints the header block from "# Usage:" down to the closing line of hashes.
# The end marker is the full hash line, not "^#+$": a bare "#" separator line also matches that pattern and would cut the help short.
print_help() {
  awk '/^# Usage:/{show=1} show && /^##########/{exit} show' "$0" | sed 's/^# \{0,1\}//'
  exit 0
}

# Value of an option that takes one: stops when it is missing
option_argument() {
  if [ -z "$2" ]; then
    echo "Error: $1 needs a value. Use --help for usage." >&2
    exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      print_help
      ;;
    --post-process-only|-p)
      POST_PROCESS_ONLY=1
      shift
      ;;
    --summary-plots-only|-s)
      SUMMARY_ONLY=1
      shift
      ;;
    --tag|-t)
      option_argument "$1" "$2"
      RUN_TAG="$2"
      shift 2
      ;;
    --only|-o)
      option_argument "$1" "$2"
      ONLY_RUNS+=("$2")
      shift 2
      ;;
    --jobs|-j)
      option_argument "$1" "$2"
      if ! [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
        echo "Error: $1 needs a positive integer, got '$2'" >&2
        exit 1
      fi
      JOBS="$2"
      shift 2
      ;;
    --runs)
      option_argument "$1" "$2"
      RUNS_FILE="$2"
      shift 2
      ;;
    --datasets)
      option_argument "$1" "$2"
      DATASETS_FILE="$2"
      shift 2
      ;;
    --cleanup)
      CLEANUP_ONLY=1
      shift
      ;;
    *)
      echo "Unknown option: $1. Use --help for usage." >&2
      exit 1
      ;;
  esac
done

if [ $POST_PROCESS_ONLY -eq 1 ] && [ $SUMMARY_ONLY -eq 1 ]; then
  echo "Error: -p and -s exclude each other (-p already reruns SummaryPlots). Use --help for usage." >&2
  exit 1
fi

if [ $CLEANUP_ONLY -eq 1 ]; then
  echo "Removing compiled macro executables..."
  rm -f "$PRODUCER_EXE" "$PERCONFIG_EXE" "$SUMMARY_EXE"
  echo "Done."
  exit 0
fi

# ==============================================================================
# CAMPAIGN FOLDER AND LOGGING
# ==============================================================================
CAMPAIGN_DIR="${OUT_BASE}/${RUN_TAG}"
# SummaryPlots gets its files as comma-separated lists
if [[ "$CAMPAIGN_DIR" == *,* ]]; then
  echo "Error: the campaign folder ${CAMPAIGN_DIR} contains a comma" >&2
  exit 1
fi
mkdir -p "$CAMPAIGN_DIR" || exit 1

# # Appended rather than overwritten: a post-process-only rerun keeps the log of the run that made the producer files
# exec > >(tee -a -i "${CAMPAIGN_DIR}/logHydroTruthLab.log") 2>&1
# Overwritten rather than appended -- way cleaner!
exec > >(tee -i "${CAMPAIGN_DIR}/logHydroTruthLab.log") 2>&1
echo ""
echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"

# ==============================================================================
# SHARED RUN STATE
# ==============================================================================
# Failure and skip logs live in files, not bash arrays: the runs go in background subshells, which cannot append to the parent's arrays.
# Each entry is one short line written with ">>", which the kernel appends atomically. Format: "RUN | STEP".
RUN_TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/runHydroTruthLab.XXXXXX")
FAILURES_FILE="${RUN_TMP_DIR}/failures.txt"
SKIPPED_FILE="${RUN_TMP_DIR}/skipped.txt"
: > "$FAILURES_FILE"
: > "$SKIPPED_FILE"

LANE_PIDS=()

record_failure() { echo "$1" >> "$FAILURES_FILE"; }
record_skip()    { echo "$1" >> "$SKIPPED_FILE"; }

# Prints one "RUN | STEP" table.
#   $1        : heading for the second column ("FAILED STEP", "SKIPPED STEP")
#   $2 .. $n  : the entries
print_report_table() {
  local STEP_HEADER="$1"
  shift
  printf "  %-32s  %-60s\n" "RUN" "$STEP_HEADER"
  printf "  %-32s  %-60s\n" "--------------------------------" "------------------------------------------------------------"
  local ENTRY RUN_COL STEP_COL
  for ENTRY in "$@"; do
    RUN_COL=$( echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
    STEP_COL=$(echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
    printf "  %-32s  %-60s\n" "$RUN_COL" "$STEP_COL"
  done
}

# ==============================================================================
# SIGNAL HANDLING
# ==============================================================================
# Background jobs of a non-interactive script ignore SIGINT, so Ctrl+C does not reach the runs by themselves. Each run is a subshell with its
# producers in subshells of their own, so the whole tree below each run is terminated. Every process is stopped before its children are
# listed: a run whose producer had just been killed would otherwise move on and start PerConfig, out of the list. A stopped process keeps the
# SIGTERM pending until SIGCONT
kill_tree() {
  local CHILD
  kill -STOP "$1" 2>/dev/null
  for CHILD in $(pgrep -P "$1"); do
    kill_tree "$CHILD"
  done
  kill -TERM "$1" 2>/dev/null
  kill -CONT "$1" 2>/dev/null
}

handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  if [ ${#LANE_PIDS[@]} -gt 0 ]; then
    echo "    Stopping the runs..."
    local LANE_PID
    for LANE_PID in "${LANE_PIDS[@]}"; do
      kill_tree "$LANE_PID"
    done
    wait "${LANE_PIDS[@]}" 2>/dev/null
  fi
  # Whatever was running when the signal came left only "<file>.part" behind: never a final file
  find "$CAMPAIGN_DIR" -name '*.root.part' -delete 2>/dev/null
  echo "    Stopping. Partial failure report:"
  echo ""
  mapfile -t FAILURES < "$FAILURES_FILE"
  if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "  No failures recorded before interrupt."
  else
    print_report_table "FAILED STEP" "${FAILURES[@]}"
  fi
  echo "========================================================"
  exit 130
}

cleanup() {
  if [ -n "$RUN_TMP_DIR" ] && [ -d "$RUN_TMP_DIR" ]; then
    rm -rf "$RUN_TMP_DIR"
  fi
}

trap cleanup EXIT
trap handle_interrupt INT TERM

# ==============================================================================
# HELPERS
# ==============================================================================

# Path of the registries: relative to HYDRO_BASE unless absolute
hydro_path() {
  if [[ "$1" = /* ]]; then echo "$1"; else echo "${HYDRO_BASE}/$1"; fi
}

# Reads integers, one per line, and prints them as ranges ("28-29, 33-34, 49-199"), as the producer does in its log
compress_ids() {
  sort -n | awk '
    NR == 1 { first = $1; last = $1; next }
    $1 == last + 1 { last = $1; next }
    { text = text (text == "" ? "" : ", ") (first == last ? first : first "-" last); first = $1; last = $1 }
    END { if (NR > 0) print text (text == "" ? "" : ", ") (first == last ? first : first "-" last) }'
}

# True when the named array of --key=value words has option KEY
has_option() {
  local -n OPTION_WORDS="$1"
  local WORD
  for WORD in "${OPTION_WORDS[@]}"; do
    [[ "$WORD" == "--$2="* ]] && return 0
  done
  return 1
}

# Value of option KEY in the named array of --key=value words (empty when absent)
option_value() {
  local -n OPTION_WORDS="$1"
  local WORD
  for WORD in "${OPTION_WORDS[@]}"; do
    if [[ "$WORD" == "--$2="* ]]; then
      echo "${WORD#*=}"
      return
    fi
  done
}

# Replaces the option of the same name in the named array by KEY=VALUE ($2); nothing when the array has no such option
replace_option() {
  local -n OPTION_WORDS="$1"
  local INDEX
  for INDEX in "${!OPTION_WORDS[@]}"; do
    [[ "${OPTION_WORDS[$INDEX]}" == "--${2%%=*}="* ]] && OPTION_WORDS[INDEX]="--$2"
  done
}

# Effective options of a run, in RUN_PRODUCER_OPTS and RUN_PERCONFIG_OPTS: the defaults, changed by the overrides of the run (each applies to
# the macro that has an option of that name), and mixing off for the rp trigger. The seed is set per run elsewhere
resolve_run_options() {
  local OVERRIDES=() OVERRIDE
  RUN_PRODUCER_OPTS=("${PRODUCER_OPTS[@]}")
  RUN_PERCONFIG_OPTS=("${PERCONFIG_OPTS[@]}")
  read -r -a OVERRIDES <<< "${RUN_OVERRIDES[$1]}"
  for OVERRIDE in "${OVERRIDES[@]}"; do
    [[ "${OVERRIDE%%=*}" == "seed" ]] && continue
    replace_option RUN_PRODUCER_OPTS "$OVERRIDE"
    replace_option RUN_PERCONFIG_OPTS "$OVERRIDE"
  done
  if [ "$(option_value RUN_PRODUCER_OPTS trigger)" = "rp" ]; then
    replace_option RUN_PRODUCER_OPTS "mixDepth=0"
  fi
}

# Last event of a sample folder, from FIRST_EVENT on, that holds INPUT_NAME (empty when none)
last_event_with_input() {
  local EVENT_DIR ID LAST=""
  for EVENT_DIR in "$1"/*/; do
    ID="${EVENT_DIR%/}"
    ID="${ID##*/}" # the folder name, without a basename process per event
    [[ "$ID" =~ ^[0-9]+$ ]] || continue
    ID=$((10#$ID))
    (( ID >= FIRST_EVENT )) || continue
    [ -f "${EVENT_DIR}results/${INPUT_NAME}" ] || continue
    if [ -z "$LAST" ] || (( ID > LAST )); then
      LAST=$ID
    fi
  done
  echo "$LAST"
}

# Events of FIRST_EVENT .. $2 without INPUT_NAME in sample folder $1, as ranges (empty when none)
missing_events() {
  local ID
  for (( ID = FIRST_EVENT; ID <= $2; ++ID )); do
    [ -f "$1/${ID}/results/${INPUT_NAME}" ] || echo "$ID"
  done | compress_ids
}

# ==============================================================================
# REGISTRIES
# ==============================================================================
for REQUIRED in "$PRODUCER_MACRO" "$PERCONFIG_MACRO" "$SUMMARY_MACRO" "$COMMON_HEADER" "$DATASETS_FILE" "$RUNS_FILE"; do
  if [ ! -f "$REQUIRED" ]; then
    echo "Error: not found: ${REQUIRED}"
    exit 1
  fi
done
# Labels: no "_", which separates the dataset from the variation in the run label
LABEL_PATTERN='^[A-Za-z0-9.+-]+$'

# Datasets: label centrality initialConditions bullet noBullet jetFile
declare -A DATASET_CENTRALITY DATASET_INITIAL_CONDITIONS DATASET_BULLET DATASET_NO_BULLET DATASET_JET_FILE
LINE_NUMBER=0
while IFS= read -r LINE || [ -n "$LINE" ]; do
  LINE_NUMBER=$((LINE_NUMBER + 1))
  read -r -a FIELDS <<< "${LINE%%#*}"
  [ ${#FIELDS[@]} -eq 0 ] && continue
  if [ ${#FIELDS[@]} -ne 6 ] || ! [[ "${FIELDS[0]}" =~ $LABEL_PATTERN ]] || [ -n "${DATASET_BULLET[${FIELDS[0]}]}" ]; then
    echo "Error: ${DATASETS_FILE}, line ${LINE_NUMBER}: expected 'label centrality initialConditions bullet noBullet jetFile' with a new label of letters, digits, '.', '+' or '-'"
    exit 1
  fi
  DATASET_CENTRALITY[${FIELDS[0]}]="${FIELDS[1]}"
  DATASET_INITIAL_CONDITIONS[${FIELDS[0]}]="${FIELDS[2]}"
  DATASET_BULLET[${FIELDS[0]}]="${FIELDS[3]}"
  DATASET_NO_BULLET[${FIELDS[0]}]="${FIELDS[4]}"
  DATASET_JET_FILE[${FIELDS[0]}]="${FIELDS[5]}"
done < "$DATASETS_FILE"

# Runs: "run dataset variation [option=value ...]", and the comparisons for SummaryPlots, "compare name x run ...", only checked here
declare -A RUN_DATASET RUN_OVERRIDES
RUN_LABELS=()
COMPARE_LINES=()
LINE_NUMBER=0
while IFS= read -r LINE || [ -n "$LINE" ]; do
  LINE_NUMBER=$((LINE_NUMBER + 1))
  read -r -a FIELDS <<< "${LINE%%#*}"
  [ ${#FIELDS[@]} -eq 0 ] && continue
  if [ "${FIELDS[0]}" = "compare" ]; then
    COMPARE_LINES+=("${LINE_NUMBER}|${FIELDS[*]}")
    continue
  fi
  if [ "${FIELDS[0]}" != "run" ] || [ ${#FIELDS[@]} -lt 3 ] || ! [[ "${FIELDS[2]}" =~ $LABEL_PATTERN ]]; then
    echo "Error: ${RUNS_FILE}, line ${LINE_NUMBER}: expected 'run dataset variation [option=value ...]', variation of letters, digits, '.', '+' or '-'"
    exit 1
  fi
  if [ -z "${DATASET_BULLET[${FIELDS[1]}]}" ]; then
    echo "Error: ${RUNS_FILE}, line ${LINE_NUMBER}: dataset '${FIELDS[1]}' is not in ${DATASETS_FILE}"
    exit 1
  fi
  LABEL="${FIELDS[1]}_${FIELDS[2]}"
  if [ -n "${RUN_DATASET[$LABEL]}" ]; then
    echo "Error: ${RUNS_FILE}, line ${LINE_NUMBER}: run ${LABEL} appears twice"
    exit 1
  fi
  # Every override must name an option of a macro, or the seed: a typo must stop here, not run with the default
  for OVERRIDE in "${FIELDS[@]:3}"; do
    KEY="${OVERRIDE%%=*}"
    if [[ "$OVERRIDE" != *=?* ]] || { [ "$KEY" != "seed" ] && ! has_option PRODUCER_OPTS "$KEY" && ! has_option PERCONFIG_OPTS "$KEY"; }; then
      echo "Error: ${RUNS_FILE}, line ${LINE_NUMBER}: '${OVERRIDE}' is not option=value for an option of PRODUCER_OPTS, PERCONFIG_OPTS or the seed"
      exit 1
    fi
    # A seed of digits is a fixed seed, and the producer reads it as a 32-bit unsigned integer. Anything else names a seed group
    if [ "$KEY" = "seed" ] && [[ "${OVERRIDE#seed=}" =~ ^[0-9]+$ ]] && { [ ${#OVERRIDE} -gt 15 ] || [ "${OVERRIDE#seed=}" -gt 4294967295 ]; }; then
      echo "Error: ${RUNS_FILE}, line ${LINE_NUMBER}: '${OVERRIDE}' is above 4294967295, the largest 32-bit seed"
      exit 1
    fi
  done
  RUN_LABELS+=("$LABEL")
  RUN_DATASET[$LABEL]="${FIELDS[1]}"
  RUN_OVERRIDES[$LABEL]="${FIELDS[*]:3}"
done < "$RUNS_FILE"

# Comparisons can only name runs of the registry, along an x that every PerConfig file records in Config/: an option of either macro, the
# centrality, the initial conditions or the run itself (categorical, in the order listed)
for ENTRY in "${COMPARE_LINES[@]}"; do
  read -r -a FIELDS <<< "${ENTRY#*|}"
  if [ ${#FIELDS[@]} -lt 5 ]; then
    echo "Error: ${RUNS_FILE}, line ${ENTRY%%|*}: expected 'compare name x run run ...' with two runs at least"
    exit 1
  fi
  if ! has_option PRODUCER_OPTS "${FIELDS[2]}" && ! has_option PERCONFIG_OPTS "${FIELDS[2]}" && [[ ! "${FIELDS[2]}" =~ ^(centralityLabel|initialConditions|runLabel|seed)$ ]]; then
    echo "Error: ${RUNS_FILE}, line ${ENTRY%%|*}: x '${FIELDS[2]}' is neither an option of PRODUCER_OPTS or PERCONFIG_OPTS nor centralityLabel, initialConditions, runLabel or seed"
    exit 1
  fi
  for LABEL in "${FIELDS[@]:3}"; do
    if [ -z "${RUN_DATASET[$LABEL]}" ]; then
      echo "Error: ${RUNS_FILE}, line ${ENTRY%%|*}: comparison ${FIELDS[1]} names run ${LABEL}, which is not in the registry"
      exit 1
    fi
  done
done

# Runs of this invocation
SELECTED_RUNS=()
if [ ${#ONLY_RUNS[@]} -eq 0 ]; then
  SELECTED_RUNS=("${RUN_LABELS[@]}")
else
  for LABEL in "${ONLY_RUNS[@]}"; do
    if [ -z "${RUN_DATASET[$LABEL]}" ]; then
      echo "Error: --only ${LABEL}: no such run in ${RUNS_FILE}"
      exit 1
    fi
    SELECTED_RUNS+=("$LABEL")
  done
fi
if [ ${#SELECTED_RUNS[@]} -eq 0 ]; then
  echo "Error: no run in ${RUNS_FILE}"
  exit 1
fi

cp "${BASH_SOURCE[0]}" "${CAMPAIGN_DIR}/backupOf_runHydroTruthLab.sh"
cp "$RUNS_FILE" "${CAMPAIGN_DIR}/backupOf_hydroTruthLab_runs.conf"
cp "$DATASETS_FILE" "${CAMPAIGN_DIR}/backupOf_hydroTruthLab_datasets.conf"

# ==============================================================================
# PLAN: event range and seed of every run
# ==============================================================================
# Seeds: a run without seed= draws its own from the clock, seed=<integer> fixes it, and seed=<name> shares one clock-drawn seed among every
# run with that name in this invocation (the same smudge draws, rescaled by each width). Both samples of a run always share theirs, which
# keeps their smudged triggers paired (README, key hypothesis)
declare -A RUN_FIRST_EVENT RUN_N_EVENTS_BULLET RUN_N_EVENTS_NO_BULLET RUN_SEED SEED_OF_GROUP
draw_seed() { echo $(( $(date +%s%N) % 4294967296 )); }

# The seed= value of a run (empty when it has none)
seed_spec() {
  local OVERRIDES=() OVERRIDE
  read -r -a OVERRIDES <<< "${RUN_OVERRIDES[$1]}"
  for OVERRIDE in "${OVERRIDES[@]}"; do
    [[ "$OVERRIDE" == seed=* ]] && echo "${OVERRIDE#seed=}"
  done
}

# Last event of the range of one sample folder ($1), with the jet file of the run ($2): its last event with iSS output, capped at the last
# entry of the jet file when the jet trigger or mixing needs one per event, and at MAX_EVENTS. Prints "<last event>|<caps applied>", or an
# error on stderr and a nonzero status
range_end() {
  local SAMPLE_DIR="$1" JET_FILE="$2" LAST LAST_JET CAPS=""
  LAST=$(last_event_with_input "$SAMPLE_DIR")
  if [ -z "$LAST" ]; then
    echo "Error: no event of ${SAMPLE_DIR} from ${FIRST_EVENT} on holds ${INPUT_NAME}" >&2
    return 1
  fi
  if [ "$(option_value RUN_PRODUCER_OPTS trigger)" = "jet" ] || [ "$(option_value RUN_PRODUCER_OPTS mixDepth)" != "0" ]; then
    LAST_JET=$(awk '!/^[[:space:]]*#/ && NF {print int($1)}' "$JET_FILE" | sort -n | tail -1)
    if [ -z "$LAST_JET" ] || [ "$LAST_JET" -lt "$FIRST_EVENT" ]; then
      echo "Error: ${JET_FILE} has no event from ${FIRST_EVENT} on" >&2
      return 1
    fi
    if [ "$LAST" -gt "$LAST_JET" ]; then
      LAST=$LAST_JET
      CAPS=", capped by the jet file"
    fi
  fi
  if [ "$MAX_EVENTS" -gt 0 ] && [ $((LAST - FIRST_EVENT + 1)) -gt "$MAX_EVENTS" ]; then
    LAST=$((FIRST_EVENT + MAX_EVENTS - 1))
    CAPS="${CAPS}, capped by MAX_EVENTS"
  fi
  echo "${LAST}|${CAPS}"
}

# Checks the inputs of a run and sets the event range of each sample and the seed. Stops the script on a missing input: nothing has run yet
plan_run() {
  local LABEL="$1"
  local DATASET="${RUN_DATASET[$LABEL]}"
  local BULLET_DIR NO_BULLET_DIR="" JET_FILE SAMPLE_DIR RANGE LAST_BULLET LAST_NO_BULLET CAPS_BULLET CAPS_NO_BULLET="" SEED_SPEC MISSING
  resolve_run_options "$LABEL"
  BULLET_DIR=$(hydro_path "${DATASET_BULLET[$DATASET]}")
  [ "${DATASET_NO_BULLET[$DATASET]}" != "none" ] && NO_BULLET_DIR=$(hydro_path "${DATASET_NO_BULLET[$DATASET]}")
  JET_FILE=$(hydro_path "${DATASET_JET_FILE[$DATASET]}")
  for SAMPLE_DIR in "$BULLET_DIR" "$NO_BULLET_DIR"; do
    if [ -n "$SAMPLE_DIR" ] && [ ! -d "$SAMPLE_DIR" ]; then
      echo "Error: ${LABEL}: sample folder not found: ${SAMPLE_DIR}"
      exit 1
    fi
  done
  if [ ! -f "$JET_FILE" ]; then
    echo "Error: ${LABEL}: jet file not found: ${JET_FILE}"
    exit 1
  fi

  # Event ranges. Mixing chooses the partners of an event among the events of the range (partner rule), so with mixing no_bullet takes the
  # range of the bullet sample, which keeps the partners of each pair the same. Without mixing each sample runs to its own last event: an extra
  # no_bullet event enters the no_bullet curves, and the difference uses the pairs only
  RANGE=$(range_end "$BULLET_DIR" "$JET_FILE") || exit 1
  LAST_BULLET=${RANGE%%|*}
  CAPS_BULLET=${RANGE#*|}
  LAST_NO_BULLET=$LAST_BULLET
  if [ -n "$NO_BULLET_DIR" ]; then
    if [ "$(option_value RUN_PRODUCER_OPTS mixDepth)" != "0" ]; then
      CAPS_NO_BULLET=" (the bullet range, for mixing)"
    else
      RANGE=$(range_end "$NO_BULLET_DIR" "$JET_FILE") || exit 1
      LAST_NO_BULLET=${RANGE%%|*}
      CAPS_NO_BULLET=${RANGE#*|}
    fi
  fi
  RUN_FIRST_EVENT[$LABEL]=$FIRST_EVENT
  RUN_N_EVENTS_BULLET[$LABEL]=$((LAST_BULLET - FIRST_EVENT + 1))
  RUN_N_EVENTS_NO_BULLET[$LABEL]=$((LAST_NO_BULLET - FIRST_EVENT + 1))

  # Seed
  SEED_SPEC=$(seed_spec "$LABEL")
  if [ -z "$SEED_SPEC" ]; then
    RUN_SEED[$LABEL]=$(draw_seed)
  elif [[ "$SEED_SPEC" =~ ^[0-9]+$ ]]; then
    RUN_SEED[$LABEL]=$SEED_SPEC
  else
    [ -z "${SEED_OF_GROUP[$SEED_SPEC]}" ] && SEED_OF_GROUP[$SEED_SPEC]=$(draw_seed)
    RUN_SEED[$LABEL]=${SEED_OF_GROUP[$SEED_SPEC]}
  fi

  echo "  ${LABEL}: events ${FIRST_EVENT}-${LAST_BULLET}${CAPS_BULLET}${NO_BULLET_DIR:+, no_bullet ${FIRST_EVENT}-${LAST_NO_BULLET}${CAPS_NO_BULLET}}, seed ${RUN_SEED[$LABEL]}${SEED_SPEC:+ (seed=${SEED_SPEC})}, overrides: ${RUN_OVERRIDES[$LABEL]:-none}"
  [ "$(option_value RUN_PRODUCER_OPTS trigger)" = "rp" ] && echo "    [NOTE] rp trigger: mixing turned off (all mixed triggers would coincide with phi_t = 0)"
  MISSING=$(missing_events "$BULLET_DIR" "$LAST_BULLET")
  [ -n "$MISSING" ] && echo "    [NOTE] ${BULLET_DIR#"${HYDRO_BASE}"/}: no ${INPUT_NAME} in events ${MISSING}, which the producer skips"
  if [ -n "$NO_BULLET_DIR" ]; then
    MISSING=$(missing_events "$NO_BULLET_DIR" "$LAST_NO_BULLET")
    [ -n "$MISSING" ] && echo "    [NOTE] ${NO_BULLET_DIR#"${HYDRO_BASE}"/}: no ${INPUT_NAME} in events ${MISSING}, which the producer skips"
  fi
}

if [ $SUMMARY_ONLY -eq 1 ]; then
  echo "SummaryPlots only: ${#COMPARE_LINES[@]} comparisons, on the PerConfig files in ${CAMPAIGN_DIR}"
elif [ $POST_PROCESS_ONLY -eq 1 ]; then
  echo "Post-processing ${#SELECTED_RUNS[@]} runs: PerConfig on their producer files in ${CAMPAIGN_DIR}, then SummaryPlots"
else
  echo "Planning ${#SELECTED_RUNS[@]} runs (${RUNS_FILE}, ${DATASETS_FILE})"
  for LABEL in "${SELECTED_RUNS[@]}"; do
    plan_run "$LABEL"
  done
  # A seed group shares its seed within one invocation only: rerunning part of it (--only) gives that part a new seed
  for GROUP in "${!SEED_OF_GROUP[@]}"; do
    N_IN_REGISTRY=0
    N_SELECTED=0
    for LABEL in "${RUN_LABELS[@]}"; do
      [ "$(seed_spec "$LABEL")" = "$GROUP" ] && N_IN_REGISTRY=$((N_IN_REGISTRY + 1))
    done
    for LABEL in "${SELECTED_RUNS[@]}"; do
      [ "$(seed_spec "$LABEL")" = "$GROUP" ] && N_SELECTED=$((N_SELECTED + 1))
    done
    if [ "$N_SELECTED" -lt "$N_IN_REGISTRY" ]; then
      echo "  [NOTE] seed group ${GROUP}: ${N_SELECTED} of its ${N_IN_REGISTRY} runs selected. They get a new seed, which the others of an earlier invocation do not share"
    fi
  done
fi

# ==============================================================================
# COMPILE THE MACROS AHEAD OF TIME
# ==============================================================================
echo ""
echo "Compiling ROOT macros into native executables..."

COMPILE_WARN_FLAGS="-Wall -Wextra -Wpedantic -Wshadow"
# Same as the main analysis chain, without -DNDEBUG: assert() stays active in the macros of this folder
ROOT_AWARE_FLAGS="$(root-config --cflags)"
# -fno-math-errno: no errno bookkeeping for sqrt and friends. -fno-trapping-math: no floating-point traps are relied upon anywhere
OPTIMIZATION_FLAGS="-O3 -march=native -fno-math-errno -flto=auto -fno-trapping-math"
# -lGenVector for the ROOT::Math vectors and boosts: header-only in practice, but linked explicitly in case a ROOT build moves them into the library
ROOT_LIBS="$(root-config --glibs) -lGenVector"

# Compiles macro $1 into executable $2
compile_macro() {
  echo -n "Compiling $(basename "$1")..."
  # shellcheck disable=SC2086 # the flag strings must split into separate words
  g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$2" "$1" $ROOT_LIBS || exit 1
  echo " done!"
}

[ $POST_PROCESS_ONLY -eq 0 ] && [ $SUMMARY_ONLY -eq 0 ] && compile_macro "$PRODUCER_MACRO" "$PRODUCER_EXE"
[ $SUMMARY_ONLY -eq 0 ] && compile_macro "$PERCONFIG_MACRO" "$PERCONFIG_EXE"
[ ${#COMPARE_LINES[@]} -gt 0 ] && compile_macro "$SUMMARY_MACRO" "$SUMMARY_EXE"
echo "Compilation successful!"
echo ""

# ==============================================================================
# HEADER
# ==============================================================================
echo "========================================================"
echo "  runHydroTruthLab.sh"
echo "  Campaign     : ${CAMPAIGN_DIR}"
echo "  Hydro base   : ${HYDRO_BASE}"
echo "  Runs         : ${#SELECTED_RUNS[@]} of ${#RUN_LABELS[@]}, ${JOBS} at a time"
echo "  Producer     : ${PRODUCER_OPTS[*]}"
echo "  PerConfig    : ${PERCONFIG_OPTS[*]}"
echo "  Comparisons  : ${#COMPARE_LINES[@]}"
if [ $SUMMARY_ONLY -eq 1 ]; then
  echo "  Mode         : SUMMARY PLOTS ONLY"
elif [ $POST_PROCESS_ONLY -eq 1 ]; then
  echo "  Mode         : POST-PROCESS ONLY"
else
  echo "  Mode         : FULL CHAIN"
fi
echo "========================================================"
echo ""

# ==============================================================================
# RUNS
# ==============================================================================
# Progress is printed as one complete line per step, prefixed with the run, never with "echo -n" plus a later "-> OK": partial lines from
# concurrent runs would interleave.
# The macros open their output file at the start, so a failed or interrupted run leaves a partial file behind. They therefore write to
# "<file>.part", renamed only on success, and the file of an earlier invocation is removed first: PerConfig and later -p runs can only ever see
# complete files.

# Arguments: $1 RUN, $2 SAMPLE KEY (bullet or noBullet), $3 SAMPLE FOLDER, $4 NUMBER OF EVENTS. Uses RUN_PRODUCER_OPTS of the calling run
run_producer() {
  local LABEL="$1" KEY="$2" SAMPLE_DIR="$3" N_EVENTS="$4"
  local OUT="${CAMPAIGN_DIR}/${LABEL}/hydroTruthLab_${LABEL}_${KEY}.root"
  local PART="${OUT}.part"
  local LOG="${OUT%.root}.log"
  rm -f "$OUT" "$PART"
  # </dev/null: no child may read this script's stdin
  if "$PRODUCER_EXE" "$SAMPLE_DIR" "$PART" "${RUN_PRODUCER_OPTS[@]}" --jetFile="$(hydro_path "${DATASET_JET_FILE[${RUN_DATASET[$LABEL]}]}")" --inputName="${INPUT_NAME}" --firstEvent="${RUN_FIRST_EVENT[$LABEL]}" --nEvents="${N_EVENTS}" --seed="${RUN_SEED[$LABEL]}" > "$LOG" 2>&1 < /dev/null && [ -f "$PART" ]; then
    mv "$PART" "$OUT"
    # The executable only ever saw the temporary name: point its log at the final file. A failed run keeps ".part", which is what it wrote
    sed -i 's/\.root\.part/.root/g' "$LOG"
    echo "  [${LABEL}] producer ${KEY} -> OK  (warning lines in its log: $(grep -c '^\[WARNING\]' "$LOG"))"
  else
    rm -f "$PART"
    echo "  [${LABEL}] producer ${KEY} -> FAILED  (log: ${LOG})"
    record_failure "${LABEL} | producer ${KEY}"
  fi
}

# Arguments: $1 RUN. Uses RUN_PERCONFIG_OPTS of the calling run
run_perconfig() {
  local LABEL="$1"
  local DATASET="${RUN_DATASET[$LABEL]}"
  local RUN_FOLDER="${CAMPAIGN_DIR}/${LABEL}"
  local BULLET_FILE="${RUN_FOLDER}/hydroTruthLab_${LABEL}_bullet.root"
  local NO_BULLET_ARG="${RUN_FOLDER}/hydroTruthLab_${LABEL}_noBullet.root"
  local OUT="${RUN_FOLDER}/hydroTruthLabPerConfig_${LABEL}.root"
  local PART="${OUT}.part"
  local LOG="${OUT%.root}.log"
  if [ ! -f "$BULLET_FILE" ]; then
    echo "  [${LABEL}] PerConfig -> SKIPPED (no bullet producer file)"
    record_skip "${LABEL} | PerConfig (missing bullet producer file)"
    return
  fi
  # A dataset without no_bullet runs on the bullet sample alone. With one, a missing no_bullet file is reported, and PerConfig still runs
  if [ "${DATASET_NO_BULLET[$DATASET]}" = "none" ]; then
    NO_BULLET_ARG="none"
  elif [ ! -f "$NO_BULLET_ARG" ]; then
    echo "  [${LABEL}] PerConfig: no no_bullet producer file, running on the bullet sample alone (no Difference/)"
    record_skip "${LABEL} | PerConfig Difference/ (missing no_bullet producer file)"
    NO_BULLET_ARG="none"
  fi
  rm -f "$OUT" "$PART"
  if "$PERCONFIG_EXE" "$BULLET_FILE" "$NO_BULLET_ARG" "$PART" "${RUN_PERCONFIG_OPTS[@]}" --centralityLabel="${DATASET_CENTRALITY[$DATASET]}" --initialConditions="${DATASET_INITIAL_CONDITIONS[$DATASET]}" --runLabel="${LABEL}" > "$LOG" 2>&1 < /dev/null && [ -f "$PART" ]; then
    mv "$PART" "$OUT"
    sed -i 's/\.root\.part/.root/g' "$LOG" # as for the producers
    echo "  [${LABEL}] PerConfig -> OK  (warning lines in its log: $(grep -c '^\[WARNING\]' "$LOG"))"
  else
    rm -f "$PART"
    echo "  [${LABEL}] PerConfig -> FAILED  (log: ${LOG})"
    record_failure "${LABEL} | PerConfig"
  fi
}

# One run: its producers at once, then PerConfig. Runs in a background subshell, so its option arrays are its own
run_one() {
  local LABEL="$1"
  local DATASET="${RUN_DATASET[$LABEL]}"
  resolve_run_options "$LABEL"
  if [ $POST_PROCESS_ONLY -eq 0 ]; then
    mkdir -p "${CAMPAIGN_DIR}/${LABEL}" || { record_failure "${LABEL} | run folder"; return; }
    local PRODUCER_PIDS=()
    run_producer "$LABEL" bullet "$(hydro_path "${DATASET_BULLET[$DATASET]}")" "${RUN_N_EVENTS_BULLET[$LABEL]}" &
    PRODUCER_PIDS+=($!)
    if [ "${DATASET_NO_BULLET[$DATASET]}" != "none" ]; then
      run_producer "$LABEL" noBullet "$(hydro_path "${DATASET_NO_BULLET[$DATASET]}")" "${RUN_N_EVENTS_NO_BULLET[$LABEL]}" &
      PRODUCER_PIDS+=($!)
    fi
    wait "${PRODUCER_PIDS[@]}"
  fi
  run_perconfig "$LABEL"
}

if [ $SUMMARY_ONLY -eq 0 ]; then
  for LABEL in "${SELECTED_RUNS[@]}"; do
    # At most JOBS runs at once: wait for one to finish before starting the next
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
      wait -n
    done
    run_one "$LABEL" &
    LANE_PIDS+=($!)
  done
  wait "${LANE_PIDS[@]}"
  LANE_PIDS=()
  echo ""
fi

# ==============================================================================
# SUMMARY PLOTS
# ==============================================================================
# One call for every compare line of the registry, with the PerConfig file of each run of the line, whichever invocation made it (--only
# included). The macro leaves out the runs without one, with a warning
run_summary_plots() {
  local OUT="${CAMPAIGN_DIR}/hydroTruthLabSummaryPlots.root"
  local PART="${OUT}.part"
  local LOG="${OUT%.root}.log"
  local ENTRY SPEC LABEL FIELDS=() COMPARE_ARGS=()
  for ENTRY in "${COMPARE_LINES[@]}"; do
    read -r -a FIELDS <<< "${ENTRY#*|}"
    SPEC="${FIELDS[1]},${FIELDS[2]}"
    for LABEL in "${FIELDS[@]:3}"; do
      SPEC="${SPEC},${CAMPAIGN_DIR}/${LABEL}/hydroTruthLabPerConfig_${LABEL}.root"
    done
    COMPARE_ARGS+=("--compare=${SPEC}")
  done
  rm -f "$OUT" "$PART"
  if "$SUMMARY_EXE" "$PART" "${COMPARE_ARGS[@]}" > "$LOG" 2>&1 < /dev/null && [ -f "$PART" ]; then
    mv "$PART" "$OUT"
    sed -i 's/\.root\.part/.root/g' "$LOG" # as for the producers
    echo "  SummaryPlots -> OK  (warning lines in its log: $(grep -c '^\[WARNING\]' "$LOG"))"
  else
    rm -f "$PART"
    echo "  SummaryPlots -> FAILED  (log: ${LOG})"
    record_failure "comparisons | SummaryPlots"
  fi
}

if [ ${#COMPARE_LINES[@]} -eq 0 ]; then
  echo "  [NOTE] no compare line in ${RUNS_FILE}, so no SummaryPlots"
else
  echo "Running SummaryPlots on ${#COMPARE_LINES[@]} comparisons..."
  run_summary_plots
fi
echo ""

# ==============================================================================
# FAILURE SUMMARY TABLE
# ==============================================================================
echo "========================================================"
mapfile -t FAILURES < "$FAILURES_FILE"
mapfile -t SKIPPED  < "$SKIPPED_FILE"
if [ ${#FAILURES[@]} -eq 0 ]; then
  echo "  All steps completed successfully."
else
  echo "  FAILURES (${#FAILURES[@]} total)"
  echo ""
  print_report_table "FAILED STEP" "${FAILURES[@]}"
  echo "  (logs: next to each output, in ${CAMPAIGN_DIR}/<run>/ and, for SummaryPlots, in ${CAMPAIGN_DIR}/)"
fi
# Printed separately and only when non-empty, so it never competes with the failure table above
if [ ${#SKIPPED[@]} -gt 0 ]; then
  echo "--------------------------------------------------------"
  echo "  SKIPPED (${#SKIPPED[@]} total) -- nothing to process, not an error per se."
  echo ""
  print_report_table "SKIPPED STEP" "${SKIPPED[@]}"
fi
echo "========================================================"

# A nonzero exit status when anything failed, so that the run can be chained or checked by other scripts
[ ${#FAILURES[@]} -eq 0 ] || exit 1