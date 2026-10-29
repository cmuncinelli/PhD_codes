#!/usr/bin/env bash

###############################################################################
# Script name:
#   runHydroTruthLab.sh
#
# Purpose:
#   Runs the hydroTruthLab chain on Vitor's lumpy 40-50% hydro events:
#   (you can generalize to other datasets: just make sure to grab the appropriate inputs!)
#     1. hydroTruthLab.cxx (one process per sample) -> results_producer/hydroTruthLab_<sample>.root
#     2. hydroTruthLabSummary.cxx (both samples)    -> results_summary/hydroTruthLabSummary.root
#
#   Every configurable of both macros lives in the CONFIGURATION block below, including the output folder.
#   Each sample uses every event it has (consecutive event folders from FIRST_EVENT), except that the jet trigger caps both at the last
#   event of the jet file. The summary pairs the two samples over the events they share.
#   For design choices, some definitions and the output reference, refer to the README.md in this folder.
#
# Usage:
#   ./runHydroTruthLab.sh [OPTIONS]
#
# Options:
#   -h, --help                 Show this help message.
#   -p, --post-process-only    Skip Step 1 (producers) and rerun Step 2 (summary) on the existing producer files.
#   -t, --tag NAME             Name of the run subfolder under OUT_BASE (default: "default").
#       --serial               Run the two producers one after the other (default: both at once).
#       --cleanup              Remove the compiled executables and exit. g++ rebuilds them on every run anyway,
#                              so a stale binary can never be picked up.
#
# Output (inside ${OUT_BASE}/${RUN_TAG}/):
#   logHydroTruthLab.log            this script's log (appended, one block per run)
#   backupOf_runHydroTruthLab.sh    copy of this script as it was run: the complete configuration of the run.
#                                   A record only: it is not meant to be run from there (the macros are found next to the original).
#                                   A -p run overwrites it, but the producer files keep their own options in Config/
#   *.root.part                     never left behind: outputs are written under this name and renamed only on success
#   results_producer/               producer files and their logs
#   results_summary/                summary file and its log
#
# Examples:
#   ./runHydroTruthLab.sh
#   ./runHydroTruthLab.sh --tag smudge0p05
#   ./runHydroTruthLab.sh -p --tag smudge0p05
#
###############################################################################

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# The macros sit next to this script, so the repository path is never hardcoded
# (if needed though, it should be "PhD_codes/vorticity_hydro/hydroTruthLab/" unless some structural changes are made to hydroTruthLab)
LAB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRODUCER_MACRO="${LAB_DIR}/hydroTruthLab.cxx"
SUMMARY_MACRO="${LAB_DIR}/hydroTruthLabSummary.cxx"
PRODUCER_EXE="${LAB_DIR}/hydroTruthLab.exe"
SUMMARY_EXE="${LAB_DIR}/hydroTruthLabSummary.exe"

# Output base folder: the only place where it is set
OUT_BASE="${HOME}/RingPol/hydroTruthLab"

# Hydro inputs (Vitor's events in jarvis4). Bullet event i and no_bullet event i share the same initial condition (README, key hypothesis)
HYDRO_BASE="/storage1/vribeiro/lumpy_events_40_50"
BULLET_SAMPLE="bullet"
NO_BULLET_SAMPLE="no_bullet"
SAMPLES=("$BULLET_SAMPLE" "$NO_BULLET_SAMPLE")
INPUT_NAME="Smu_dpTdphidy_Thermal_rapidity_3122.dat"
# The copy in no_bullet/ is identical, so this one file serves both samples: no_bullet event i takes the phi_J of bullet event i
JET_FILE="${HYDRO_BASE}/random_parameters.dat"

# jet: phi_t = phi_J (default). rp: phi_t = 0, along the impact parameter (README, "Trigger modes")
TRIGGER="jet"

# Event counts are found per sample (today: 250 bullet, 300 no_bullet), not set here: every consecutive event folder from FIRST_EVENT that
# holds INPUT_NAME. The jet trigger caps them at the last event of the jet file (249), since later no_bullet events have no phi_J
FIRST_EVENT=0
MAX_EVENTS=0 # 0 = every event found. A positive value caps each sample (quick tests)

# Mixing partners per event. Mixing only exists for the jet trigger (with phi_t = 0 for every event, all mixed triggers coincide), so
# TRIGGER=rp turns it off, with a note in the log
MIX_DEPTH=10

# One seed per run, from the current time, shared by both samples so that their smudged triggers stay paired (README, key hypothesis).
# It is printed in the log and stored in Config/ of each producer file: set a fixed value here to repeat a run exactly
SEED=$(( $(date +%s%N) % 4294967296 )) # 2^32 (current UTC time is of order 2^30~2^31)

# Constants copied from the consumer (lambdaJetPolarizationIonsDerived.cxx): LambdaMass = o2::constants::physics::MassLambda,
# LambdaWeakDecayConstant (DPG 2025 update). Keep them in sync with it by hand: this chain cannot read the O2 headers.
LAMBDA_MASS=1.115683
ALPHA_LAMBDA=0.749

# Producer options shared by both samples (hydroTruthLab.exe --help lists them all). --nEvents and --mixDepth are added per sample in Step 1
PRODUCER_OPTS=(
  --jetFile="${JET_FILE}"
  --trigger="${TRIGGER}"
  --firstEvent="${FIRST_EVENT}"
  --inputName="${INPUT_NAME}"
  # Cuts of Vitor's PRC paper (PhysRevC.109.014905), repeated as a closure test of our code
  --ptMin=0.5
  --ptMax=1.5
  --absYMax=0.5
  --lambdaMass="${LAMBDA_MASS}"
  --deltaScanSteps=100 # one step per phi cell (3.6 deg)
  --smudgeSigma=0.1 # Higher than the consumer's forceJetDirectionSmudge width, which would be 5% * jetR ~ 0.02, with jetR = 0.4
  --smudgeMode=cone
  --smudgeDraws=10
  --seed="${SEED}"
  --mapSubdivisions=3 # PolarizationMaps/: each cell spread over 3x3x3 points of its momentum volume (1 = nodes only, with empty bins at high p)
)

# Summary options (hydroTruthLabSummary.exe --help lists them all)
SUMMARY_OPTS=(
  --paired=1 # key hypothesis: paired initial conditions (README)
  --alphaLambda="${ALPHA_LAMBDA}"
  --nLambdaPrecision=0 # 0: Precision/ holds sigma x sqrt(N) only. N > 0 adds sigma and the expected significance for N Lambdas
  --centralityLabel="40-50%" # label of the Fig. 11 point: the centrality class of HYDRO_BASE
)

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================
POST_PROCESS_ONLY=0
SERIAL=0
CLEANUP_ONLY=0
RUN_TAG="default"

# Prints the header block from "# Usage:" down to the closing line of hashes.
# The end marker is the full hash line, not "^#+$": a bare "#" separator line also matches that pattern and would cut the help short.
print_help() {
  awk '/^# Usage:/{show=1} show && /^##########/{exit} show' "$0" | sed 's/^# \{0,1\}//'
  exit 0
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
    --tag|-t)
      if [ -z "$2" ]; then
        echo "Error: $1 needs a name. Use --help for usage." >&2
        exit 1
      fi
      RUN_TAG="$2"
      shift 2
      ;;
    --serial)
      SERIAL=1
      shift
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

if [ $CLEANUP_ONLY -eq 1 ]; then
  echo "Removing compiled macro executables..."
  rm -f "$PRODUCER_EXE" "$SUMMARY_EXE"
  echo "Done."
  exit 0
fi

# ==============================================================================
# RUN FOLDER AND LOGGING
# ==============================================================================
RUN_DIR="${OUT_BASE}/${RUN_TAG}"
PRODUCER_DIR="${RUN_DIR}/results_producer"
SUMMARY_DIR="${RUN_DIR}/results_summary"
mkdir -p "${PRODUCER_DIR}/logs" "${SUMMARY_DIR}/logs" || exit 1

# Appended rather than overwritten: a post-process-only rerun keeps the log of the run that made the producer files
exec > >(tee -a -i "${RUN_DIR}/logHydroTruthLab.log") 2>&1
echo ""
echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"

cp "${BASH_SOURCE[0]}" "${RUN_DIR}/backupOf_runHydroTruthLab.sh"

# ==============================================================================
# SHARED RUN STATE
# ==============================================================================
# Failure and skip logs live in files, not bash arrays: the producers run as background subshells, which cannot append to the parent's arrays.
# Each entry is one short line written with ">>", which the kernel appends atomically. Format: "SAMPLE | STEP".
RUN_TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/runHydroTruthLab.XXXXXX")
FAILURES_FILE="${RUN_TMP_DIR}/failures.txt"
SKIPPED_FILE="${RUN_TMP_DIR}/skipped.txt"
: > "$FAILURES_FILE"
: > "$SKIPPED_FILE"

PRODUCER_PIDS=()

record_failure() { echo "$1" >> "$FAILURES_FILE"; }
record_skip()    { echo "$1" >> "$SKIPPED_FILE"; }

# Prints one "SAMPLE | STEP" table.
#   $1        : heading for the second column ("FAILED STEP", "SKIPPED STEP")
#   $2 .. $n  : the entries
print_report_table() {
  local STEP_HEADER="$1"
  shift
  printf "  %-20s  %-60s\n" "SAMPLE" "$STEP_HEADER"
  printf "  %-20s  %-60s\n" "--------------------" "------------------------------------------------------------"
  local ENTRY SAMPLE_COL STEP_COL
  for ENTRY in "$@"; do
    SAMPLE_COL=$(echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$1); print $1}')
    STEP_COL=$(  echo "$ENTRY" | awk -F'|' '{gsub(/^ +| +$/,"",$2); print $2}')
    printf "  %-20s  %-60s\n" "$SAMPLE_COL" "$STEP_COL"
  done
}

# ==============================================================================
# SIGNAL HANDLING
# ==============================================================================
# Background jobs of a non-interactive script ignore SIGINT, so Ctrl+C does not reach the producers by itself: they are SIGTERMed explicitly.
handle_interrupt() {
  echo ""
  echo "!!! INTERRUPT DETECTED (Ctrl+C) !!!"
  if [ ${#PRODUCER_PIDS[@]} -gt 0 ]; then
    echo "    Stopping the producers..."
    # Each producer runs inside a subshell (run_producer &): its children are signalled too, so the executables stop, not only the subshells
    pkill -TERM -P "$(echo "${PRODUCER_PIDS[@]}" | tr ' ' ',')" 2>/dev/null
    kill -TERM "${PRODUCER_PIDS[@]}" 2>/dev/null
    wait "${PRODUCER_PIDS[@]}" 2>/dev/null
  fi
  # Whatever was running when the signal came left only "<file>.part" behind: never a final file
  rm -f "${PRODUCER_DIR}"/*.root.part "${SUMMARY_DIR}"/*.root.part
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
# PRE-FLIGHT CHECKS
# ==============================================================================
if [ ! -f "$PRODUCER_MACRO" ]; then
  echo "Error: producer macro not found: ${PRODUCER_MACRO}"
  exit 1
fi

# The summary macro is optional while it does not exist yet: Step 2 is then skipped, and reported as such
HAS_SUMMARY=1
if [ ! -f "$SUMMARY_MACRO" ]; then
  HAS_SUMMARY=0
  echo "[WARNING] summary macro not found (${SUMMARY_MACRO}): Step 2 will be skipped."
fi

if [ "$TRIGGER" != "jet" ] && [ "$TRIGGER" != "rp" ]; then
  echo "Error: TRIGGER must be jet or rp, got '${TRIGGER}'"
  exit 1
fi

# Per-sample event counts and mixing depth, passed in Step 1
declare -A N_EVENTS_FOR
EFFECTIVE_MIX_DEPTH=$MIX_DEPTH
if [ $POST_PROCESS_ONLY -eq 0 ]; then
  if [ ! -f "$JET_FILE" ]; then
    echo "Error: jet file not found: ${JET_FILE}"
    exit 1
  fi
  if [ "$TRIGGER" = "rp" ] && [ "$MIX_DEPTH" -gt 0 ]; then
    echo "[NOTE] TRIGGER=rp: mixing turned off (all mixed triggers would coincide with phi_t = 0)"
    EFFECTIVE_MIX_DEPTH=0
  fi
  # Last event of the jet file (the producer checks every entry and the atan2 consistency anyway)
  LAST_JET_EVENT=$(awk '!/^[[:space:]]*#/ && NF {print int($1)}' "$JET_FILE" | sort -n | tail -1)

  for SAMPLE in "${SAMPLES[@]}"; do
    # Consecutive event folders from FIRST_EVENT that hold the input file
    N_FOUND=0
    while [ -f "${HYDRO_BASE}/${SAMPLE}/$((FIRST_EVENT + N_FOUND))/results/${INPUT_NAME}" ]; do
      N_FOUND=$((N_FOUND + 1))
    done
    if [ $N_FOUND -eq 0 ]; then
      echo "Error: no event file found at ${HYDRO_BASE}/${SAMPLE}/${FIRST_EVENT}/results/${INPUT_NAME}"
      exit 1
    fi
    # Numbered folders beyond a gap would be dropped without a word otherwise
    N_FOLDERS=$(find "${HYDRO_BASE}/${SAMPLE}" -mindepth 1 -maxdepth 1 -type d -regex '.*/[0-9]+' | wc -l)
    if [ "$N_FOLDERS" -gt $((FIRST_EVENT + N_FOUND)) ]; then
      echo "[WARNING] ${SAMPLE}: ${N_FOLDERS} numbered event folders exist, but only events ${FIRST_EVENT} to $((FIRST_EVENT + N_FOUND - 1)) are consecutive with ${INPUT_NAME}: the rest is not used"
    fi

    N_USED=$N_FOUND
    CAP_NOTE=""
    if [ "$TRIGGER" = "jet" ] && [ $N_USED -gt $((LAST_JET_EVENT - FIRST_EVENT + 1)) ]; then
      N_USED=$((LAST_JET_EVENT - FIRST_EVENT + 1))
      CAP_NOTE=" (capped by the jet file, whose last event is ${LAST_JET_EVENT})"
    fi
    if [ "$MAX_EVENTS" -gt 0 ] && [ $N_USED -gt "$MAX_EVENTS" ]; then
      N_USED=$MAX_EVENTS
      CAP_NOTE="${CAP_NOTE} (capped by MAX_EVENTS)"
    fi
    echo "  ${SAMPLE}: ${N_FOUND} events found from ${FIRST_EVENT}, ${N_USED} used${CAP_NOTE}"
    N_EVENTS_FOR[$SAMPLE]=$N_USED
  done
fi

# ==============================================================================
# COMPILE THE MACROS AHEAD OF TIME
# ==============================================================================
echo "Compiling ROOT macros into native executables..."

COMPILE_WARN_FLAGS="-Wall -Wextra -Wpedantic -Wshadow"
# Same as the main analysis chain, without -DNDEBUG: assert() stays active in the macros of this folder
ROOT_AWARE_FLAGS="$(root-config --cflags)"
# -fno-math-errno: no errno bookkeeping for sqrt and friends. -fno-trapping-math: no floating-point traps are relied upon anywhere
OPTIMIZATION_FLAGS="-O3 -march=native -fno-math-errno -flto=auto -fno-trapping-math"
# -lGenVector for the ROOT::Math vectors and boosts: header-only in practice, but linked explicitly in case a ROOT build moves them into the library
ROOT_LIBS="$(root-config --glibs) -lGenVector"

if [ $POST_PROCESS_ONLY -eq 0 ]; then
  echo -n "Compiling $(basename "$PRODUCER_MACRO")..."
  # shellcheck disable=SC2086 # the flag strings must split into separate words
  g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$PRODUCER_EXE" "$PRODUCER_MACRO" $ROOT_LIBS || exit 1
  echo " done!"
fi
if [ $HAS_SUMMARY -eq 1 ]; then
  echo -n "Compiling $(basename "$SUMMARY_MACRO")..."
  # shellcheck disable=SC2086 # the flag strings must split into separate words
  g++ $COMPILE_WARN_FLAGS $OPTIMIZATION_FLAGS $ROOT_AWARE_FLAGS -o "$SUMMARY_EXE" "$SUMMARY_MACRO" $ROOT_LIBS || exit 1
  echo " done!"
fi
echo "Compilation successful!"
echo ""

# ==============================================================================
# SUMMARY HEADER
# ==============================================================================
echo "========================================================"
echo "  runHydroTruthLab.sh"
echo "  Run folder   : ${RUN_DIR}"
echo "  Hydro base   : ${HYDRO_BASE}"
for SAMPLE in "${SAMPLES[@]}"; do
  if [ $POST_PROCESS_ONLY -eq 0 ]; then
    printf "  %-13s: %s\n" "$SAMPLE" "events ${FIRST_EVENT} to $((FIRST_EVENT + N_EVENTS_FOR[$SAMPLE] - 1))"
  fi
done
echo "  Trigger      : ${TRIGGER}"
echo "  Mixing depth : ${EFFECTIVE_MIX_DEPTH}"
echo "  Seed         : ${SEED}"
echo "  Jet file     : ${JET_FILE}"
echo "  Producer     : ${PRODUCER_OPTS[*]}"
echo "  Summary      : ${SUMMARY_OPTS[*]}"
echo "  Mode         : $( [ $POST_PROCESS_ONLY -eq 1 ] && echo 'POST-PROCESS ONLY' || echo 'FULL CHAIN' )$( [ $SERIAL -eq 1 ] && echo ' [SERIAL PRODUCERS]' )"
echo "========================================================"
echo ""

# ==============================================================================
# STEP 1: PRODUCERS
# ==============================================================================
# One single-threaded process per sample. Progress is printed as one complete line per producer, never with "echo -n" plus a later "-> OK":
# partial lines from two concurrent producers would interleave.
# The producer opens its output file at the start, so a failed or interrupted run leaves a partial file behind. It therefore writes to
# "<file>.part", renamed only on success, and the file of an earlier run with this tag is removed first: Step 2 and later -p runs can only
# ever see complete files from this run.
# Arguments: $1 SAMPLE
run_producer() {
  local SAMPLE="$1"
  local PRODUCER_OUT="${PRODUCER_DIR}/hydroTruthLab_${SAMPLE}.root"
  local PRODUCER_PART="${PRODUCER_OUT}.part"
  local PRODUCER_LOG="${PRODUCER_DIR}/logs/producer_${SAMPLE}.log"
  rm -f "$PRODUCER_OUT" "$PRODUCER_PART"
  # </dev/null: no child may read this script's stdin
  if "$PRODUCER_EXE" "${HYDRO_BASE}/${SAMPLE}" "$PRODUCER_PART" "${PRODUCER_OPTS[@]}" --nEvents="${N_EVENTS_FOR[$SAMPLE]}" --mixDepth="${EFFECTIVE_MIX_DEPTH}" > "$PRODUCER_LOG" 2>&1 < /dev/null && [ -f "$PRODUCER_PART" ]; then
    mv "$PRODUCER_PART" "$PRODUCER_OUT"
    # The executable only ever saw the temporary name: point its log at the final file. A failed run keeps ".part", which is what it wrote
    sed -i 's/\.root\.part/.root/g' "$PRODUCER_LOG"
    echo "  [1/2] producer : ${SAMPLE}  -> OK  (warning lines in its log: $(grep -c '^\[WARNING\]' "$PRODUCER_LOG"))"
  else
    rm -f "$PRODUCER_PART"
    echo "  [1/2] producer : ${SAMPLE}  -> FAILED  (log: ${PRODUCER_LOG})"
    record_failure "${SAMPLE} | producer"
  fi
}

if [ $POST_PROCESS_ONLY -eq 0 ]; then
  for SAMPLE in "${SAMPLES[@]}"; do
    if [ $SERIAL -eq 1 ]; then
      run_producer "$SAMPLE"
    else
      run_producer "$SAMPLE" &
      PRODUCER_PIDS+=($!)
    fi
  done
  if [ ${#PRODUCER_PIDS[@]} -gt 0 ]; then
    wait "${PRODUCER_PIDS[@]}"
    PRODUCER_PIDS=()
  fi
  echo ""
else
  echo "  [1/2] producers: SKIPPED (post-process only)"
  echo ""
fi

# ==============================================================================
# STEP 2: SUMMARY
# ==============================================================================
BULLET_FILE="${PRODUCER_DIR}/hydroTruthLab_${BULLET_SAMPLE}.root"
NO_BULLET_FILE="${PRODUCER_DIR}/hydroTruthLab_${NO_BULLET_SAMPLE}.root"
SUMMARY_OUT="${SUMMARY_DIR}/hydroTruthLabSummary.root"
SUMMARY_LOG="${SUMMARY_DIR}/logs/summary.log"

if [ $HAS_SUMMARY -eq 0 ]; then
  echo "  [2/2] summary  : SKIPPED (no summary macro yet)"
  record_skip "both | summary (macro not found)"
elif [ ! -f "$BULLET_FILE" ]; then
  echo "  [2/2] summary  : SKIPPED (no bullet producer file: ${BULLET_FILE})"
  record_skip "${BULLET_SAMPLE} | summary (missing producer file)"
else
  # Without the no_bullet file, the summary still runs on the bullet sample alone and skips Difference/
  SUMMARY_NO_BULLET_ARG="$NO_BULLET_FILE"
  if [ ! -f "$NO_BULLET_FILE" ]; then
    echo "  [2/2] summary  : no no_bullet producer file, running on the bullet sample alone (no Difference/)"
    record_skip "${NO_BULLET_SAMPLE} | summary Difference/ (missing producer file)"
    SUMMARY_NO_BULLET_ARG="none"
  fi
  # Same "<file>.part" scheme as the producers
  SUMMARY_PART="${SUMMARY_OUT}.part"
  rm -f "$SUMMARY_OUT" "$SUMMARY_PART"
  if "$SUMMARY_EXE" "$BULLET_FILE" "$SUMMARY_NO_BULLET_ARG" "$SUMMARY_PART" "${SUMMARY_OPTS[@]}" > "$SUMMARY_LOG" 2>&1 < /dev/null && [ -f "$SUMMARY_PART" ]; then
    mv "$SUMMARY_PART" "$SUMMARY_OUT"
    sed -i 's/\.root\.part/.root/g' "$SUMMARY_LOG" # as for the producers
    echo "  [2/2] summary  : -> OK"
  else
    rm -f "$SUMMARY_PART"
    echo "  [2/2] summary  : -> FAILED  (log: ${SUMMARY_LOG})"
    record_failure "both | summary"
  fi
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
  echo "  (logs: ${PRODUCER_DIR}/logs/ and ${SUMMARY_DIR}/logs/)"
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