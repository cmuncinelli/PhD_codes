#!/usr/bin/env bash
# =============================================================================
# runHelicityToyModel.sh
# =============================================================================
#
# PURPOSE
# -------
# Coordinator script for the Lambda helicity efficiency toy model study.
# Runs helicityEfficiencyToyModel.cxx and plotHelicityEfficiency.cxx over
# several families of parameter variations designed to isolate the two
# detector-induced fake-polarization effects:
#
#   (1) The helicity efficiency effect (forward-backward / cos theta* asymmetry)
#       driven by the minimum-pT threshold applied to daughter tracks.
#
#   (2) The DCA-cut left-right asymmetry (phi* modulation) driven by the
#       interplay between the magnetic field direction and the minimum-DCA
#       requirement applied to daughter tracks.
#
# PARALLELISM DESIGN
# ------------------
# All jobs from ALL families are collected into a single queue first, then
# dispatched simultaneously up to MAX_PARALLEL concurrent processes.
# On a 192-core machine with 24 total runs, the default MAX_PARALLEL=24
# launches every job at once and the total wall time equals the slowest
# single run.
#
# Dispatch strategy (in order of preference):
#   1. GNU parallel (parallel --jobs N): preferred; handles job tracking,
#      per-job logging via --joblog, and clean failure reporting.
#   2. Pure bash background pool: fallback if parallel is not installed;
#      uses the same subshell structure with a lightweight semaphore.
#
# The family structure is kept purely for the --family filter convenience:
# it controls which jobs enter the queue, not when they run.
# run_family_header/footer are informational only and contain NO wait calls.
#
# PARAMETER FAMILIES
# ------------------
# The script runs the following scan families, each producing an independent
# ROOT file and a dedicated plots subfolder:
#
#   Family 0: BASELINE
#       One reference run with ALICE O-O default parameters.
#       Every other family uses this as the comparison point.
#
#   Family 1: MAGNETIC FIELD SIGN
#       Bz = +0.5 T  vs  Bz = -0.5 T
#       The DCA-asymmetry sign must flip with the field.  The cos(theta*)
#       helicity asymmetry is field-independent and must NOT flip.
#
#   Family 2: DAUGHTER pT THRESHOLD
#       pT_min = 0.10, 0.15, 0.20, 0.30 GeV/c (applied to both daughters)
#       The helicity efficiency effect grows with the pT threshold because
#       a higher threshold removes a larger fraction of the soft-pion
#       configurations.  For a threshold of 0 (no cut), the cos(theta*)
#       distribution must be flat.  This scan verifies the monotonic
#       dependence and identifies the threshold at which the effect becomes
#       significant relative to the expected physics signal (~10^-3).
#       Reference: ALICE track quality cuts typically require pT > 0.15 GeV/c
#       for TPC tracks; the ITS-only regime extends down to ~0.08 GeV/c.
#
#   Family 3: DCA THRESHOLD
#       DCA_min = 0.0/0.0 (no cut), 0.02/0.05, 0.05/0.10, 0.10/0.20 cm
#       (proton / pion)
#       The left-right phi* asymmetry grows with the DCA threshold because
#       a tighter cut selects more displaced decay geometries where the
#       magnetic bending has a larger differential effect on the two daughters.
#       The 0.0/0.0 baseline must produce a flat phi* distribution.
#       Reference: ALICE V0 topological cuts for Lambda selection typically
#       require DCA(pion to PV) > 0.1 cm and DCA(proton to PV) > 0.05 cm
#
#   Family 4: LAMBDA KINEMATIC WINDOW
#       Window A: pT in [0.3, 10.0] GeV/c, |y| < 0.9  (full acceptance)
#       Window B: pT in [0.5,  1.5] GeV/c, |y| < 0.5  (ring analysis cuts)
#       Window C: pT in [1.5,  4.0] GeV/c, |y| < 0.5  (harder Lambda regime)
#       The ring observable analysis uses the kinematic window [0.5, 1.5] GeV/c
#       and |y| < 0.5 to select the Lambda population most sensitive to
#       vorticity-induced polarization (as established in PhysRevC.109.014905).
#       Restricting to this window also changes which part of the daughter pT
#       spectrum is populated, potentially affecting both fake-signal magnitudes.
#
#   Family 5: THERMAL TEMPERATURE (Lambda pT spectrum shape)
#       T = 0.25, 0.30, 0.35, 0.45 GeV
#       The Boltzmann mT temperature controls the Lambda pT spectrum.
#       Values of T ~ 0.25-0.35 GeV are typical for Lambda production in
#       Pb-Pb at LHC energies (measured by ALICE: JHEP 07 (2015) 116).
#       T ~ 0.45 GeV is an upper estimate including harder production channels.
#       Because softer Lambdas produce softer daughters, the effective pT and
#       DCA fake-signal magnitudes are spectrum-dependent.
#
#   Family 6: STATISTICS SCALING
#       N = 100k, 500k, 2M, 5M Lambdas
#       Verifies that the measured <R_proxy> values are statistically stable
#       (do not drift with statistics) and confirms the expected 1/sqrt(N)
#       scaling of the error bars.  Also establishes the minimum N needed
#       for a statistically significant measurement of the fake signal size,
#       which informs how many events are needed in the data analysis.
#
# USAGE
# -----
#   chmod +x runHelicityToyModel.sh
#   ./runHelicityToyModel.sh                    # Run all families (all at once)
#   ./runHelicityToyModel.sh --family 1         # Register and run only Family 1
#   ./runHelicityToyModel.sh --family 2 3       # Families 2 and 3 only
#   ./runHelicityToyModel.sh --jobs 8           # Override parallelism limit
#   ./runHelicityToyModel.sh --dry-run          # Print jobs without executing
#   ./runHelicityToyModel.sh --list             # List all registered jobs and exit
#   ./runHelicityToyModel.sh --help             # Show this help
#
# OUTPUT STRUCTURE
# ----------------
#   BASE_DIR/
#     baseline/        helicity_baseline.root   plots/
#     field_pos/       field_neg/
#     ptcut_000/ ... ptcut_030/
#     dcacut_none/ ... dcacut_tight/
#     window_full/     window_ring/     window_hard/
#     temp_025/ ... temp_045/
#     stats_100k/ ... stats_50M/
#     logs/
#       <run_name>.log           (one per job: generator + plotter combined)
#       parallel_joblog.tsv      (GNU parallel timing/status log, if used)
#
# REQUIREMENTS
# ------------
#   - ROOT 6.14+ with C++14 support (root must be in PATH)
#   - helicityEfficiencyToyModel.cxx and plotHelicityEfficiency.cxx in the
#     same directory as this script
#   - GNU parallel (optional but recommended; apt install parallel)
#   - Bash 4.0+
#
# DISK ESTIMATE
# -------------
#   Each ROOT file (histograms only, no TTrees): ~2-4 MB after ROOT zlib
#   compression (independent of N; all data is binned into fixed histograms).
#   Each plots/ folder: ~5-15 MB of PDFs (13 files).
#   Total for 24 runs: ~200-450 MB.
#   Log files: ~50-200 KB per run; negligible.
#
# =============================================================================

set -euo pipefail   # Exit on error, undefined variable, or pipe failure

# =============================================================================
# CONFIGURATION -- edit these paths as needed
# =============================================================================

# Directory containing this script and the .cxx macros
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base output directory for all ROOT files and plot folders
BASE_DIR="/home/users/cicerodm/RingPol/HelicityToyModel"

# Log directory (on SSD: no buffering concerns; each log is written
# sequentially within its own subshell, so concurrent writes go to
# separate files and there is no interleaving risk)
LOG_DIR="${BASE_DIR}/logs"

# Maximum concurrent ROOT processes.
# Default: 30 (all jobs at once on a 192-core machine with 24 total runs).
# Set to 0 to let GNU parallel use all available cores (parallel -j 0),
# or to any positive integer to cap concurrency explicitly.
MAX_PARALLEL=30

# Default Lambda count per run.
# At ~90 s per 10M Lambdas, 1B = ~900 s (~15 min) per run.
# With 24 runs in parallel on 192 cores, total wall time ~ 15 min.
DEFAULT_N=1000000000  # Bumped from 1.000.000 to 100.000.000 (10.000.000 takes about 1m30s per run, and we have only 23 runs)

# ROOT executable (set to full path if not in PATH, e.g. /opt/root/bin/root)
ROOT_EXE="root"

# Generator macro and plotting macro (relative to SCRIPT_DIR)
GEN_MACRO="${SCRIPT_DIR}/helicityEfficiencyToyModel.cxx"
PLT_MACRO="${SCRIPT_DIR}/plotHelicityEfficiency.cxx"


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

DRY_RUN=0
LIST_ONLY=0
FAMILIES_TO_RUN=()

print_help() {
    # Extract the USAGE block from this file's header comment
    awk '/^# USAGE/,/^# OUTPUT/' "$0" | sed 's/^# *//' | head -n -1
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            print_help
            ;;
        --dry-run|-n)
            DRY_RUN=1
            shift
            ;;
        --list|-l)
            LIST_ONLY=1
            shift
            ;;
        --jobs|-j)
            if [[ -z "${2:-}" || ! "${2}" =~ ^[0-9]+$ ]]; then
                echo "ERROR: --jobs requires a non-negative integer argument." >&2
                exit 1
            fi
            MAX_PARALLEL="$2"
            shift 2
            ;;
        --family|-f)
            shift
            while [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; do
                FAMILIES_TO_RUN+=("$1")
                shift
            done
            ;;
        *)
            echo "Unknown argument: $1.  Use --help for usage." >&2
            exit 1
            ;;
    esac
done

# Default: run all families
if [[ ${#FAMILIES_TO_RUN[@]} -eq 0 ]]; then
    FAMILIES_TO_RUN=(0 1 2 3 4 5 6)
fi


# =============================================================================
# SETUP
# =============================================================================

# Create base and log directories
mkdir -p "${BASE_DIR}"
mkdir -p "${LOG_DIR}"

# Verify macros exist (skip in list/dry-run mode so the user can check
# the job table even before copying the .cxx files into place)
if [[ ${DRY_RUN} -eq 0 && ${LIST_ONLY} -eq 0 ]]; then
    for MACRO in "${GEN_MACRO}" "${PLT_MACRO}"; do
        if [[ ! -f "${MACRO}" ]]; then
            echo "ERROR: Macro not found: ${MACRO}" >&2
            echo "       Place helicityEfficiencyToyModel.cxx and" >&2
            echo "       plotHelicityEfficiency.cxx in: ${SCRIPT_DIR}" >&2
            exit 1
        fi
    done
fi

# Check for GNU parallel availability once upfront
USE_GNU_PARALLEL=0
if command -v parallel &>/dev/null; then
    USE_GNU_PARALLEL=1
fi


# =============================================================================
# JOB QUEUE
# Each call to register_job() appends one record to JOB_QUEUE[].
# Format: colon-delimited, 13 fields matching helicityEfficiencyToyModel args:
#   NAME:SUBDIR:N:BZ:PTMIN_LAM:PTMAX_LAM:RAPMAX:T:PTMIN_P:PTMIN_PI:DCAMIN_P:DCAMIN_PI:SEED
# No spaces in any field.
# =============================================================================

declare -a JOB_QUEUE=()
declare -a JOB_FAMILIES=()   # Parallel array: family number for each job

# family_in_scope: returns 0 (true) if family NUM is in the run list
family_in_scope() {
    local NUM="$1"
    local FAM
    for FAM in "${FAMILIES_TO_RUN[@]}"; do
        [[ "${FAM}" == "${NUM}" ]] && return 0
    done
    return 1
}

# register_job: same 13-argument signature as the old run_one,
# but only adds to the queue -- does NOT launch anything.
#
#   $1  NAME       short identifier, used for output subdir and log name
#   $2  SUBDIR     subdirectory name under BASE_DIR  (usually same as NAME)
#   $3  N          number of Lambdas
#   $4  BZ         magnetic field [T]
#   $5  PTMIN_LAM  Lambda min pT [GeV/c]
#   $6  PTMAX_LAM  Lambda max pT [GeV/c]
#   $7  RAPMAX     Lambda max |rapidity|
#   $8  T          Boltzmann temperature [GeV]
#   $9  PTMIN_P    proton min pT [GeV/c]
#   $10 PTMIN_PI   pion min pT [GeV/c]
#   $11 DCAMIN_P   proton min DCA_xy [cm]
#   $12 DCAMIN_PI  pion min DCA_xy [cm]
#   $13 SEED       TRandom3 seed
register_job() {
    local PACKED="$1:$2:$3:$4:$5:$6:$7:$8:$9:${10}:${11}:${12}:${13}"
    JOB_QUEUE+=("${PACKED}")
    JOB_FAMILIES+=("${CURRENT_FAMILY}")
}

# run_family_header: cosmetic only -- no wait, no throttle.
# Records the current family number so register_job can tag each entry.
CURRENT_FAMILY=0
run_family_header() {
    CURRENT_FAMILY="$1"
    echo ""
    echo "------------------------------------------------------------"
    echo "  Registering Family ${1}: ${2}"
    echo "------------------------------------------------------------"
}


# =============================================================================
# JOB EXECUTOR (called once per job, in a subshell by the dispatcher)
# =============================================================================
# This function is exported so GNU parallel can call it in a child process.
# It receives one colon-delimited packed argument string.

execute_job() {
    # Unpack the 13 fields
    IFS=':' read -r \
        NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX T \
        PTMIN_P PTMIN_PI DCAMIN_P DCAMIN_PI SEED \
        <<< "$1"

    local OUTDIR="${BASE_DIR}/${SUBDIR}"
    local ROOT_FILE="${OUTDIR}/helicity_${NAME}.root"
    local PLOTS_DIR="${OUTDIR}/plots"
    local LOG_FILE="${LOG_DIR}/${NAME}.log"

    mkdir -p "${OUTDIR}"
    mkdir -p "${PLOTS_DIR}"

    # Build ROOT call strings.
    # Single quotes around the macro path protect spaces; the inner double
    # quotes around the ROOT file path protect spaces in BASE_DIR.
    local GEN_CALL="${N},\"${ROOT_FILE}\",${BZ},${PTMIN_LAM},${PTMAX_LAM},${RAPMAX},${T},${PTMIN_P},${PTMIN_PI},${DCAMIN_P},${DCAMIN_PI},${SEED}"
    local PLT_CALL="\"${ROOT_FILE}\",\"${PLOTS_DIR}\""

    # Write log header
    {
        echo "=== ${NAME} started at $(date '+%Y-%m-%d %H:%M:%S') ==="
        echo "Host      : $(hostname)"
        echo "PID       : $$"
        echo "Generator : ${ROOT_EXE} -l -b -q '${GEN_MACRO}(${GEN_CALL})'"
        echo "Plotter   : ${ROOT_EXE} -l -b -q '${PLT_MACRO}(${PLT_CALL})'"
        echo ""
    } > "${LOG_FILE}"

    # Run generator
    local T0; T0=$(date +%s)
    if eval "${ROOT_EXE} -l -b -q '${GEN_MACRO}(${GEN_CALL})'" \
            >> "${LOG_FILE}" 2>&1; then
        local T1; T1=$(date +%s)
        echo "Generator OK  ($(( T1 - T0 )) s)" >> "${LOG_FILE}"
    else
        local STATUS=$?
        echo "Generator FAILED (exit ${STATUS})" >> "${LOG_FILE}"
        echo "ERROR [${NAME}]: generator failed (exit ${STATUS}). Log: ${LOG_FILE}" >&2
        return ${STATUS}
    fi

    echo "" >> "${LOG_FILE}"

    # Run plotter
    local T2; T2=$(date +%s)
    if eval "${ROOT_EXE} -l -b -q '${PLT_MACRO}(${PLT_CALL})'" \
            >> "${LOG_FILE}" 2>&1; then
        local T3; T3=$(date +%s)
        echo "Plotter  OK  ($(( T3 - T2 )) s)" >> "${LOG_FILE}"
    else
        local STATUS=$?
        echo "Plotter FAILED (exit ${STATUS})" >> "${LOG_FILE}"
        echo "ERROR [${NAME}]: plotter failed (exit ${STATUS}). Log: ${LOG_FILE}" >&2
        return ${STATUS}
    fi

    local TOTAL=$(( $(date +%s) - T0 ))
    echo "" >> "${LOG_FILE}"
    echo "=== ${NAME} finished at $(date '+%Y-%m-%d %H:%M:%S') (total: ${TOTAL} s) ===" \
        >> "${LOG_FILE}"

    # Brief status line to stdout (visible in the parent terminal)
    printf "  [DONE] %-28s  %d s\n" "${NAME}" "${TOTAL}"
}

# Export so GNU parallel can call it across forked subshells
export -f execute_job
# Also export the path variables that execute_job needs
export BASE_DIR LOG_DIR ROOT_EXE GEN_MACRO PLT_MACRO


# =============================================================================
# FAMILY DEFINITIONS
# Each block is guarded by family_in_scope(); it only calls register_job.
# The order of register_job calls within a block is irrelevant because all
# jobs are dispatched simultaneously by the pool at the end.
# =============================================================================

# -----------------------------------------------------------------------------
# FAMILY 0: BASELINE
# Default ALICE O-O parameters.  All other families should be compared
# against this reference run.
# Parameters: N=1M, Bz=+0.5T, pT_lam=[0.3,10], |y|<0.9, T=0.30 GeV,
#             pT_min_daughters=0.15 GeV/c, DCA_min=[0.05,0.10] cm
# -----------------------------------------------------------------------------
if family_in_scope 0; then
    run_family_header 0 "BASELINE (ALICE O-O defaults)"

    #              NAME         SUBDIR       N            Bz   pTmin pTmax rap    T     pTp   pTpi  dcaP  dcaPi seed
    register_job  "baseline"  "baseline"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
fi

# -----------------------------------------------------------------------------
# Family 1: MAGNETIC FIELD SIGN
# Swap Bz from +0.5 T to -0.5 T (same as baseline otherwise).
# Expected result: the phi* left-right asymmetry and the eta-antisymmetric
# ring proxy must flip sign exactly.  The cos(theta*) helicity asymmetry
# is field-independent (it only depends on pT thresholds) and must NOT flip.
#
# The ALICE solenoid operated at both polarities during Run 2; in Run 3 the
# nominal polarity is +0.5 T.  The data-driven DCA-based correction must be
# validated separately for each polarity.
# -----------------------------------------------------------------------------
if family_in_scope 1; then
    run_family_header 1 "MAGNETIC FIELD SIGN"

    register_job  "field_pos"  "field_pos"  ${DEFAULT_N}  +0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    register_job  "field_neg"  "field_neg"  ${DEFAULT_N}  -0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
fi

# -----------------------------------------------------------------------------
# Family 2: DAUGHTER pT THRESHOLD SCAN
# Vary pT_min applied to BOTH daughters simultaneously.
# pT = 0.10 GeV/c: lower edge of TPC tracking efficiency in ALICE Run 3
# pT = 0.15 GeV/c: ALICE standard V0 daughter cut (baseline)
# pT = 0.20 GeV/c: tighter cut sometimes used in high-multiplicity Pb-Pb
# pT = 0.30 GeV/c: upper limit explored in systematic studies
#
# A run with pT = 0.00 (no cut) is included as the zero-asymmetry anchor.
# It must produce a flat cos(theta*) distribution.
#
# Reference: ALICE-PUBLIC-2017-005 (V0 selection);
#            ALICE strangeness PWG internal note on Lambda selection
# -----------------------------------------------------------------------------
if family_in_scope 2; then
    run_family_header 2 "DAUGHTER pT THRESHOLD SCAN"

    # 0.001 GeV/c instead of 0.000 avoids a divide-by-zero guard in the helix
    register_job  "ptcut_000"  "ptcut_000"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.001 0.001 0.050 0.100  42
    register_job  "ptcut_010"  "ptcut_010"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.100 0.100 0.050 0.100  42
    register_job  "ptcut_015"  "ptcut_015"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    register_job  "ptcut_020"  "ptcut_020"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.200 0.200 0.050 0.100  42
    register_job  "ptcut_030"  "ptcut_030"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.300 0.300 0.050 0.100  42
fi

# -----------------------------------------------------------------------------
# Family 3: DCA THRESHOLD SCAN
# DCA cuts are specified as (proton DCA min, pion DCA min).
# The pion cut is typically tighter than the proton cut because the pion
# is lighter and therefore subject to larger multiple-scattering DCA smearing.
#
# dcacut_none  [0.00, 0.00]: no DCA requirement -- phi* must be flat
# dcacut_loose [0.05, 0.05]: symmetric DCA cut
# dcacut_loose [0.02, 0.05]: very loose, near the smearing limit
# dcacut_std   [0.05, 0.10]: ALICE standard for Lambda (baseline)
# dcacut_tight [0.10, 0.20]: tighter, sometimes used in central Pb-Pb
# -----------------------------------------------------------------------------
if family_in_scope 3; then
    run_family_header 3 "DCA THRESHOLD SCAN"

    register_job  "dcacut_none"            "dcacut_none"            ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.000 0.000  42
    register_job  "dcacut_loose"           "dcacut_loose"           ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.020 0.050  42
    register_job  "dcacut_loose_symmetric" "dcacut_loose_symmetric" ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.050 0.050  42
    register_job  "dcacut_std"             "dcacut_std"             ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.050 0.100  42
    register_job  "dcacut_tight"           "dcacut_tight"           ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.100 0.200  42
fi

# -----------------------------------------------------------------------------
# Family 4: LAMBDA KINEMATIC WINDOW
# Three kinematic windows motivated by the ring analysis and by the Lambda
# polarization measurement strategy in ALICE.
#
# window_full:  pT in [0.3, 10.0] GeV/c, |y| < 0.9  (full TPC acceptance)
# window_ring:  pT in [0.5,  1.5] GeV/c, |y| < 0.5  (ring analysis window
#               from PhysRevC.109.014905 -- vorticity-sensitive population)
# window_hard:  pT in [1.5,  4.0] GeV/c, |y| < 0.5  (harder fragmentation
#               regime, expected larger fragmentation polarization baseline)
#
# The ring-analysis window is the most physically important: the fake-signal
# size in this specific kinematic range is what must be corrected in the data.
# A smaller fake signal here would reduce the needed correction factor.
# -----------------------------------------------------------------------------
if family_in_scope 4; then
    run_family_header 4 "LAMBDA KINEMATIC WINDOW"

    register_job  "window_full"  "window_full"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    register_job  "window_ring"  "window_ring"  ${DEFAULT_N}  0.5  0.5   1.5  0.5  0.30  0.150 0.150 0.050 0.100  42
    register_job  "window_hard"  "window_hard"  ${DEFAULT_N}  0.5  1.5   4.0  0.5  0.30  0.150 0.150 0.050 0.100  42
fi

# -----------------------------------------------------------------------------
# Family 5: THERMAL TEMPERATURE SCAN
# The Boltzmann mT exponential slope parameter T controls which Lambda pT
# values are populated.  Softer spectra (lower T) produce more soft Lambdas
# whose daughters are more likely to be affected by the pT threshold.
#
# T = 0.25 GeV: lower bound of Lambda inverse slope at LHC energies
#               (ALICE, arXiv:1910.07678, Pb-Pb 0-10% at 2.76 TeV)
# T = 0.30 GeV: baseline / mid-range estimate for O-O and Pb-Pb
# T = 0.35 GeV: upper estimate for Pb-Pb at LHC energies
# T = 0.45 GeV: represents harder production channels (e.g., high-pT jets)
#               or a Tsallis tail approximation at moderate pT
# -----------------------------------------------------------------------------
if family_in_scope 5; then
    run_family_header 5 "THERMAL TEMPERATURE SCAN"

    register_job  "temp_025"  "temp_025"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.25  0.150 0.150 0.050 0.100  42
    register_job  "temp_030"  "temp_030"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    register_job  "temp_035"  "temp_035"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.35  0.150 0.150 0.050 0.100  42
    register_job  "temp_045"  "temp_045"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.45  0.150 0.150 0.050 0.100  42
fi

# -----------------------------------------------------------------------------
# Family 6: STATISTICS SCALING
# Different seeds so all N-values are statistically independent samples.
# Verifies 1/sqrt(N) error scaling and statistical stability.
# Uses the baseline cuts throughout.  Runs are seeded with different values
# so that repeated runs with the same N are statistically independent.
#
# N = 100k: ~10 seconds per run; useful for rapid iteration during debugging
# N = 500k: ~1 minute; enough to see the asymmetries clearly
# N =   2M: ~5 minutes; statistical precision comparable to a short GRID job
# N =   5M: ~10 minutes; sufficient to resolve the phi* modulation at <1% level
#
# For reference, the ring analysis signal is O(10^-3) in the observable.
# To resolve a fake signal of that size at 5 sigma with statistical error
# alone, one needs N > (5 * sigma_R / signal)^2 * N_1, where sigma_R is the
# RMS of the ring proxy distribution (~1 for unpolarized Lambdas).
# Estimate: N > (5 / 10^-3)^2 = 25 * 10^6 for 5-sigma, but since we are
# studying the FAKE signal (which can be much larger than 10^-3), the
# statistics needed here are much more modest.
# -----------------------------------------------------------------------------
if family_in_scope 6; then
    run_family_header 6 "STATISTICS SCALING"

    register_job  "stats_100k"  "stats_100k"     100000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  101
    register_job  "stats_500k"  "stats_500k"     500000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  201
    register_job  "stats_2M"    "stats_2M"       2000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  301
    register_job  "stats_5M"    "stats_5M"       5000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  401
    register_job  "stats_50M"   "stats_50M"     50000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  501
fi


# =============================================================================
# QUEUE SUMMARY
# =============================================================================

N_JOBS=${#JOB_QUEUE[@]}

BATCH_START=$(date '+%Y-%m-%d %H:%M:%S')

echo ""
echo "============================================================"
echo "  runHelicityToyModel.sh"
echo "  Batch started  : ${BATCH_START}"
echo "  Script dir     : ${SCRIPT_DIR}"
echo "  Output dir     : ${BASE_DIR}"
echo "  Log dir        : ${LOG_DIR}"
echo "  Jobs registered: ${N_JOBS}"
echo "  Max parallel   : ${MAX_PARALLEL}  (0 = unlimited / all cores)"
echo "  GNU parallel   : $([ ${USE_GNU_PARALLEL} -eq 1 ] && echo 'available' || echo 'NOT found -- using bash pool')"
echo "  Dry run        : ${DRY_RUN}"
echo "============================================================"

# Print the full job table (useful for --list and --dry-run)
echo ""
echo "  Registered jobs:"
printf "  %-4s  %-28s  %-10s  %-5s  %-14s  %-14s  %-5s\n" \
    "#" "Name" "N" "Bz" "pTmin[p,pi]" "DCA[p,pi]" "T"
printf "  %-4s  %-28s  %-10s  %-5s  %-14s  %-14s  %-5s\n" \
    "---" "----------------------------" "----------" "-----" \
    "--------------" "--------------" "-----"

IDX=0
for JOB in "${JOB_QUEUE[@]}"; do
    IFS=':' read -r NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX T \
                      PTMIN_P PTMIN_PI DCAMIN_P DCAMIN_PI SEED <<< "${JOB}"
    printf "  %-4d  %-28s  %-10s  %+.2f  [%.3f, %.3f]   [%.3f, %.3f]   %.2f\n" \
        $(( IDX + 1 )) "${NAME}" "${N}" "${BZ}" \
        "${PTMIN_P}" "${PTMIN_PI}" "${DCAMIN_P}" "${DCAMIN_PI}" "${T}"
    IDX=$(( IDX + 1 ))
done
echo ""

# --list stops here
if [[ ${LIST_ONLY} -eq 1 ]]; then
    echo "  (--list mode: exiting without running)"
    exit 0
fi

# --dry-run prints the execute_job call for each job then exits
if [[ ${DRY_RUN} -eq 1 ]]; then
    echo "  [DRY-RUN] Commands that would be executed:"
    for JOB in "${JOB_QUEUE[@]}"; do
        echo "    execute_job '${JOB}'"
    done
    echo ""
    echo "  (--dry-run mode: exiting without running)"
    exit 0
fi

if [[ ${N_JOBS} -eq 0 ]]; then
    echo "  No jobs registered.  Nothing to do."
    exit 0
fi


# =============================================================================
# DISPATCH
# =============================================================================

echo "  Dispatching ${N_JOBS} jobs (max ${MAX_PARALLEL} concurrent)..."
echo ""

if [[ ${USE_GNU_PARALLEL} -eq 1 ]]; then
    # ------------------------------------------------------------------
    # GNU parallel dispatch
    # ------------------------------------------------------------------
    # --jobs N       : limit concurrency
    # --line-buffer  : print output line-by-line (avoids interleaving)
    # --joblog FILE  : write per-job timing and status to a TSV file
    # --halt soon,fail=1 : stop launching new jobs if any job fails
    #                      but let already-running jobs finish
    # --env execute_job,BASE_DIR,LOG_DIR,...: re-export into parallel env
    # --quote prevents parallel from mangling colons in the packed string;
    #   instead we pass each job as a single argument via null-delimiter.
    #
    # Note: we use printf '%s\0' to produce NUL-delimited input so that
    # colons and other special characters in packed strings are safe.

    JOBLOG="${LOG_DIR}/parallel_joblog.tsv"

    printf '%s\0' "${JOB_QUEUE[@]}" \
    | parallel \
        --null \
        --jobs "${MAX_PARALLEL}" \
        --line-buffer \
        --joblog "${JOBLOG}" \
        --halt soon,fail=1 \
        --env execute_job \
        --env BASE_DIR \
        --env LOG_DIR \
        --env ROOT_EXE \
        --env GEN_MACRO \
        --env PLT_MACRO \
        execute_job {}

    DISPATCH_STATUS=$?

    echo ""
    if [[ ${DISPATCH_STATUS} -eq 0 ]]; then
        echo "  All jobs completed successfully."
        echo "  GNU parallel job log: ${JOBLOG}"
    else
        echo "  WARNING: One or more jobs failed (parallel exit ${DISPATCH_STATUS})."
        echo "  Check individual logs in: ${LOG_DIR}/"
        echo "  GNU parallel job log   : ${JOBLOG}"
    fi

else
    # ------------------------------------------------------------------
    # Pure bash background pool (fallback when GNU parallel is absent)
    # ------------------------------------------------------------------
    # Strategy: launch all jobs as background subshells immediately.
    # On a 192-core machine with 24 jobs this is unconditionally fine.
    # If MAX_PARALLEL is set below N_JOBS (e.g. for testing), a simple
    # slot-counter loop throttles the launch rate.
    #
    # We collect all PIDs so we can wait for them individually and
    # detect failures without aborting surviving jobs.

    declare -a PIDS=()
    declare -a NAMES=()

    IDX=0
    for JOB in "${JOB_QUEUE[@]}"; do
        # Throttle: wait until a slot is free
        while true; do
            # Count running background jobs
            RUNNING=0
            for PID in "${PIDS[@]+"${PIDS[@]}"}"; do
                if kill -0 "${PID}" 2>/dev/null; then
                    RUNNING=$(( RUNNING + 1 ))
                fi
            done
            if [[ "${MAX_PARALLEL}" -eq 0 || "${RUNNING}" -lt "${MAX_PARALLEL}" ]]; then
                break
            fi
            sleep 1
        done

        # Launch job in background subshell
        execute_job "${JOB}" &
        PIDS+=($!)

        # Extract name for the status summary
        IFS=':' read -r JOB_NAME REST <<< "${JOB}"
        NAMES+=("${JOB_NAME}")

        IDX=$(( IDX + 1 ))
    done

    echo "  All ${N_JOBS} jobs launched.  Waiting for completion..."
    echo ""

    # Wait for each job individually and collect exit statuses
    DISPATCH_STATUS=0
    for I in "${!PIDS[@]}"; do
        PID="${PIDS[${I}]}"
        NAME="${NAMES[${I}]}"
        if wait "${PID}"; then
            printf "  [OK]   %-28s  (PID %d)\n" "${NAME}" "${PID}"
        else
            STAT=$?
            printf "  [FAIL] %-28s  (PID %d, exit %d)\n" "${NAME}" "${PID}" "${STAT}"
            DISPATCH_STATUS=1
        fi
    done

    echo ""
    if [[ ${DISPATCH_STATUS} -eq 0 ]]; then
        echo "  All jobs completed successfully."
    else
        echo "  WARNING: One or more jobs failed."
        echo "  Check individual logs in: ${LOG_DIR}/"
    fi
fi


# =============================================================================
# FINAL SUMMARY
# =============================================================================

BATCH_END=$(date '+%Y-%m-%d %H:%M:%S')

echo ""
echo "============================================================"
echo "  Batch started : ${BATCH_START}"
echo "  Batch finished: ${BATCH_END}"
echo ""
echo "  Output tree:"
echo "    ROOT files : ${BASE_DIR}/<name>/helicity_<name>.root"
echo "    Plot PDFs  : ${BASE_DIR}/<name>/plots/*.pdf"
echo "    Run logs   : ${LOG_DIR}/<name>.log"
echo ""
echo "  Quick diagnostics:"
echo "    # Check all jobs succeeded:"
echo "    grep -l 'FAILED\|ERROR' ${LOG_DIR}/*.log 2>/dev/null || echo 'All clean'"
echo ""
echo "    # Print integrated R_proxy summary from every log:"
echo "    grep 'R_proxy summary' -A6 ${LOG_DIR}/*.log"
echo ""
echo "    # Check timings:"
echo "    grep 'total:' ${LOG_DIR}/*.log | sort -t: -k3 -n"
echo ""
echo "    # Open a result in ROOT:"
echo "    root -l ${BASE_DIR}/baseline/helicity_baseline.root"
echo "============================================================"
echo ""

exit ${DISPATCH_STATUS}