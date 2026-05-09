#!/usr/bin/env bash
# =============================================================================
# runHelicityToyModel.sh
# =============================================================================
#
# PURPOSE
# -------
# Coordinator script for the Lambda helicity efficiency toy model study.
# Runs helicityEfficiencyToyModel.cxx and plotHelicityEfficiency.cxx over
# ten families of parameter variations designed to isolate and quantify the
# two largest detector-induced fake-polarization effects we could think of:
#
#   (1) The Helicity Efficiency Effect (HEE): forward-backward asymmetry in
#       cos(theta*) driven by the minimum-pT threshold on daughter tracks.
#       The soft pion from a backward decay (in the Lambda rest frame) is
#       preferentially removed, biasing the distribution toward negative
#       cos(theta*). This effect should be field-independent and grows with the
#       pT threshold relative to the daughter pT spectrum.
#
#   (2) The Azimuthal Efficiency Effect (AEE): left-right asymmetry in phi*
#       driven by the interplay between the magnetic field and the minimum-DCA
#       cut applied to daughter tracks. The proton and pion, being opposite-
#       charge, curve in opposite directions; the DCA cut selects decay
#       geometries where one daughter curves away from the PV. The resulting
#       phi* modulation flips sign with the field polarity.
#       Notice that we don't have a DCA between daughters in this Toy Model,
#       so if we do manage to see the AEE effect, then it is indeed caused
#       mostly be the DCA to PV, not the DCA between daughters!
#
# To be more precise: we aim to understand what exactly causes these two effects
# in a simple Toy Model. Are their mechanisms truly the ones we believe to be?
#
# DESIGN PRINCIPLE
# ----------------
# Each family isolates ONE effect by silencing the other:
#   - AEE families (1, 2, 3):   DCA cuts active, pT cuts = 0.
#   - HEE families (5, 6):      pT cuts active, DCA cuts = 0.
#   - Context families (4, 7):  standard-ish ALICE pT+DCA cuts to see combined
#                               signal across varying kinematic conditions.
#   - Families 8, 9:            vary spectrum/window with standard-ish ALICE cuts
#                               so the fake signal is visible and comparable.
#   - Family 10 (Realistic):    full ALICE-like cuts; the combined estimate needed
#                               for the actual data comparison.
#
# OUTPUT STRUCTURE
# ----------------
# Each run produces a ROOT file with two parallel histogram families:
#   WithoutEtaGate/  -- legacy, no daughter eta requirement (kept for reference)
#   WithEtaGate/     -- physically consistent set (BOTH daughters in acceptance)
# Each family contains NoCuts / pTCutOnly / DCACutOnly / BothCuts scenarios,
# each split into EtaPos / EtaNeg / All sub-directories.
#
# PARALLELISM DESIGN
# ------------------
# All jobs from ALL families are collected into a single queue first, then
# dispatched simultaneously up to MAX_PARALLEL concurrent processes.
# On a 192-core machine with ~59 total runs, the default MAX_PARALLEL=60
# launches every job at once and the total wall time equals the slowest
# single run.
#
# Dispatch strategy (in order of preference):
#   1. GNU parallel (parallel --jobs N): preferred; handles job tracking,
#      per-job logging via --joblog, and clean failure reporting.
#   2. Pure bash background pool: fallback if parallel is not installed;
#      uses the same subshell structure with a lightweight semaphore.
#
# PARAMETER FAMILIES
# ------------------
#   Family 0: BASELINE
#       One reference run with NO cuts at all (pT_min = DCA_min = 0).
#       All other families should be compared against this flat reference.
#       In this family, we should not see any difference between the "pTcut"
#       or "DCAcut" TDirectories: the "cut" values are zeroed when calling the
#       Toy Model, after all!
#
#   Family 1: ASYMMETRIC DCA CUTS  [AEE probe]
#       Vary DCA cuts asymmetrically between proton and pion.
#       pT cuts = 0. Tests whether the AEE sign depends on which daughter's
#       DCA dominates, and whether a proton-only vs pion-only cut gives
#       opposite phi* asymmetry directions.
#
#   Family 2: SYMMETRIC DCA CUTS  [AEE strength scan]
#       Equal DCA cuts on both daughters, ranging from 0.05 to 0.30 cm.
#       pT cuts = 0. Measures AEE magnitude as a clean function of cut strength.
#
#   Family 3: MAGNETIC FIELD  [AEE field dependence]
#       Fixed standard DCA cuts (0.05 / 0.05 cm). pT cuts = 0.
#       Sign flip test (+0.5 vs -0.5 T) plus strength scan from near-zero
#       (Bz = 0.0001 T, where AEE must vanish) to 1.0 T.
#       Note: Bz = 0.0001 T uses the straight-track geometric DCA (no bending),
#       which is charge-independent, so the phi* asymmetry correctly goes to zero.
#       That small field is a trick to not use Bz = 0, which would cause a division
#       by zero in the code!
#
#   Family 4: LAMBDA GENERATOR PT MINIMUM  [kinematic context for HEE+AEE]
#       Standard-ish ALICE cuts (pT = 0.15 GeV/c, DCA = 0.05/0.05 cm).
#       Varies the minimum Lambda pT used in generation to probe how the
#       combined fake signal depends on the Lambda pT regime.
#
#   Family 5: DAUGHTER PT CUTS  [HEE primary probe]
#       DCA cuts = 0. Varies pT thresholds symmetrically and asymmetrically.
#       Pion-only cuts give the cleanest HEE isolation (the pion is always the
#       softer daughter). Proton-only cuts serve as a control: because the
#       proton is harder, the threshold is less likely to bite, so the HEE
#       should be much weaker. Have to thoroughly test and check this anyways.
#
#   Family 6: ETA ACCEPTANCE WINDOW  [both effects vs acceptance]
#       Standard-ish ALICE cuts (pT = 0.15 GeV/c, DCA = 0.05/0.05 cm).
#       Varies etaMaxDetector from 0.5 to 3.0. Shows how the fake signal
#       magnitude changes with the detector acceptance window.
#
#   Family 7: TEMPERATURE SCAN  [HEE+AEE vs Lambda pT spectrum shape]
#       Standard-ish ALICE cuts (pT = 0.15 GeV/c, DCA = 0.05/0.05 cm).
#       Varies the Boltzmann mT temperature from 0.20 to 0.50 GeV.
#       Lower T -> softer Lambdas -> softer daughters -> stronger HEE.
#
#   Family 8: RING KINEMATIC WINDOWS  [fake signal per pT window]
#       Standard-ish ALICE cuts (pT = 0.15 GeV/c, DCA = 0.05/0.05 cm).
#       Splits the Lambda pT spectrum into successive windows to identify
#       where the fake signal is largest (relevant for the ring observable).
#       The kinematic cuts are those of Vitor's PRC paper (PRC 109, 014905 (2024)).
#
#   Family 9: REALISTIC ALICE CUTS  [combined HEE+AEE estimate]
#       Mimics what I did in experimental V0 selection conditions, to some degree.
#       Combines pT and DCA cuts at loose / standard / tight working points, with
#       both field polarities and in the ring analysis kinematic window.
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
#     0_Baseline/          baseline/
#     1_AsymDCA/           dca_p_only_005/ dca_pi_only_005/ dca_pi_only_010/
#                          dca_asym_std/ dca_asym_wide/ dca_asym_rev/
#     2_SymDCA/            dca_sym_005/ ... dca_sym_030/
#     3_BField/            field_vlow/ field_b010/ ... field_b100/ field_b050neg/
#     4_LamPtMin/          lam_ptmin_000/ ... lam_ptmin_200/
#     5_DaughterPt/        pt_sym_005/ ... pt_sym_300/
#                          pt_pi_010/ pt_pi_150/ pt_pi_200/ pt_p_150/
#     6_EtaMax/            eta_050/ eta_070/ eta_090/ eta_120/ eta_150/
#                          eta_200/ eta_300/
#     7_Temperature/       temp_020/ ... temp_050/
#     8_KinWindow/         win_inclusive/ win_vsoft/ win_soft/ win_ring/
#                          win_mid/ win_hard/ win_vhard/
#     9_RealisticAlice/    alice_loose/ alice_std/ alice_std_neg/
#                          alice_tight/ alice_ring/
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
#   Each ROOT file (histograms only, no TTrees): ~4-8 MB after ROOT zlib
#   compression (two histogram families: WithEtaGate + WithoutEtaGate).
#   Total for ~59 runs: ~600 MB. (EYEBALLED! Didn't measure these.)
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
# Default: 64 (all jobs at once on a 192-core machine with ~63 total runs).
# Set to 0 to let GNU parallel use all available cores (parallel -j 0),
# or to any positive integer to cap concurrency explicitly.
MAX_PARALLEL=64

# Default Lambda count per run.
# At ~90 s per 10M Lambdas, 1B = ~900 s (~15 min) per run.
# With ~59 runs in parallel on 192 cores, total wall time ~ 15 min.
DEFAULT_N=10000000000
# DEFAULT_N=10000000 # For testing only (~30 s per run)

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
            echo "Unknown argument: $1. Use --help for usage." >&2
            exit 1
            ;;
    esac
done

# Default: run all families
if [[ ${#FAMILIES_TO_RUN[@]} -eq 0 ]]; then
    FAMILIES_TO_RUN=(0 1 2 3 4 5 6 7 8 9 10)
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
# Format: colon-delimited, 14 fields matching helicityEfficiencyToyModel args:
#   NAME:SUBDIR:N:BZ:PTMIN_LAM:PTMAX_LAM:RAPMAX:ETAMAX:T:PTMIN_P:PTMIN_PI:DCAMIN_P:DCAMIN_PI:SEED
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

# register_job: same 14-argument signature as the old run_one,
# but only adds to the queue -- does NOT launch anything.
#
#   $1  NAME       short identifier, used for output subdir and log name
#   $2  SUBDIR     subdirectory name under BASE_DIR  (usually same as NAME)
#   $3  N          number of Lambdas
#   $4  BZ         magnetic field [T]
#   $5  PTMIN_LAM  Lambda min pT [GeV/c]
#   $6  PTMAX_LAM  Lambda max pT [GeV/c]
#   $7  RAPMAX     Lambda max |rapidity|
#   $8  ETAMAX     Daughter maximum detectable |eta|
#   $9  T          Boltzmann temperature [GeV]
#   $10 PTMIN_P    proton min pT [GeV/c]
#   $11 PTMIN_PI   pion min pT [GeV/c]
#   $12 DCAMIN_P   proton min DCA_xy [cm]
#   $13 DCAMIN_PI  pion min DCA_xy [cm]
#   $14 SEED       TRandom3 seed
register_job() {
    local PACKED="$1:$2:$3:$4:$5:$6:$7:$8:$9:${10}:${11}:${12}:${13}:${14}"
    JOB_QUEUE+=("${PACKED}")
    JOB_FAMILIES+=("${CURRENT_FAMILY}")
}

# run_family_header: cosmetic only -- no wait, no throttle.
# Records the current family number so register_job can tag each entry.
CURRENT_FAMILY=0
run_family_header() {
    CURRENT_FAMILY="$1"
    # echo ""
    # echo "------------------------------------------------------------"
    # echo "  Registering Family ${1}: ${2}"
    # echo "------------------------------------------------------------"
}


# =============================================================================
# JOB EXECUTOR (called once per job, in a subshell by the dispatcher)
# =============================================================================
# This function is exported so GNU parallel can call it in a child process.
# It receives one colon-delimited packed argument string.

execute_job() {
    # Unpack the 14 fields
    IFS=':' read -r \
        NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX ETAMAX T \
        PTMIN_P PTMIN_PI DCAMIN_P DCAMIN_PI SEED \
        <<< "$1"

    local OUTDIR="${BASE_DIR}/${SUBDIR}"
    local ROOT_FILE="${OUTDIR}/helicity_${NAME}.root"
    local LOG_FILE="${LOG_DIR}/${NAME}.log"

    mkdir -p "${OUTDIR}"

    # Build ROOT call strings.
    # Single quotes around the macro path protect spaces; the inner double
    # quotes around the ROOT file path protect spaces in BASE_DIR.
    local GEN_CALL="${N},\"${ROOT_FILE}\",${BZ},${PTMIN_LAM},${PTMAX_LAM},${RAPMAX},${ETAMAX},${T},${PTMIN_P},${PTMIN_PI},${DCAMIN_P},${DCAMIN_PI},${SEED}"
    local PLT_CALL="\"${ROOT_FILE}\""

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
#
# Column guide (all register_job calls):
#   NAME  SUBDIR  N  Bz  pTminLam  pTmaxLam  rapMax  etaMax  T  pTp  pTpi  dcaP  dcaPi  seed
#
# Shared baseline parameter values (always used unless the family varies them):
#   Bz         = +0.5 T
#   pTminLam   = 0.0 GeV/c
#   pTmaxLam   = 10.0 GeV/c
#   rapMax     = 5.0          (wide enough to include all daughters inside etaMax=3.0)
#   etaMax     = 0.9          (ALICE inner barrel)
#   T          = 0.30 GeV
#   pTp/pTpi   = 0.0 GeV/c    (no pT cut)
#   dcaP/dcaPi = 0.0 cm       (no DCA cut)
#   seed       = 0
#
# Standard-ish ALICE working point (used by context families 4, 6, 7, 8, 9):
# (the actual minimum for Lambdas would actually be around 0.3. Less than that has really bad efficiency!)
#   pTp = pTpi = 0.15 GeV/c
#   dcaP = 0.05 cm,  dcaPi = 0.05 cm
# =============================================================================

# -----------------------------------------------------------------------------
# FAMILY 0: BASELINE
# All cuts disabled. Every other family is compared against this run.
# The WithEtaGate/NoCuts histograms must be flat in both cos(theta*) and phi*.
# Any asymmetry in any other family beyond what appears here is artefactual.
# -----------------------------------------------------------------------------
if family_in_scope 0; then
    run_family_header 0 "BASELINE (no cuts, all defaults)"

    #              NAME        SUBDIR                    N             Bz    pTminL pTmaxL  rap   eta    T    pTp   pTpi  dcaP  dcaPi seed
    register_job  "baseline"  "0_Baseline/baseline"  ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.0   0.0   0
fi

# -----------------------------------------------------------------------------
# FAMILY 1: ASYMMETRIC DCA CUTS  [AEE isolation -- daughter-level asymmetry]
# DCA cuts applied to ONE or both daughters with unequal thresholds.
# pT cuts = 0 throughout to isolate the AEE from the HEE.
#
# Rationale:
#   The AEE arises because the proton (+) and pion (-) bend in opposite
#   directions and the DCA cut selects different phi* configurations for each.
#   Applying the DCA cut to ONLY the proton or ONLY the pion reveals which
#   daughter drives the phi* asymmetry and whether the effect is additive.
#   The "reversed" run (proton cut tighter than pion) tests whether swapping
#   the cut ratio flips the phi* asymmetry direction.
#
# Expected results:
#   dca_p_only_005:  non-zero phi* asymmetry from proton DCA alone
#   dca_pi_only_005: non-zero phi* asymmetry from pion DCA alone; sign opposite
#                    to proton-only because pion has opposite charge and thus
#                    curves the opposite way in the same field
#   dca_pi_only_010: same sign as dca_pi_only_005, larger magnitude
#   dca_asym_std:    ALICE standard (0.05 / 0.05 cm) -- reference for Family 3
#   dca_asym_wide:   larger pion cut; stronger AEE
#   dca_asym_rev:    proton cut > pion cut; should partially cancel or invert
#                    the standard asymmetry
# -----------------------------------------------------------------------------
if family_in_scope 1; then
    run_family_header 1 "ASYMMETRIC DCA CUTS  [AEE: daughter-level asymmetry]"

    #               NAME                   SUBDIR                                  N             Bz    pTminL pTmaxL  rap   eta    T    pTp   pTpi  dcaP   dcaPi  seed
    register_job   "dca_p_only_005"       "1_AsymDCA/dca_p_only_005"           ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.000  0
    register_job   "dca_pi_only_005"      "1_AsymDCA/dca_pi_only_005"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.000  0.050  0
    register_job   "dca_pi_only_010"      "1_AsymDCA/dca_pi_only_010"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.000  0.100  0
    register_job   "dca_asym_std"         "1_AsymDCA/dca_asym_std"             ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.100  0
    register_job   "dca_asym_wide"        "1_AsymDCA/dca_asym_wide"            ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.200  0
    register_job   "dca_asym_rev"         "1_AsymDCA/dca_asym_rev"             ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.100  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 2: SYMMETRIC DCA CUTS  [AEE strength scan]
# Equal DCA cuts on both daughters. pT cuts = 0.
#
# Rationale:
#   Removing the asymmetry between daughters lets us study the AEE magnitude
#   as a clean function of cut strength without confounding from different
#   effective geometries. Starting at 0.05 cm (the minimum ALICE cut) and
#   progressing to 0.30 cm.
#
# Expected result:
#   A monotonically increasing phi* asymmetry magnitude with tighter cuts,
#   since tighter DCA cuts select more displaced decays where magnetic bending
#   has a larger differential effect on the two opposite-charge daughters.
# -----------------------------------------------------------------------------
if family_in_scope 2; then
    run_family_header 2 "SYMMETRIC DCA CUTS  [AEE: cut strength scan]"

    #               NAME              SUBDIR                          N             Bz    pTminL pTmaxL  rap   eta    T    pTp   pTpi  dcaP   dcaPi  seed
    register_job   "dca_sym_005"     "2_SymDCA/dca_sym_005"       ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "dca_sym_010"     "2_SymDCA/dca_sym_010"       ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.100  0.100  0
    register_job   "dca_sym_015"     "2_SymDCA/dca_sym_015"       ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.150  0.150  0
    register_job   "dca_sym_020"     "2_SymDCA/dca_sym_020"       ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.200  0.200  0
    register_job   "dca_sym_030"     "2_SymDCA/dca_sym_030"       ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.300  0.300  0
fi

# -----------------------------------------------------------------------------
# FAMILY 3: MAGNETIC FIELD  [AEE field dependence]
# Fixed standard DCA cuts (0.05 / 0.05 cm). pT cuts = 0.
#
# Rationale:
#   (a) Sign flip: comparing field_b050 and field_b050neg tests the fundamental
#       prediction that the phi* asymmetry sign flips with the field polarity.
#       The cos(theta*) HEE (which is field-independent) must NOT flip.
#       ALICE has operated at both +0.5 T and -0.5 T in Run 2, but not for OO
#       in Run 3 (AFAIK), so this is just a test of the AEE phenomenon.
#
#   (b) Strength scan: the DCA-based AEE depends on the helix curvature
#       radius R = pT / (q * kBConv * Bz). A stronger field means tighter
#       helices and larger differential DCA between the two daughters for the
#       same decay geometry. Larger |Bz| -> stronger AEE.
#
#   (c) Near-zero field (Bz = 0.0001 T): the helix radius is ~km-scale,
#       so both daughters travel in straight lines. The DCA reduces to the
#       charge-independent geometric impact parameter (perpendicular distance
#       from origin to the straight-line track). Since this is identical for
#       both daughters (same vertex, symmetric to charge), there is no phi*
#       asymmetry. The AEE must go to zero.
#       Numerically: DCA_xy -> |xv*(py/pT) - yv*(px/pT)|, which is well-
#       defined and non-zero in general (decay geometry is not symmetric with
#       respect to the origin), just charge-independent.
# -----------------------------------------------------------------------------
if family_in_scope 3; then
    run_family_header 3 "MAGNETIC FIELD  [AEE: field sign and strength]"

    #               NAME               SUBDIR                          N             Bz       pTminL pTmaxL  rap   eta    T    pTp   pTpi  dcaP   dcaPi  seed
    register_job   "field_vlow"       "3_BField/field_vlow"        ${DEFAULT_N}   0.0001   0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b010"       "3_BField/field_b010"        ${DEFAULT_N}   0.10     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b020"       "3_BField/field_b020"        ${DEFAULT_N}   0.20     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b030"       "3_BField/field_b030"        ${DEFAULT_N}   0.30     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b050"       "3_BField/field_b050"        ${DEFAULT_N}   0.50     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b075"       "3_BField/field_b075"        ${DEFAULT_N}   0.75     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b100"       "3_BField/field_b100"        ${DEFAULT_N}   1.00     0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b050neg"    "3_BField/field_b050neg"     ${DEFAULT_N}   -0.50    0.0   10.0   5.0   0.9   0.30  0.0   0.0   0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 4: LAMBDA GENERATOR PT MINIMUM  [kinematic context for HEE+AEE]
# Standard-ish ALICE cuts (pTp=pTpi=0.15, dcaP=0.05, dcaPi=0.05).
#
# Rationale:
#   Both the HEE and AEE depend on how soft the Lambda daughters are in the
#   lab frame. The daughter pT scales roughly with the Lambda pT (through the
#   Lorentz boost). Raising the minimum generated Lambda pT shifts the
#   population to harder Lambdas with harder daughters, where:
#     - The pT threshold is less likely to bite: HEE weaker.
#     - The decay vertices are more displaced (longer c*tau * beta*gamma): AEE
#       might be stronger because the daughters have more time to separate in
#       DCA space before hitting the acceptance boundary.
###         This is actually one of the nicest things we aim to check in this Toy
###         Model!!! Is AEE an effect that comes from reconstruction algorithms
###         prioritizing smaller DCA between the daughters (which is not something
###         covered at all by this Toy Model: DCAdau = 0 by definition. The daughters
###         come exactly from the same point and we have the equivalent to infinite
###         momentum and spatial resolution in the Toy Model), or is AEE an effect
###         that comes from the DCA to PV of the two daughters? If we see AEE here,
###         then it comes (or at least part of it does) from the DCA to the PV, not
###         from the DCA between the daughters of the Lambda!
#       
#   With standard-ish ALICE cuts active, this family shows the combined HEE+AEE
#   as a function of the Lambda pT regime.
# -----------------------------------------------------------------------------
if family_in_scope 4; then
    run_family_header 4 "LAMBDA GENERATOR PT MINIMUM  [kinematic context, standard ALICE cuts]"

    #               NAME                SUBDIR                            N             Bz    pTminL  pTmaxL  rap   eta    T    pTp    pTpi   dcaP   dcaPi  seed
    register_job   "lam_ptmin_000"     "4_LamPtMin/lam_ptmin_000"     ${DEFAULT_N}   0.5   0.000   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_030"     "4_LamPtMin/lam_ptmin_030"     ${DEFAULT_N}   0.5   0.300   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_050"     "4_LamPtMin/lam_ptmin_050"     ${DEFAULT_N}   0.5   0.500   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_100"     "4_LamPtMin/lam_ptmin_100"     ${DEFAULT_N}   0.5   1.000   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_150"     "4_LamPtMin/lam_ptmin_150"     ${DEFAULT_N}   0.5   1.500   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_200"     "4_LamPtMin/lam_ptmin_200"     ${DEFAULT_N}   0.5   2.000   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 5: DAUGHTER PT CUTS  [HEE primary probe]
# DCA cuts = 0 throughout to isolate the HEE from the AEE.
#
# Rationale:
#   The HEE grows when the pT threshold removes a significant fraction of
#   backward-decay pions (soft in lab) while leaving forward-decay protons
#   intact. The key observables are:
#     - Symmetric cuts (both daughters): combined threshold effect.
#     - Pion-only cuts: cleanest HEE isolation. The pion is ALWAYS softer
#       than the proton in the Lambda rest frame (mass hierarchy), so cutting
#       on pion pT alone is the purest driver of the cos(theta*) asymmetry.
#     - Proton-only cuts: control case. The proton carries most of the
#       Lambda momentum, so it is harder; the threshold bites less and the
#       HEE should be significantly weaker than for pion-only cuts at the
#       same threshold value.
#
# Expected results:
#   Symmetric: monotonically increasing |<cos(theta*)>| effects with tighter cuts.
#   Pion-only: similar sign and magnitude to symmetric, but slightly larger
#              because only the pion asymmetry contributes without dilution.
#   Proton-only: same sign (or zero) but much weaker HEE.
# -----------------------------------------------------------------------------
if family_in_scope 5; then
    run_family_header 5 "DAUGHTER PT CUTS  [HEE: primary probe]"

    # -- Symmetric cuts --
    #               NAME              SUBDIR                          N             Bz    pTminL pTmaxL  rap   eta    T    pTp    pTpi   dcaP  dcaPi  seed
    register_job   "pt_sym_005"      "5_DaughterPt/pt_sym_005"    ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.050  0.050  0.0   0.0   0
    register_job   "pt_sym_010"      "5_DaughterPt/pt_sym_010"    ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.100  0.100  0.0   0.0   0
    register_job   "pt_sym_150"      "5_DaughterPt/pt_sym_150"    ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.150  0.150  0.0   0.0   0
    register_job   "pt_sym_200"      "5_DaughterPt/pt_sym_200"    ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.200  0.200  0.0   0.0   0
    register_job   "pt_sym_300"      "5_DaughterPt/pt_sym_300"    ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.300  0.300  0.0   0.0   0
    # -- Pion-only cuts (purest HEE isolation) --
    register_job   "pt_pi_010"       "5_DaughterPt/pt_pi_010"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.000  0.100  0.0   0.0   0
    register_job   "pt_pi_150"       "5_DaughterPt/pt_pi_150"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.000  0.150  0.0   0.0   0
    register_job   "pt_pi_200"       "5_DaughterPt/pt_pi_200"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.000  0.200  0.0   0.0   0
    # -- Proton-only cut (control: expect weaker HEE than pion-only) --
    register_job   "pt_p_150"        "5_DaughterPt/pt_p_150"      ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.150  0.000  0.0   0.0   0
fi

# -----------------------------------------------------------------------------
# FAMILY 6: ETA ACCEPTANCE WINDOW  [HEE+AEE vs detector acceptance]
# Standard-ish ALICE cuts (pTp=pTpi=0.15, dcaP=0.05, dcaPi=0.05).
#
# Rationale:
#   etaMaxDetector controls both the eta gate applied inside WithEtaGate and
#   the kinematic range of Lambda daughters that are accepted. A narrower
#   window (small etaMax) accepts mostly central (|eta|~0) daughters, which
#   tend to come from higher-pT Lambdas with harder daughters -- both effects
#   should be weaker. A wider window (large etaMax) includes softer forward
#   daughters and more extreme decay geometries -- both effects should grow.
#
#   The comparison between eta_090 (ALICE inner barrel) and wider values is
#   important for understanding the systematics of forward-tracking.
#   eta_090 duplicates the baseline kinematics but with cuts active; it serves
#   as the reference point within this family.
# -----------------------------------------------------------------------------
if family_in_scope 6; then
    run_family_header 6 "ETA ACCEPTANCE WINDOW  [HEE+AEE vs detector acceptance, standard ALICE cuts]"

    #               NAME           SUBDIR                        N             Bz    pTminL pTmaxL  rap   etaMax  T    pTp    pTpi   dcaP   dcaPi  seed
    register_job   "eta_050"      "6_EtaMax/eta_050"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.50    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_070"      "6_EtaMax/eta_070"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.70    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_090"      "6_EtaMax/eta_090"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_120"      "6_EtaMax/eta_120"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   1.20    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_150"      "6_EtaMax/eta_150"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   1.50    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_200"      "6_EtaMax/eta_200"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   2.00    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_300"      "6_EtaMax/eta_300"          ${DEFAULT_N}   0.5   0.0   10.0   5.0   3.00    0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 7: TEMPERATURE SCAN  [HEE+AEE vs Lambda pT spectrum shape]
# Standard ALICE cuts (pTp=pTpi=0.15, dcaP=0.05, dcaPi=0.05).
#
# Rationale:
#   The Boltzmann mT temperature T controls the shape of the Lambda pT spectrum.
#   Lower T -> softer Lambdas -> softer daughters in the lab frame.
#   Softer daughters are more likely to fail the pT cut -> stronger HEE.
#   They also produce smaller decay radii (shorter c*tau*betagamma) but with
#   smaller DCA from weaker bending -> the net effect on AEE is non-trivial.
#   Values of T ~ 0.25-0.35 GeV bracket the measured Lambda inverse slope at
#   LHC energies. T ~ 0.50 GeV represents a harder production scenario
#   (e.g. fragmentation-dominated high-pT region).
#   This is mostly out of curiosity. I expect no real insight out of this haha
# -----------------------------------------------------------------------------
if family_in_scope 7; then
    run_family_header 7 "TEMPERATURE SCAN  [HEE+AEE vs Lambda pT spectrum, standard ALICE cuts]"

    #               NAME           SUBDIR                          N             Bz    pTminL pTmaxL  rap   eta    T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "temp_020"     "7_Temperature/temp_020"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.20  0.150  0.150  0.050  0.050  0
    register_job   "temp_025"     "7_Temperature/temp_025"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.25  0.150  0.150  0.050  0.050  0
    register_job   "temp_030"     "7_Temperature/temp_030"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "temp_035"     "7_Temperature/temp_035"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.35  0.150  0.150  0.050  0.050  0
    register_job   "temp_040"     "7_Temperature/temp_040"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.40  0.150  0.150  0.050  0.050  0
    register_job   "temp_050"     "7_Temperature/temp_050"     ${DEFAULT_N}   0.5   0.0   10.0   5.0   0.9   0.50  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 8: RING KINEMATIC WINDOWS  [fake signal magnitude per Lambda pT window]
# Standard-ish ALICE cuts (pTp=pTpi=0.15, dcaP=0.05, dcaPi=0.05).
#
# Rationale:
#   The Hydro ring observable analysis suggests a specific kinematic window. By
#   restricting the generated Lambdas to successive pT bands we can identify
#   in which band the fake signal (AEE + HEE) is largest, and whether the
#   ring-analysis window [0.5, 1.5] GeV/c is particularly contaminated.
#
#   win_inclusive: full spectrum, same kinematic parameters as baseline;
#                  serves as the reference point within this family.
#   win_vsoft:     [0.0, 0.5] GeV/c -- very soft Lambdas, daughters near pT cut
#   win_soft:      [0.0, 1.0] GeV/c -- includes the soft regime
#   win_ring:      [0.5, 1.5] GeV/c -- ring analysis window (PhysRevC.109.014905)
#   win_mid:       [1.0, 3.0] GeV/c -- moderate pT
#   win_hard:      [1.5, 4.0] GeV/c -- hard regime
#   win_vhard:     [3.0,10.0] GeV/c -- very hard; daughters well above any cut
# -----------------------------------------------------------------------------
if family_in_scope 8; then
    run_family_header 8 "RING KINEMATIC WINDOWS  [fake signal per Lambda pT band, standard ALICE cuts]"

    #               NAME                SUBDIR                              N             Bz    pTminL  pTmaxL   rap   eta    T    pTp    pTpi   dcaP   dcaPi  seed
    register_job   "win_inclusive"     "8_KinWindow/win_inclusive"      ${DEFAULT_N}   0.5   0.000   10.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_vsoft"         "8_KinWindow/win_vsoft"          ${DEFAULT_N}   0.5   0.000    0.500   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_soft"          "8_KinWindow/win_soft"           ${DEFAULT_N}   0.5   0.000    1.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_ring"          "8_KinWindow/win_ring"           ${DEFAULT_N}   0.5   0.500    1.500   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_mid"           "8_KinWindow/win_mid"            ${DEFAULT_N}   0.5   1.000    3.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_hard"          "8_KinWindow/win_hard"           ${DEFAULT_N}   0.5   1.500    4.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_vhard"         "8_KinWindow/win_vhard"          ${DEFAULT_N}   0.5   3.000   10.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 9: REALISTIC ALICE CUTS  [combined HEE+AEE estimate for data correction]
# Mimics actual ALICE V0 selection conditions across working points and field
# polarities. This is the family whose <R_proxy> numbers matter for understanding
# the fake polarization effects in the real ring observable analysis.
#
# alice_loose:    Similar to ITS-track regime; loose DCA, lower pT threshold
# alice_std:      ALICE standard-ish V0 selection
# alice_std_neg:  Same as alice_std but Bz = -0.5 T; the AEE must flip, the
#                 HEE must not. Maybe interesting for a Run 2 combined-polarity dataset.
# alice_tight:    Tighter selection for a central Pb-Pb; larger
#                 pT and DCA thresholds produce larger combined fake signal.
# alice_ring:     Standard-ish selection restricted to the ring analysis kinematic
#                 window [0.5, 1.5] GeV/c; the most directly relevant number
#                 for the ring polarization fake signal study.
# -----------------------------------------------------------------------------
if family_in_scope 9; then
    run_family_header 9 "REALISTIC ALICE CUTS  [combined HEE+AEE estimate for data correction]"

    #               NAME                  SUBDIR                                N             Bz     pTminL  pTmaxL   rap   eta    T    pTp    pTpi   dcaP   dcaPi  seed
    register_job   "alice_loose"         "9_RealisticAlice/alice_loose"      ${DEFAULT_N}   0.50   0.000   10.000   5.0   0.9   0.30  0.100  0.100  0.020  0.020  0
    register_job   "alice_std"           "9_RealisticAlice/alice_std"        ${DEFAULT_N}   0.50   0.000   10.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "alice_std_neg"       "9_RealisticAlice/alice_std_neg"    ${DEFAULT_N}   -0.50  0.000   10.000   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
    register_job   "alice_tight"         "9_RealisticAlice/alice_tight"      ${DEFAULT_N}   0.50   0.000   10.000   5.0   0.9   0.30  0.200  0.200  0.100  0.100  0
    register_job   "alice_ring"          "9_RealisticAlice/alice_ring"       ${DEFAULT_N}   0.50   0.500    1.500   5.0   0.9   0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 10: MAGNETIC FIELD  [AEE field dependence + HEE with field dependence]
# Fixed standard DCA cuts (0.05 / 0.05 cm). pT cuts = 0.2 GeV/c
# -----------------------------------------------------------------------------
if family_in_scope 10; then
    run_family_header 10 "MAGNETIC FIELD  [AEE: field sign and strength, with HEE tests using min pT cuts]"

    #               NAME               SUBDIR                          N             Bz       pTminL pTmaxL  rap   eta    T    pTp   pTpi  dcaP   dcaPi  seed
    register_job   "field_vlow"       "10_BField/field_vlow"        ${DEFAULT_N}   0.0001   0.0   10.0   5.0   0.9   0.30  0.200   0.200   0.050  0.050  0
    register_job   "field_b050"       "10_BField/field_b050"        ${DEFAULT_N}   0.50     0.0   10.0   5.0   0.9   0.30  0.200   0.200   0.050  0.050  0
    register_job   "field_b100"       "10_BField/field_b100"        ${DEFAULT_N}   1.00     0.0   10.0   5.0   0.9   0.30  0.200   0.200   0.050  0.050  0
    register_job   "field_b050neg"    "10_BField/field_b050neg"     ${DEFAULT_N}   -0.50    0.0   10.0   5.0   0.9   0.30  0.200   0.200   0.050  0.050  0
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
printf "  %-4s  %-28s  %-10s  %-6s  %-14s  %-14s  %-5s\n" \
    "#" "Name" "N" "Bz" "pTmin[p,pi]" "DCA[p,pi]" "T"
printf "  %-4s  %-28s  %-10s  %-6s  %-14s  %-14s  %-5s\n" \
    "---" "----------------------------" "----------" "------" \
    "--------------" "--------------" "-----"

IDX=0
for JOB in "${JOB_QUEUE[@]}"; do
    IFS=':' read -r NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX ETAMAX T \
                      PTMIN_P PTMIN_PI DCAMIN_P DCAMIN_PI SEED <<< "${JOB}"
    printf "  %-4d  %-28s  %-10s  %+.4f  [%.3f, %.3f]   [%.3f, %.3f]   %.2f\n" \
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
    echo "  No jobs registered. Nothing to do."
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
    # On a 192-core machine with ~59 jobs this is unconditionally fine.
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

    echo "  All ${N_JOBS} jobs launched. Waiting for completion..."
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
echo "    ROOT files : ${BASE_DIR}/<family>/<name>/helicity_<name>.root"
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
echo "    # Compare AEE sign flip between positive and negative field:"
echo "    grep -A4 'WithEtaGate.*DCACutOnly' ${LOG_DIR}/field_b050.log"
echo "    grep -A4 'WithEtaGate.*DCACutOnly' ${LOG_DIR}/field_b050neg.log"
echo ""
echo "    # Compare HEE for pion-only vs proton-only pT cuts:"
echo "    grep 'pTCutOnly' ${LOG_DIR}/pt_pi_150.log ${LOG_DIR}/pt_p_150.log"
echo ""
echo "    # Open a result in ROOT:"
echo "    root -l ${BASE_DIR}/0_Baseline/baseline/helicity_baseline.root"
echo "============================================================"
echo ""

exit ${DISPATCH_STATUS}