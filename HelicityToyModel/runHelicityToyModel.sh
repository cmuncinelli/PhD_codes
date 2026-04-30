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
#       The DCA-asymmetry sign must flip with the field.  This is the most
#       direct test of the mechanism: a positive field should give a positive
#       <R_proxy> at eta>0 and negative at eta<0 (or vice versa), and the
#       signs must swap when B is reversed.  Any other behavior indicates
#       that the asymmetry has a different origin.
#       Reference: The ALICE solenoidal magnet operated at both polarities
#       during Run 2/3, making this comparison experimentally accessible.
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
#       (see ALICE-PUBLIC-2017-005 and the strangeness PWG recommendations).
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
#   ./runHelicityToyModel.sh                   # Run all families
#   ./runHelicityToyModel.sh --family 1        # Run only Family 1
#   ./runHelicityToyModel.sh --family 2 3      # Run Families 2 and 3
#   ./runHelicityToyModel.sh --dry-run         # Print commands without running
#   ./runHelicityToyModel.sh --help            # Show this help
#
# PARALLELISM
# -----------
# Runs within the same family are launched in the background (parallel) by
# default, limited to MAX_PARALLEL simultaneous ROOT processes to avoid
# overloading the machine.  Runs from different families are sequential
# (the script waits for a family to finish before starting the next).
# Set MAX_PARALLEL=1 to force fully sequential execution.
#
# OUTPUT STRUCTURE
# ----------------
#   BASE_DIR/
#     baseline/
#       helicity_baseline.root
#       plots/
#     field_pos/     field_neg/
#     ptcut_010/     ptcut_015/     ptcut_020/     ptcut_030/
#     dcacut_none/   dcacut_loose/  dcacut_std/    dcacut_tight/
#     window_full/   window_ring/   window_hard/
#     temp_025/      temp_030/      temp_035/      temp_045/
#     stats_100k/    stats_500k/    stats_2M/      stats_5M/
#     logs/
#       <run_name>.log
#
# REQUIREMENTS
# ------------
#   - ROOT 6.14+ with C++14 support
#   - The two .cxx files must be in the same directory as this script
#     (or in SCRIPT_DIR as set below)
#   - Sufficient disk space: each ROOT file is ~5-20 MB depending on N
#
# =============================================================================

set -euo pipefail   # Exit on error, undefined variable, or pipe failure

# =============================================================================
# CONFIGURATION -- edit these paths as needed
# =============================================================================

# Directory where this script and the .cxx files live
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Base output directory for all ROOT files and plot folders
BASE_DIR="/home/users/cicerodm/RingPol/HelicityToyModel"

# Log directory
LOG_DIR="${BASE_DIR}/logs"

# Maximum number of simultaneous ROOT processes within one family
MAX_PARALLEL=8

# Default number of Lambdas per run (overridden by Family 6)
DEFAULT_N=100000000 # Bumped from 1.000.000 to 100.000.000 (10.000.000 takes about 1m30s per run, and we have only 23 runs)

# ROOT executable (set to full path if not in PATH, e.g. /opt/root/bin/root)
ROOT_EXE="root"

# Generator macro and plotting macro (relative to SCRIPT_DIR)
GEN_MACRO="${SCRIPT_DIR}/helicityEfficiencyToyModel.cxx"
PLT_MACRO="${SCRIPT_DIR}/plotHelicityEfficiency.cxx"

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

DRY_RUN=0
FAMILIES_TO_RUN=()   # Empty = run all

print_help() {
    sed -n '/^# USAGE/,/^# OUTPUT/p' "$0" | head -n -2 | sed 's/^# //'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)      print_help ;;
        --dry-run|-n)   DRY_RUN=1; shift ;;
        --family|-f)
            shift
            while [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; do
                FAMILIES_TO_RUN+=("$1")
                shift
            done
            ;;
        *) echo "Unknown argument: $1.  Use --help for usage." >&2; exit 1 ;;
    esac
done

# If no family specified, run all (0 through 6)
if [[ ${#FAMILIES_TO_RUN[@]} -eq 0 ]]; then
    FAMILIES_TO_RUN=(0 1 2 3 4 5 6)
fi

# =============================================================================
# SETUP
# =============================================================================

# Create base and log directories
mkdir -p "${BASE_DIR}"
mkdir -p "${LOG_DIR}"

# Verify the macros exist
for MACRO in "${GEN_MACRO}" "${PLT_MACRO}"; do
    if [[ ! -f "${MACRO}" ]]; then
        echo "ERROR: Macro not found: ${MACRO}" >&2
        echo "       Make sure helicityEfficiencyToyModel.cxx and"   >&2
        echo "       plotHelicityEfficiency.cxx are in: ${SCRIPT_DIR}" >&2
        exit 1
    fi
done

# Timestamp for the start of this batch
BATCH_START=$(date '+%Y-%m-%d %H:%M:%S')
echo ""
echo "============================================================"
echo "  runHelicityToyModel.sh"
echo "  Batch started : ${BATCH_START}"
echo "  Script dir    : ${SCRIPT_DIR}"
echo "  Output dir    : ${BASE_DIR}"
echo "  Log dir       : ${LOG_DIR}"
echo "  Max parallel  : ${MAX_PARALLEL}"
echo "  Families      : ${FAMILIES_TO_RUN[*]}"
echo "  Dry run       : ${DRY_RUN}"
echo "============================================================"
echo ""

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# run_one  NAME  OUTDIR  N  BZ  PTMIN_LAM  PTMAX_LAM  RAPMAX  T  \
#               PTMIN_P  PTMIN_PI  DCAMIN_P  DCAMIN_PI  SEED
# -----------------------------------------------------------------------------
# Launches ONE generator + plotter pair.
# All numeric arguments correspond exactly to the helicityEfficiencyToyModel
# function signature (in order).
# The function is non-blocking (appends & to background).
# It returns the PID of the root process via the global LAST_PID variable.
# -----------------------------------------------------------------------------
LAST_PID=0

run_one() {
    local NAME="$1"
    local OUTDIR="${BASE_DIR}/$2"
    local N="$3"
    local BZ="$4"
    local PTMIN_LAM="$5"
    local PTMAX_LAM="$6"
    local RAPMAX="$7"
    local T="$8"
    local PTMIN_P="$9"
    local PTMIN_PI="${10}"
    local DCAMIN_P="${11}"
    local DCAMIN_PI="${12}"
    local SEED="${13}"

    local ROOT_FILE="${OUTDIR}/helicity_${NAME}.root"
    local PLOTS_DIR="${OUTDIR}/plots"
    local LOG_FILE="${LOG_DIR}/${NAME}.log"

    mkdir -p "${OUTDIR}"
    mkdir -p "${PLOTS_DIR}"

    # Build the ROOT command for the generator
    # ROOT macro arguments are passed as a comma-separated list in the
    # function call string.  Strings must be double-quoted inside single quotes.
    local GEN_CALL="${N},\"${ROOT_FILE}\",${BZ},${PTMIN_LAM},${PTMAX_LAM},${RAPMAX},${T},${PTMIN_P},${PTMIN_PI},${DCAMIN_P},${DCAMIN_PI},${SEED}"
    local PLT_CALL="\"${ROOT_FILE}\",\"${PLOTS_DIR}\""

    local CMD_GEN="${ROOT_EXE} -l -b -q '${GEN_MACRO}(${GEN_CALL})'"
    local CMD_PLT="${ROOT_EXE} -l -b -q '${PLT_MACRO}(${PLT_CALL})'"

    # Print a concise summary line regardless of dry-run mode
    printf "  %-28s  N=%-7d  Bz=%+.1f  pTmin=[%.2f,%.2f]  dca=[%.3f,%.3f]  T=%.2f\n" \
        "${NAME}" "${N}" "${BZ}" "${PTMIN_P}" "${PTMIN_PI}" \
        "${DCAMIN_P}" "${DCAMIN_PI}" "${T}"

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "    [DRY] ${CMD_GEN}"
        echo "    [DRY] ${CMD_PLT}"
        LAST_PID=0
        return 0
    fi

    # Run generator then plotter in a subshell, redirect both to the log file
    (
        echo "=== ${NAME} started at $(date '+%Y-%m-%d %H:%M:%S') ===" > "${LOG_FILE}"
        echo "Generator call: ${CMD_GEN}"                              >> "${LOG_FILE}"
        echo ""                                                         >> "${LOG_FILE}"

        # Run generator
        eval "${CMD_GEN}"  >> "${LOG_FILE}" 2>&1
        local GEN_STATUS=$?

        echo ""                                                         >> "${LOG_FILE}"
        echo "Generator exit status: ${GEN_STATUS}"                   >> "${LOG_FILE}"

        if [[ ${GEN_STATUS} -ne 0 ]]; then
            echo "ERROR: Generator failed for ${NAME}. See ${LOG_FILE}" >&2
            exit 1
        fi

        echo "Plotter call: ${CMD_PLT}"                               >> "${LOG_FILE}"
        echo ""                                                        >> "${LOG_FILE}"

        # Run plotter
        eval "${CMD_PLT}" >> "${LOG_FILE}" 2>&1
        local PLT_STATUS=$?

        echo ""                                                        >> "${LOG_FILE}"
        echo "Plotter exit status: ${PLT_STATUS}"                     >> "${LOG_FILE}"
        echo "=== ${NAME} finished at $(date '+%Y-%m-%d %H:%M:%S') ===" >> "${LOG_FILE}"

        if [[ ${PLT_STATUS} -ne 0 ]]; then
            echo "ERROR: Plotter failed for ${NAME}. See ${LOG_FILE}" >&2
            exit 1
        fi
    ) &

    LAST_PID=$!
}


# -----------------------------------------------------------------------------
# wait_for_slots  MAX
# Waits until the number of background jobs is below MAX.
# Used to implement the parallel throttle.
# -----------------------------------------------------------------------------
wait_for_slots() {
    local MAX="$1"
    while true; do
        local RUNNING
        RUNNING=$(jobs -r | wc -l)
        if [[ "${RUNNING}" -lt "${MAX}" ]]; then
            break
        fi
        sleep 2
    done
}


# -----------------------------------------------------------------------------
# run_family  FAMILY_NUM  "FAMILY_NAME"  (then calls to run_one follow)
# Prints a section header and waits for all background jobs at the end.
# -----------------------------------------------------------------------------
run_family_header() {
    local NUM="$1"
    local NAME="$2"
    echo ""
    echo "------------------------------------------------------------"
    echo "  Family ${NUM}: ${NAME}"
    echo "------------------------------------------------------------"
}

run_family_footer() {
    echo "  Waiting for Family ${1} to finish..."
    wait   # Wait for ALL background jobs from this family
    echo "  Family ${1} done."
}


# =============================================================================
# FAMILY 0: BASELINE
# Default ALICE O-O parameters.  All other families should be compared
# against this reference run.
# Parameters: N=1M, Bz=+0.5T, pT_lam=[0.3,10], |y|<0.9, T=0.30 GeV,
#             pT_min_daughters=0.15 GeV/c, DCA_min=[0.05,0.10] cm
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 0 "* ]]; then
    run_family_header 0 "BASELINE (ALICE O-O defaults)"

    #         NAME           SUBDIR         N       Bz    pTmin  pTmax  rap    T     pTp   pTpi  dcaP  dcaPi  seed
    run_one  "baseline"     "baseline"  ${DEFAULT_N}  0.5   0.3   10.0  0.9  0.30  0.150 0.150 0.050 0.100   42

    run_family_footer 0
fi


# =============================================================================
# FAMILY 1: MAGNETIC FIELD SIGN
# Swap Bz from +0.5 T to -0.5 T (same as baseline otherwise).
# Expected result: the phi* left-right asymmetry and the eta-antisymmetric
# ring proxy must flip sign exactly.  The cos(theta*) helicity asymmetry
# is field-independent (it only depends on pT thresholds) and must NOT flip.
#
# The ALICE solenoid operated at both polarities during Run 2; in Run 3 the
# nominal polarity is +0.5 T.  The data-driven DCA-based correction must be
# validated separately for each polarity.
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 1 "* ]]; then
    run_family_header 1 "MAGNETIC FIELD SIGN"

    #         NAME           SUBDIR        N       Bz     pTmin pTmax  rap    T     pTp   pTpi  dcaP  dcaPi  seed
    run_one  "field_pos"    "field_pos"  ${DEFAULT_N}  +0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100   42
    wait_for_slots ${MAX_PARALLEL}
    run_one  "field_neg"    "field_neg"  ${DEFAULT_N}  -0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100   42

    run_family_footer 1
fi


# =============================================================================
# FAMILY 2: DAUGHTER pT THRESHOLD SCAN
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
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 2 "* ]]; then
    run_family_header 2 "DAUGHTER pT THRESHOLD SCAN"

    # Zero-cut anchor: no pT requirement on daughters
    run_one  "ptcut_000"  "ptcut_000"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.001 0.001 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "ptcut_010"  "ptcut_010"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.100 0.100 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    # 0.15 is the baseline -- run with the same seed for a direct comparison
    run_one  "ptcut_015"  "ptcut_015"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "ptcut_020"  "ptcut_020"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.200 0.200 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "ptcut_030"  "ptcut_030"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.300 0.300 0.050 0.100  42

    run_family_footer 2
fi


# =============================================================================
# FAMILY 3: DCA THRESHOLD SCAN
# DCA cuts are specified as (proton DCA min, pion DCA min).
# The pion cut is typically tighter than the proton cut because the pion
# is lighter and therefore subject to larger multiple-scattering DCA smearing.
#
# dcacut_none  [0.00, 0.00]: no DCA requirement -- phi* must be flat
# dcacut_loose [0.05, 0.05]: symmetric DCA cut
# dcacut_loose [0.02, 0.05]: very loose, near the smearing limit
# dcacut_std   [0.05, 0.10]: ALICE standard for Lambda (baseline)
# dcacut_tight [0.10, 0.20]: tighter, sometimes used in central Pb-Pb
#
# Reference: ALICE strangeness PWG recommendations;
#            arXiv:1910.07678 (ALICE Lambda polarization, Run 2)
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 3 "* ]]; then
    run_family_header 3 "DCA THRESHOLD SCAN"

    # Zero anchor: phi* must be flat
    run_one  "dcacut_none"   "dcacut_none"   ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.000 0.000  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "dcacut_loose"  "dcacut_loose"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.020 0.050  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "dcacut_loose_symmetric"  "dcacut_loose_symmetric"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.050 0.050  42
    wait_for_slots ${MAX_PARALLEL}

    # Standard baseline values (same as baseline run)
    run_one  "dcacut_std"    "dcacut_std"    ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "dcacut_tight"  "dcacut_tight"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150  0.100 0.200  42

    run_family_footer 3
fi


# =============================================================================
# FAMILY 4: LAMBDA KINEMATIC WINDOW
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
#
# Reference: PhysRevC.109.014905 (Serone et al., the 3rd ring paper)
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 4 "* ]]; then
    run_family_header 4 "LAMBDA KINEMATIC WINDOW"

    #         NAME             SUBDIR          N       Bz   pTmin pTmax  rap    T     pTp   pTpi  dcaP  dcaPi  seed
    run_one  "window_full"   "window_full"   ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "window_ring"   "window_ring"   ${DEFAULT_N}  0.5  0.5   1.5  0.5  0.30  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "window_hard"   "window_hard"   ${DEFAULT_N}  0.5  1.5   4.0  0.5  0.30  0.150 0.150 0.050 0.100  42

    run_family_footer 4
fi


# =============================================================================
# FAMILY 5: THERMAL TEMPERATURE (Lambda pT spectrum)
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
#
# Reference: ALICE, Phys.Lett.B 728 (2014) 25; ALICE, JHEP 07 (2015) 116
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 5 "* ]]; then
    run_family_header 5 "THERMAL TEMPERATURE SCAN"

    #         NAME         SUBDIR       N       Bz   pTmin pTmax  rap    T     pTp   pTpi  dcaP  dcaPi  seed
    run_one  "temp_025"  "temp_025"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.25  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "temp_030"  "temp_030"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "temp_035"  "temp_035"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.35  0.150 0.150 0.050 0.100  42
    wait_for_slots ${MAX_PARALLEL}

    run_one  "temp_045"  "temp_045"  ${DEFAULT_N}  0.5  0.3  10.0  0.9  0.45  0.150 0.150 0.050 0.100  42

    run_family_footer 5
fi


# =============================================================================
# FAMILY 6: STATISTICS SCALING
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
# =============================================================================
if [[ " ${FAMILIES_TO_RUN[*]} " == *" 6 "* ]]; then
    run_family_header 6 "STATISTICS SCALING"

    #         NAME          SUBDIR       N       Bz   pTmin pTmax  rap    T     pTp   pTpi  dcaP  dcaPi  seed
    run_one  "stats_100k"  "stats_100k"  100000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  101
    wait_for_slots ${MAX_PARALLEL}

    run_one  "stats_500k"  "stats_500k"  500000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  201
    wait_for_slots ${MAX_PARALLEL}

    run_one  "stats_2M"    "stats_2M"    2000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  301
    wait_for_slots ${MAX_PARALLEL}

    run_one  "stats_5M"    "stats_5M"    5000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  401

    run_one  "stats_50M"   "stats_50M"   50000000  0.5  0.3  10.0  0.9  0.30  0.150 0.150 0.050 0.100  401

    run_family_footer 6
fi


# =============================================================================
# FINAL SUMMARY
# =============================================================================
BATCH_END=$(date '+%Y-%m-%d %H:%M:%S')

echo ""
echo "============================================================"
echo "  All requested families complete."
echo "  Batch started : ${BATCH_START}"
echo "  Batch finished: ${BATCH_END}"
echo ""
echo "  Output tree:"
echo "    ROOT files  : ${BASE_DIR}/<run_name>/helicity_<run_name>.root"
echo "    Plot PDFs   : ${BASE_DIR}/<run_name>/plots/*.pdf"
echo "    Logs        : ${LOG_DIR}/<run_name>.log"
echo ""
echo "  To inspect logs:"
echo "    cat ${LOG_DIR}/<run_name>.log"
echo "    grep 'R_proxy\|ERROR\|exit' ${LOG_DIR}/*.log"
echo ""
echo "  To open a result in ROOT:"
echo "    root -l ${BASE_DIR}/baseline/helicity_baseline.root"
echo "============================================================"
echo ""
