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
#       Varies etaMaxDetector (and etaMin detector symetrically) from 0.5 to 3.0.
#       Shows how the fake signal magnitude changes with the detector acceptance window.
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
#   ./runHelicityToyModel.sh --multiply 2       # Split each job into 2 NUMA-pinned shards
#   ./runHelicityToyModel.sh -m 3               # Same with short flag
#   ./runHelicityToyModel.sh --dry-run          # Print jobs without executing
#   ./runHelicityToyModel.sh --list             # List all registered jobs and exit
#   ./runHelicityToyModel.sh --plot             # Skip generation, run plotter only
#   ./runHelicityToyModel.sh --help             # Show this help
#
# NUMA SHARDING (--multiply / -m)
# --------------------------------
# On a multi-socket machine, --multiply M splits each registered job into M
# independent generator sub-processes (shards), each processing N/M events
# with a unique seed, and each pinned to a distinct NUMA node via:
#   numactl --cpunodebind=K --membind=K
# This keeps compute and memory on the same die, avoiding the cross-socket
# latency of the AMD Infinity Fabric (~2x memory latency vs local DRAM).
# After all M shards finish, hadd merges their ROOT files into the final output
# and the shard files are deleted. The plotter then runs once on the merged file.
#
# CORRECTNESS: hadd is exact for all generator output histogram types:
#   TH1D / TH2D (count histograms): hadd sums bin contents and Sumw2. EXACT.
#   TProfile (written via flat accumulators + SetBinContent/SetBinEntries):
#     hadd sums fArray (sum_wy) and fBinEntries (sum_w) per bin, giving
#     GetBinContent = (sum_wy_A + sum_wy_B) / (sum_w_A + sum_w_B), the
#     correct weighted mean. ROOT's internal Add() also rebuilds the global
#     statistics cache (fTsumwy etc.), so the statsbox is correct without
#     any additional ResetStats() call after the merge. EXACT.
#
# SEED POLICY (ensures statistical independence between shards):
#   SEED != 0  ->  shard i gets seed = SEED * 1000 + i.
#                  Jobs with neighboring base seeds never collide: they differ
#                  by >= 1000, which is always > MULTIPLY.
#   SEED == 0  ->  M seeds are drawn from /dev/urandom so that near-simultaneous
#                  processes do not risk identical TRandom3(0) auto-seeds.
#
# AMD EPYC 9634 dual-socket layout (this machine):
#   NUMA node 0: CPUs 0-83,   168-251  (84 physical cores)
#   NUMA node 1: CPUs 84-167, 252-335  (84 physical cores)
#   With ~79 registered jobs and --multiply 2: 158 generator processes, filling
#   all 168 physical cores simultaneously. Wall time ~ time of one shard.
#   With --multiply 1 (default): same behavior as before, no NUMA pinning.
#
# DOWNSIDE: with MULTIPLY > 1 the number of concurrent processes increases by
# factor M. At 79 jobs x M=2 = 158 processes on 168 physical cores this is
# near-ideal. At M=3 you get 237 processes on 168 physical cores, which starts
# to over-subscribe via hyperthreading (each physical core has 2 logical CPUs).
# HT can help for memory-bound code but our generator is compute-bound, so
# M=2 is the recommended value for this machine.
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
#     6_EtaSym/            eta_050/ eta_070/ eta_090/ eta_120/ eta_150/
#                          eta_200/ eta_300/
#     7_Temperature/       temp_020/ ... temp_050/
#     8_KinWindow/         win_inclusive/ win_vsoft/ win_soft/ win_ring/
#                          win_mid/ win_hard/ win_vhard/
#     9_RealisticAlice/    alice_loose/ alice_std/ alice_std_neg/
#                          alice_tight/ alice_ring/
#     10_BField/           fieldCuts_vlow/ fieldCuts_b050/ fieldCuts_b100/
#                          fieldCuts_b050neg/
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

# Maximum concurrent ROOT processes (one per registered job).
# Default: 80 (all jobs at once on a 192-core machine).
# Set to 0 to let GNU parallel use all available cores (parallel -j 0),
# or to any positive integer to cap concurrency explicitly.
# NOTE: with MULTIPLY > 1 each "job slot" launches MULTIPLY sub-processes
# internally, so the true peak process count is MAX_PARALLEL * MULTIPLY.
# At 80 jobs * MULTIPLY=2 = 160 processes this saturates the 168 physical
# cores of this machine near-perfectly.
MAX_PARALLEL=80

# NUMA sharding multiplier.
# MULTIPLY=1: no sharding (default, identical to previous behavior).
# MULTIPLY=2: each job spawns 2 shards on NUMA nodes 0 and 1 respectively,
#             then merges with hadd. Recommended for this dual-socket machine.
# MULTIPLY=M: cycles through N_NUMA_NODES nodes (detected automatically below).
# Override from the command line with --multiply M or -m M.
MULTIPLY=1

# NUMA topology: detected automatically from numactl --hardware.
# Falls back gracefully (no pinning, no error) if numactl is not installed.
NUMA_AVAILABLE=0
N_NUMA_NODES=1
if command -v numactl &>/dev/null; then
    NUMA_AVAILABLE=1
    # Count lines of the form "node N cpus: ..." emitted by numactl --hardware
    _DETECTED=$(numactl --hardware 2>/dev/null \
                | grep -c '^node [0-9][0-9]* cpus:' || true)
    if [[ "${_DETECTED}" -ge 1 ]]; then
        N_NUMA_NODES="${_DETECTED}"
    fi
    unset _DETECTED
fi

# hadd is ROOT's histogram merge tool, required when MULTIPLY > 1.
# It ships with ROOT and should be in PATH if root-config is available.
HADD_EXE="hadd"

# Default Lambda count per run.
DEFAULT_N=360000000000 # Takes about 4 days after the optimized generator cxx, 11.5 days with the older code
# DEFAULT_N=180000000000 # 180.000.000.000 Lambdas would take 10.5 days previously. With refactor it took 4.3 days!
# DEFAULT_N=10000000000 # 10.000.000.000 Lambdas takes about 12 hours on newer, heavier, code versions (After PGO and refactoring, went down to 6 hours)
# DEFAULT_N=100000000 # For testing only (Refactoring went from ~5 minutes per run to ~3.5 minutes. After PGO, this changed little)
# DEFAULT_N=10000000 # For testing only (Refactoring went from ~40 s per run to ~25s per run)

# ROOT executable (set to full path if not in PATH, e.g. /opt/root/bin/root)
ROOT_EXE="root"

# Generator macro and plotting macro (relative to SCRIPT_DIR)
# GEN_MACRO="${SCRIPT_DIR}/helicityEfficiencyToyModel_OptimizationRefactor.cxx"
GEN_MACRO="${SCRIPT_DIR}/helicityEfficiencyToyModel.cxx"
PLT_MACRO="${SCRIPT_DIR}/plotHelicityEfficiency.cxx"


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

DRY_RUN=0
LIST_ONLY=0
PLOT_ONLY=0
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
        --plot|-p)
            PLOT_ONLY=1
            shift
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
        --multiply|-m)
            if [[ -z "${2:-}" || ! "${2}" =~ ^[0-9]+$ || "${2}" -lt 1 ]]; then
                echo "ERROR: --multiply requires a positive integer argument (>= 1)." >&2
                exit 1
            fi
            MULTIPLY="$2"
            shift 2
            # Validate hadd availability when sharding is requested
            if [[ "${MULTIPLY}" -gt 1 ]] && ! command -v "${HADD_EXE}" &>/dev/null; then
                echo "ERROR: --multiply ${MULTIPLY} requires 'hadd' (ROOT merge tool) in PATH." >&2
                echo "       Ensure ROOT is set up: source /path/to/root/bin/thisroot.sh" >&2
                exit 1
            fi
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
    FAMILIES_TO_RUN=(0 1 2 3 4 5 6 7 8 9 10 11 12)
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
# COMPILATION
# =============================================================================
# EXE_BIN="${SCRIPT_DIR}/helicityEfficiencyToyModel_OptimizationRefactor.exe"
EXE_BIN="${SCRIPT_DIR}/helicityEfficiencyToyModel.exe"
COMPILE_WARN_FLAGS="-Wall -Wextra -Wpedantic -Wshadow" # For bug catching
ROOT_AWARE_FLAGS="-DNDEBUG" # Strips out all assert() checks. ROOT headers respect this flag!
OPTIMIZATION_FLAGS="-O3 -march=native -fno-math-errno -flto -fno-trapping-math" # Using the fno-math-errno to disable the errno verification of passing negative numbers to sqrts and the such.
                                                                                # My code should be safe, so this will never be used!
                                                                                # Also introduced -fno-trapping-math as there is no explicit usage of NaNs
                                                                                # nor any expected unsafe divisions by zero appearing anywhere in the code!
export EXE_BIN # Export so subshells and GNU parallel can see it

if [[ ${DRY_RUN} -eq 0 && ${LIST_ONLY} -eq 0 && ${PLOT_ONLY} -eq 0 ]]; then
    echo "============================================================"
    echo " Compiling generator with g++ (AOT (Ahead-Of-Time) optimization)..."
    if g++ -std=c++20 -pipe ${COMPILE_WARN_FLAGS} ${ROOT_AWARE_FLAGS} ${OPTIMIZATION_FLAGS} "${GEN_MACRO}" -o "${EXE_BIN}" $(root-config --cflags --libs); then
        echo " Compilation successful: ${EXE_BIN}"
    else
        echo " ERROR: Compilation failed!" >&2
        exit 1
    fi
    echo "============================================================"
    echo ""
fi

# Not using PGO -- Each Model Family has a very different profile of which ifs and elses it triggers, so PGO may actually worsen the binary performance:
# # =============================================================================
# # COMPILATION WITH PROFILE-GUIDED OPTIMIZATION (PGO)
# # =============================================================================
# EXE_BIN="${SCRIPT_DIR}/helicityEfficiencyToyModel_OptimizationRefactor.exe"
# COMPILE_WARN_FLAGS="-Wall -Wextra -Wpedantic -Wshadow" # For bug catching
# ROOT_AWARE_FLAGS="-DNDEBUG" # Strips out all assert() checks.
# BASE_OPT_FLAGS="-O3 -march=native -fno-math-errno -flto -fno-trapping-math" 
# export EXE_BIN # Export so subshells and GNU parallel can see it

# if [[ ${DRY_RUN} -eq 0 && ${LIST_ONLY} -eq 0 && ${PLOT_ONLY} -eq 0 ]]; then
#     echo "============================================================"
#     echo " Compiling generator with g++ (Profile-Guided Optimization)..."
    
#     # PHASE 1: Build the instrumented executable
#     echo " [1/3] Building instrumented binary for profiling..."
#     if ! g++ -std=c++20 -pipe ${COMPILE_WARN_FLAGS} ${ROOT_AWARE_FLAGS} ${BASE_OPT_FLAGS} -fprofile-generate "${GEN_MACRO}" -o "${EXE_BIN}" $(root-config --cflags --libs); then
#         echo " ERROR: PGO Instrumentation Compilation failed!" >&2
#         exit 1
#     fi
    
#     # PHASE 2: Run the training workload
#     echo " [2/3] Running PGO training workload (baseline, 10M events)..."
#     PGO_TEMP_ROOT="${SCRIPT_DIR}/pgo_training_temp.root"
    
#     # Passing 10M events and a dummy root file. 
#     # Because my C++ argv parser cascades, the rest will naturally fall back to the baseline defaults.
#     # if ! "${EXE_BIN}" 10000000 "${PGO_TEMP_ROOT}"; then
#     if ! "${EXE_BIN}" 10000000 "${PGO_TEMP_ROOT}" > /dev/null; then
#         echo " ERROR: PGO Training run failed!" >&2
#         exit 1
#     fi
    
#     # PHASE 3: Recompile using the generated profile data (.gcda files)
#     echo " [3/3] Recompiling optimized binary using profile data..."
#     # Note: -fprofile-correction is added to ensure smooth combination with LTO
#     if ! g++ -std=c++20 -pipe ${COMPILE_WARN_FLAGS} ${ROOT_AWARE_FLAGS} ${BASE_OPT_FLAGS} -fprofile-use -fprofile-correction "${GEN_MACRO}" -o "${EXE_BIN}" $(root-config --cflags --libs); then
#         echo " ERROR: PGO Optimized Compilation failed!" >&2
#         exit 1
#     fi
    
#     # CLEANUP: Remove the profiling artifacts and temporary ROOT file
#     echo " Cleaning up PGO artifacts..."
#     # The .gcda files are dumped in the working directory by GCC
#     rm -f *.gcda "${PGO_TEMP_ROOT}"
    
#     echo " PGO Compilation successful: ${EXE_BIN}"
#     echo "============================================================"
#     echo ""
# fi


# =============================================================================
# JOB QUEUE
# Each call to register_job() appends one record to JOB_QUEUE[].
# Format: colon-delimited, 14 fields matching helicityEfficiencyToyModel args:
#   NAME:SUBDIR:N:BZ:PTMIN_LAM:PTMAX_LAM:RAPMAX:ETAMIN:ETAMAX:T:PTMIN_P:PTMIN_PI:DCAMIN_P:DCAMIN_PI:SEED
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
#   $8  ETAMAX     Daughter maximum detectable eta
#   $9  ETAMIN     Daughter minimum detectable eta
#   $10 T          Boltzmann temperature [GeV]
#   $11 PTMIN_P    proton min pT [GeV/c]
#   $12 PTMIN_PI   pion min pT [GeV/c]
#   $13 DCAMIN_P   proton min DCA_xy [cm]
#   $14 DCAMIN_PI  pion min DCA_xy [cm]
#   $15 SEED       TRandom3 seed
register_job() {
    local PACKED="$1:$2:$3:$4:$5:$6:$7:$8:$9:${10}:${11}:${12}:${13}:${14}:${15}"
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
    # Unpack the 15 packed fields (colon-delimited, produced by register_job)
    IFS=':' read -r \
        NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX ETAMIN ETAMAX T \
        PTMIN_P PTMIN_PI DCAMIN_P DCAMIN_PI SEED \
        <<< "$1"

    local OUTDIR="${BASE_DIR}/${SUBDIR}"
    local ROOT_FILE="${OUTDIR}/helicity_${NAME}.root"
    local LOG_FILE="${LOG_DIR}/${NAME}.log"
    # Plotter still uses ROOT Cling (cheap, I/O-bound, so no need to pre-compile this cxx)
    local PLT_CALL="\"${ROOT_FILE}\""

    mkdir -p "${OUTDIR}"

    # Write log header
    {
        echo "=== ${NAME} started at $(date '+%Y-%m-%d %H:%M:%S') ==="
        echo "Host      : $(hostname)"
        echo "PID       : $$"
        echo "Binary    : ${EXE_BIN}"
        echo "Plotter   : ${ROOT_EXE} -l -b -q PLT_MACRO(${ROOT_FILE})"
        echo "Multiply  : ${MULTIPLY}  (NUMA nodes available: ${N_NUMA_NODES})"
        echo "NUMA pin  : $([ "${NUMA_AVAILABLE}" -eq 1 ] && echo 'yes (numactl)' || echo 'no (numactl not found)')"
        echo ""
    } > "${LOG_FILE}"

    local T0; T0=$(date +%s)

    # -------------------------------------------------------------------------
    # CONDITIONAL GENERATOR EXECUTION (if you only update the plotting function, then you may not want to re-generate the data)
    # -------------------------------------------------------------------------
    if [[ ${PLOT_ONLY} -eq 0 ]]; then

        if [[ "${MULTIPLY}" -le 1 ]]; then
            # -----------------------------------------------------------------
            # STANDARD single-process path  (MULTIPLY = 1)
            # Identical to the original behavior; no NUMA pinning, no sharding.
            # Build argument array for the compiled generator (safely handles spaces in ROOT_FILE)
            # Notice this is different from before, where we parsed a giant string
            # -----------------------------------------------------------------
            local GEN_ARGS=("${N}" "${ROOT_FILE}" "${BZ}"
                            "${PTMIN_LAM}" "${PTMAX_LAM}" "${RAPMAX}"
                            "${ETAMIN}" "${ETAMAX}" "${T}"
                            "${PTMIN_P}" "${PTMIN_PI}"
                            "${DCAMIN_P}" "${DCAMIN_PI}" "${SEED}")
            echo "Generator : ${EXE_BIN} ${GEN_ARGS[*]}" >> "${LOG_FILE}"
            echo "" >> "${LOG_FILE}"

            if "${EXE_BIN}" "${GEN_ARGS[@]}" >> "${LOG_FILE}" 2>&1; then
                local T1; T1=$(date +%s)
                echo "Generator OK  ($(( T1 - T0 )) s)" >> "${LOG_FILE}"
            else
                local STATUS=$?
                echo "Generator FAILED (exit ${STATUS})" >> "${LOG_FILE}"
                echo "ERROR [${NAME}]: generator failed (exit ${STATUS}). Log: ${LOG_FILE}" >&2
                return "${STATUS}"
            fi

        else
            # -----------------------------------------------------------------
            # NUMA-SHARDED multi-process path  (MULTIPLY > 1)
            #
            # Strategy:
            #   1. Divide N events among MULTIPLY shards. The last shard gets
            #      the remainder so that sum(shard_N) == N exactly.
            #   2. Assign each shard a unique seed (see SEED POLICY below).
            #   3. Launch each shard as a background process pinned to its
            #      NUMA node (numactl --cpunodebind=K --membind=K).
            #   4. Wait for all shards; fail the job if any shard fails.
            #   5. Merge shard ROOT files into the final output via hadd.
            #   6. Delete shard files.
            #
            # SEED POLICY (ensures independent random streams):
            #   SEED != 0 -> shard i gets SEED * 1000 + i.
            #     Neighboring job seeds (SEED and SEED+1) differ by 1000 in
            #     their base, so shard seeds never collide between jobs as long
            #     as MULTIPLY < 1000 (easily satisfied).
            #   SEED == 0 -> draw a random base from /dev/urandom so that
            #     near-simultaneous shards with seed=0 do not all auto-seed to
            #     the same value from the system clock. Each shard then gets
            #     base + i, guaranteeing distinct seeds.
            #
            # CORRECTNESS OF hadd FOR OUR HISTOGRAM TYPES:
            #   TH1D / TH2D:  hadd sums bin contents and Sumw2. EXACT.
            #   TProfile (via flat accumulators synced with SetBinContent /
            #   SetBinEntries):  hadd sums fArray (sum_wy) and fBinEntries
            #   (sum_w) per bin. GetBinContent = total_sum_wy / total_sum_w,
            #   the correct weighted mean. ROOT's Add() (called by hadd) also
            #   rebuilds the global statistics cache, so the statsbox is correct
            #   without a further ResetStats() call. EXACT.
            # -----------------------------------------------------------------

            # Per-shard event count (integer division; last shard takes remainder)
            local N_PER_SHARD=$(( N / MULTIPLY ))
            local N_LAST=$(( N - N_PER_SHARD * (MULTIPLY - 1) ))

            # Base seed for shard generation
            local BASE_SEED
            if [[ "${SEED}" -ne 0 ]]; then
                BASE_SEED=$(( SEED * 1000 ))
            else
                # /dev/urandom gives 4 random bytes -> unsigned 32-bit integer
                BASE_SEED=$(od -An -N4 -tu4 /dev/urandom | tr -d ' \n')
            fi

            echo "Sharding  : ${MULTIPLY} shards x ~${N_PER_SHARD} events (last shard: ${N_LAST})" \
                >> "${LOG_FILE}"
            echo "Base seed : ${BASE_SEED}" >> "${LOG_FILE}"
            echo "" >> "${LOG_FILE}"

            # Arrays to track shard PIDs and output files
            local -a SHARD_PIDS=()
            local -a SHARD_ROOTS=()

            local SHARD_IDX
            for (( SHARD_IDX = 0; SHARD_IDX < MULTIPLY; SHARD_IDX++ )); do

                local SHARD_ROOT="${OUTDIR}/shard_${SHARD_IDX}_${NAME}.root"
                SHARD_ROOTS+=("${SHARD_ROOT}")

                # Last shard absorbs the remainder so sum == N exactly
                local SHARD_N
                if (( SHARD_IDX == MULTIPLY - 1 )); then
                    SHARD_N=${N_LAST}
                else
                    SHARD_N=${N_PER_SHARD}
                fi

                # Unique seed for this shard
                local SHARD_SEED=$(( BASE_SEED + SHARD_IDX ))

                # NUMA node: cycle through available nodes (0, 1, 0, 1, ...)
                local NUMA_NODE=$(( SHARD_IDX % N_NUMA_NODES ))

                local GEN_ARGS=("${SHARD_N}" "${SHARD_ROOT}" "${BZ}"
                                "${PTMIN_LAM}" "${PTMAX_LAM}" "${RAPMAX}"
                                "${ETAMIN}" "${ETAMAX}" "${T}"
                                "${PTMIN_P}" "${PTMIN_PI}"
                                "${DCAMIN_P}" "${DCAMIN_PI}" "${SHARD_SEED}")

                echo "  Shard ${SHARD_IDX}: N=${SHARD_N}  seed=${SHARD_SEED}  NUMA=${NUMA_NODE}  -> $(basename "${SHARD_ROOT}")" \
                    >> "${LOG_FILE}"

                # Launch with NUMA binding if numactl is available;
                # otherwise launch without pinning (still uses separate seeds)
                if [[ "${NUMA_AVAILABLE}" -eq 1 ]]; then
                    numactl \
                        --cpunodebind="${NUMA_NODE}" \
                        --membind="${NUMA_NODE}" \
                        "${EXE_BIN}" "${GEN_ARGS[@]}" \
                        >> "${LOG_FILE}" 2>&1 &
                else
                    "${EXE_BIN}" "${GEN_ARGS[@]}" >> "${LOG_FILE}" 2>&1 &
                fi
                SHARD_PIDS+=($!)

            done  # end shard launch loop

            # Wait for all shards and collect exit statuses individually
            echo "" >> "${LOG_FILE}"
            echo "  Waiting for ${MULTIPLY} shards..." >> "${LOG_FILE}"

            local SHARD_FAIL=0
            for (( SHARD_IDX = 0; SHARD_IDX < MULTIPLY; SHARD_IDX++ )); do
                local SPID="${SHARD_PIDS[${SHARD_IDX}]}"
                if wait "${SPID}"; then
                    echo "  Shard ${SHARD_IDX} OK   (PID ${SPID})" >> "${LOG_FILE}"
                else
                    local SSTAT=$?
                    echo "  Shard ${SHARD_IDX} FAIL (PID ${SPID}, exit ${SSTAT})" >> "${LOG_FILE}"
                    SHARD_FAIL=1
                fi
            done

            if [[ ${SHARD_FAIL} -ne 0 ]]; then
                echo "" >> "${LOG_FILE}"
                echo "Generator FAILED: one or more shards exited non-zero" >> "${LOG_FILE}"
                echo "ERROR [${NAME}]: shard(s) failed. Log: ${LOG_FILE}" >&2
                return 1
            fi

            local T1; T1=$(date +%s)
            echo "" >> "${LOG_FILE}"
            echo "All ${MULTIPLY} shards OK  ($(( T1 - T0 )) s wall time)" >> "${LOG_FILE}"
            echo "" >> "${LOG_FILE}"

            # Merge all shard ROOT files into the final output file
            echo "  Running: ${HADD_EXE} -f ${ROOT_FILE} ${SHARD_ROOTS[*]}" \
                >> "${LOG_FILE}"
            if "${HADD_EXE}" -f "${ROOT_FILE}" "${SHARD_ROOTS[@]}" \
                    >> "${LOG_FILE}" 2>&1; then
                echo "  hadd OK  -> ${ROOT_FILE}" >> "${LOG_FILE}"
            else
                local HADD_STATUS=$?
                echo "  hadd FAILED (exit ${HADD_STATUS})" >> "${LOG_FILE}"
                echo "ERROR [${NAME}]: hadd failed. Log: ${LOG_FILE}" >&2
                return "${HADD_STATUS}"
            fi

            # Remove per-shard temporary files now that merge is complete
            local SR
            for SR in "${SHARD_ROOTS[@]}"; do
                rm -f "${SR}"
            done
            echo "  Shard files removed." >> "${LOG_FILE}"

        fi  # end MULTIPLY branch

        echo "" >> "${LOG_FILE}"

    fi  # end PLOT_ONLY guard

    echo "" >> "${LOG_FILE}"

    # -------------------------------------------------------------------------
    # PLOTTER  (always runs -- it is cheap and I/O-bound)
    # Runs on the final merged ROOT file regardless of MULTIPLY value.
    # -------------------------------------------------------------------------
    local T2; T2=$(date +%s)
    if eval "${ROOT_EXE} -l -b -q '${PLT_MACRO}(${PLT_CALL})'" \
            >> "${LOG_FILE}" 2>&1; then
        local T3; T3=$(date +%s)
        echo "Plotter  OK  ($(( T3 - T2 )) s)" >> "${LOG_FILE}"
    else
        local PLT_STATUS=$?
        echo "Plotter FAILED (exit ${PLT_STATUS})" >> "${LOG_FILE}"
        echo "ERROR [${NAME}]: plotter failed (exit ${PLT_STATUS}). Log: ${LOG_FILE}" >&2
        return "${PLT_STATUS}"
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
export BASE_DIR LOG_DIR ROOT_EXE EXE_BIN GEN_MACRO PLT_MACRO PLOT_ONLY
export MULTIPLY NUMA_AVAILABLE N_NUMA_NODES HADD_EXE


# =============================================================================
# FAMILY DEFINITIONS
# Each block is guarded by family_in_scope(); it only calls register_job.
# The order of register_job calls within a block is irrelevant because all
# jobs are dispatched simultaneously by the pool at the end.
#
# Column guide (all register_job calls):
#   NAME  SUBDIR  N  Bz  pTminLam  pTmaxLam  rapMax  etaMin  etaMax  T  pTp  pTpi  dcaP  dcaPi  seed
#
# Shared baseline parameter values (always used unless the family varies them):
#   Bz         = +0.5 T
#   pTminLam   = 0.0 GeV/c
#   pTmaxLam   = 10.0 GeV/c
#   rapMax     = 4.0          (wide enough to include all daughters inside etaMax=0.9)
#   etaMin     = -0.9          (ALICE inner barrel)
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

    #              NAME        SUBDIR                N              Bz    pTminL  pTmaxL  rap   etaMin   etaMax  T     pTp   pTpi  dcaP  dcaPi seed
    register_job  "baseline"  "0_Baseline/baseline"  ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9     0.9     0.30  0.0   0.0   0.0   0.0   0
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

    #               NAME                   SUBDIR                        N              Bz    pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp   pTpi  dcaP   dcaPi  seed
    register_job   "dca_p_only_005"       "1_AsymDCA/dca_p_only_005"     ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.000  0
    register_job   "dca_pi_only_005"      "1_AsymDCA/dca_pi_only_005"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.000  0.050  0
    register_job   "dca_pi_only_010"      "1_AsymDCA/dca_pi_only_010"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.000  0.100  0
    register_job   "dca_asym_std"         "1_AsymDCA/dca_asym_std"       ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.100  0
    register_job   "dca_asym_wide"        "1_AsymDCA/dca_asym_wide"      ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.200  0
    register_job   "dca_asym_rev"         "1_AsymDCA/dca_asym_rev"       ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.100  0.050  0
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

    #               NAME              SUBDIR                   N             Bz    pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp   pTpi  dcaP   dcaPi  seed
    register_job   "dca_sym_005"     "2_SymDCA/dca_sym_005"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "dca_sym_010"     "2_SymDCA/dca_sym_010"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.100  0.100  0
    register_job   "dca_sym_015"     "2_SymDCA/dca_sym_015"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.150  0.150  0
    register_job   "dca_sym_020"     "2_SymDCA/dca_sym_020"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.200  0.200  0
    register_job   "dca_sym_030"     "2_SymDCA/dca_sym_030"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.300  0.300  0
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

    #               NAME               SUBDIR                     N             Bz       pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp   pTpi  dcaP   dcaPi  seed
    register_job   "field_vlow"       "3_BField/field_vlow"       ${DEFAULT_N}  0.0001   0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b010"       "3_BField/field_b010"       ${DEFAULT_N}  0.10     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b020"       "3_BField/field_b020"       ${DEFAULT_N}  0.20     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b030"       "3_BField/field_b030"       ${DEFAULT_N}  0.30     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b050"       "3_BField/field_b050"       ${DEFAULT_N}  0.50     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b075"       "3_BField/field_b075"       ${DEFAULT_N}  0.75     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b100"       "3_BField/field_b100"       ${DEFAULT_N}  1.00     0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
    register_job   "field_b050neg"    "3_BField/field_b050neg"    ${DEFAULT_N}  -0.50    0.0     10.0    4.0   -0.9    0.9     0.30  0.0   0.0   0.050  0.050  0
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

    #               NAME                SUBDIR                       N             Bz   pTminL  pTmaxL  rap   etaMin   etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "lam_ptmin_000"     "4_LamPtMin/lam_ptmin_000"    ${DEFAULT_N}  0.5  0.000   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_030"     "4_LamPtMin/lam_ptmin_030"    ${DEFAULT_N}  0.5  0.300   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_050"     "4_LamPtMin/lam_ptmin_050"    ${DEFAULT_N}  0.5  0.500   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_100"     "4_LamPtMin/lam_ptmin_100"    ${DEFAULT_N}  0.5  1.000   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_150"     "4_LamPtMin/lam_ptmin_150"    ${DEFAULT_N}  0.5  1.500   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "lam_ptmin_200"     "4_LamPtMin/lam_ptmin_200"    ${DEFAULT_N}  0.5  2.000   10.0    4.0   -0.9     0.9     0.30  0.150  0.150  0.050  0.050  0
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
    #               NAME              SUBDIR                      N             Bz   pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp    pTpi   dcaP  dcaPi  seed
    register_job   "pt_sym_005"      "5_DaughterPt/pt_sym_005"    ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.050  0.050  0.0   0.0   0
    register_job   "pt_sym_010"      "5_DaughterPt/pt_sym_010"    ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.100  0.100  0.0   0.0   0
    register_job   "pt_sym_150"      "5_DaughterPt/pt_sym_150"    ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.150  0.150  0.0   0.0   0
    register_job   "pt_sym_200"      "5_DaughterPt/pt_sym_200"    ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.200  0.200  0.0   0.0   0
    register_job   "pt_sym_300"      "5_DaughterPt/pt_sym_300"    ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.300  0.300  0.0   0.0   0
    # -- Pion-only cuts (purest HEE isolation) --
    register_job   "pt_pi_010"       "5_DaughterPt/pt_pi_010"     ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.000  0.100  0.0   0.0   0
    register_job   "pt_pi_150"       "5_DaughterPt/pt_pi_150"     ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.000  0.150  0.0   0.0   0
    register_job   "pt_pi_200"       "5_DaughterPt/pt_pi_200"     ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.000  0.200  0.0   0.0   0
    # -- Proton-only cut (control: expect weaker HEE than pion-only) --
    register_job   "pt_p_150"        "5_DaughterPt/pt_p_150"      ${DEFAULT_N}  0.5  0.0     10.0    4.0   -0.9    0.9     0.30  0.150  0.000  0.0   0.0   0
fi

# -----------------------------------------------------------------------------
# FAMILY 6: ETA ACCEPTANCE WINDOW  [HEE+AEE vs detector acceptance]
# Standard-ish ALICE cuts (pTp=pTpi=0.15, dcaP=0.05, dcaPi=0.05).
# Now specifically a symmetric eta cut study
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

    #               NAME           SUBDIR               N              Bz    pTminL  pTmaxL  rap   etaMin   etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "eta_050"      "6_EtaSym/eta_050"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.50    0.50    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_070"      "6_EtaSym/eta_070"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.70    0.70    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_090"      "6_EtaSym/eta_090"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -0.90    0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_120"      "6_EtaSym/eta_120"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -1.20    1.20    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_150"      "6_EtaSym/eta_150"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -1.50    1.50    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_200"      "6_EtaSym/eta_200"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -2.00    2.00    0.30  0.150  0.150  0.050  0.050  0
    register_job   "eta_300"      "6_EtaSym/eta_300"    ${DEFAULT_N}   0.5   0.0     10.0    4.0   -3.00    3.00    0.30  0.150  0.150  0.050  0.050  0
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

    #               NAME           SUBDIR                     N             Bz    pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "temp_020"     "7_Temperature/temp_020"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.20  0.150  0.150  0.050  0.050  0
    register_job   "temp_025"     "7_Temperature/temp_025"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.25  0.150  0.150  0.050  0.050  0
    register_job   "temp_030"     "7_Temperature/temp_030"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "temp_035"     "7_Temperature/temp_035"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.35  0.150  0.150  0.050  0.050  0
    register_job   "temp_040"     "7_Temperature/temp_040"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.40  0.150  0.150  0.050  0.050  0
    register_job   "temp_050"     "7_Temperature/temp_050"    ${DEFAULT_N}  0.5   0.0     10.0    4.0   -0.9    0.9     0.50  0.150  0.150  0.050  0.050  0
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

    #               NAME                SUBDIR                      N             Bz   pTminL  pTmaxL  rap  etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "win_inclusive"     "8_KinWindow/win_inclusive"  ${DEFAULT_N}  0.5  0.000   10.000  4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_vsoft"         "8_KinWindow/win_vsoft"      ${DEFAULT_N}  0.5  0.000   0.500   4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_soft"          "8_KinWindow/win_soft"       ${DEFAULT_N}  0.5  0.000   1.000   4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_ring"          "8_KinWindow/win_ring"       ${DEFAULT_N}  0.5  0.500   1.500   4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_mid"           "8_KinWindow/win_mid"        ${DEFAULT_N}  0.5  1.000   3.000   4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_hard"          "8_KinWindow/win_hard"       ${DEFAULT_N}  0.5  1.500   4.000   4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "win_vhard"         "8_KinWindow/win_vhard"      ${DEFAULT_N}  0.5  3.000   10.000  4.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
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

    #               NAME                  SUBDIR                             N             Bz     pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "alice_loose"         "9_RealisticAlice/alice_loose"      ${DEFAULT_N}  0.50   0.000   10.000  4.0   -0.9    0.9     0.30  0.100  0.100  0.020  0.020  0
    register_job   "alice_std"           "9_RealisticAlice/alice_std"        ${DEFAULT_N}  0.50   0.000   10.000  4.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "alice_std_neg"       "9_RealisticAlice/alice_std_neg"    ${DEFAULT_N}  -0.50  0.000   10.000  4.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "alice_tight"         "9_RealisticAlice/alice_tight"      ${DEFAULT_N}  0.50   0.000   10.000  4.0   -0.9    0.9     0.30  0.200  0.200  0.100  0.100  0
    register_job   "alice_ring"          "9_RealisticAlice/alice_ring"       ${DEFAULT_N}  0.50   0.500   1.500   4.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 10: MAGNETIC FIELD  [AEE field dependence + HEE with field dependence]
# Fixed standard DCA cuts (0.05 / 0.05 cm). pT cuts = 0.2 GeV/c
# -----------------------------------------------------------------------------
if family_in_scope 10; then
    run_family_header 10 "MAGNETIC FIELD  [AEE: field sign and strength, with HEE tests using min pT cuts]"

    #               NAME               SUBDIR                          N             Bz      pTminL  pTmaxL  rap  etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "fieldCuts_vlow"    "10_BField/fieldCuts_vlow"      ${DEFAULT_N}  0.0001  0.0     10.0    4.0  -0.9    0.9     0.30  0.200  0.200  0.050  0.050  0
    register_job   "fieldCuts_b050"    "10_BField/fieldCuts_b050"      ${DEFAULT_N}  0.50    0.0     10.0    4.0  -0.9    0.9     0.30  0.200  0.200  0.050  0.050  0
    register_job   "fieldCuts_b100"    "10_BField/fieldCuts_b100"      ${DEFAULT_N}  1.00    0.0     10.0    4.0  -0.9    0.9     0.30  0.200  0.200  0.050  0.050  0
    register_job   "fieldCuts_b050neg" "10_BField/fieldCuts_b050neg"   ${DEFAULT_N}  -0.50   0.0     10.0    4.0  -0.9    0.9     0.30  0.200  0.200  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 11: ASYMMETRIC ETA
# Fixed standard DCA cuts (0.05 / 0.05 cm), as low as possible for ALICE algorithms.
# Softest pT cuts as well, 0.15 GeV/c (could even go lower actually, as that can usually be
# the pT cut applied on Lambdas, not on their daughters).
# 1) [-0.85. 0.91] should be equivalent to having a ZVtx shift of 10cm to the positive side,
# emulated into an eta pseudorapidity cut because I don't really have Z vertices here (field
# is constant in XYZ, particularly in Z, and the PV of the Toy Model is set as (0,0,0)).
# -----------------------------------------------------------------------------
if family_in_scope 11; then
    run_family_header 11 "ASYMMETRIC ETA"

    #               NAME                   SUBDIR                              N             Bz     pTminL  pTmaxL  rap  etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "pos10cmZvtx"           "11_EtaAsym/pos10cmZvtx"            ${DEFAULT_N}  0.50   0.0     10.0    4.0  -0.85   0.91    0.30  0.150  0.150  0.050  0.050  0
    register_job   "neg10cmZvtx"           "11_EtaAsym/neg10cmZvtx"            ${DEFAULT_N}  0.50   0.0     10.0    4.0  -0.91   0.85    0.30  0.150  0.150  0.050  0.050  0
    register_job   "pos10cmZvtxNegField"   "11_EtaAsym/pos10cmZvtxNegField"    ${DEFAULT_N}  -0.50  0.0     10.0    4.0  -0.85   0.91    0.30  0.150  0.150  0.050  0.050  0
    register_job   "EtaMin130etaMax090"    "11_EtaAsym/EtaMin130etaMax090"     ${DEFAULT_N}  0.50   0.0     10.0    4.0  -1.10   0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "EtaMin110etaMax090"    "11_EtaAsym/EtaMin110etaMax090"     ${DEFAULT_N}  0.50   0.0     10.0    4.0  -1.10   0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "EtaMin070etaMax090"    "11_EtaAsym/EtaMin070etaMax090"     ${DEFAULT_N}  0.50   0.0     10.0    4.0  -0.70   0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "EtaMin050etaMax090"    "11_EtaAsym/EtaMin050etaMax090"     ${DEFAULT_N}  0.50   0.0     10.0    4.0  -0.50   0.90    0.30  0.150  0.150  0.050  0.050  0
    register_job   "EtaMin030etaMax090"    "11_EtaAsym/EtaMin030etaMax090"     ${DEFAULT_N}  0.50   0.0     10.0    4.0  -0.30   0.90    0.30  0.150  0.150  0.050  0.050  0
fi

# -----------------------------------------------------------------------------
# FAMILY 12: LAMBDA MAX RAPIDITY DEPENDENCE
# The rapidity cut introduced when generating Lambdas is also, on itself, an
# acceptance cut introduced in the generated sample. If we stick with rap = 5.0,
# we are only ignoring 0.009% of Lambdas in cosTheta (for pTLam = 5 GeV/c).
# But from those, only the Lambdas with rapidity |yLam| <= |max(etaMin, etaMax) + 0.67|
# will pass the etaGate, and the rest will not contribute to the statistics! Therefore,
# we must establish a sweet-spot that generates unbiased withoutEtaGate histograms that
# don't have an implicit acceptance cut introduced at the yLambda sampling level, yet
# a sweet-spot that also gives as much statistics as possible to the withEtaGate histograms.
# This is also a way of measuring the strength of the polarization produced by this rapidity
# cut at generator level, which acts as if an eta cut was applied to data.
# - The first test will see if the safeguard that generates the tightest cut possible will work
# - Then we progressively increase the range of the cut
# - The last test should have absolutely no polarization at all in the withoutEtaGate-noCuts case.
#   Anything that shows in the plot should be polarization induced by statistical fluctuation.
# -----------------------------------------------------------------------------
if family_in_scope 12; then
    run_family_header 12 "LAMBDA RAPIDITY"

    #               NAME              SUBDIR                          N             Bz     pTminL  pTmaxL  rap   etaMin  etaMax  T     pTp    pTpi   dcaP   dcaPi  seed
    register_job   "yLam0_5GuardTest" "12_LambdaRap/yLam0_5GuardTest" ${DEFAULT_N}  0.50   0.0     10.0    0.5   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam2"            "12_LambdaRap/yLam2"            ${DEFAULT_N}  0.50   0.0     10.0    2.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam3"            "12_LambdaRap/yLam3"            ${DEFAULT_N}  0.50   0.0     10.0    3.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam4"            "12_LambdaRap/yLam4"            ${DEFAULT_N}  0.50   0.0     10.0    4.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam5"            "12_LambdaRap/yLam5"            ${DEFAULT_N}  0.50   0.0     10.0    5.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam6"            "12_LambdaRap/yLam6"            ${DEFAULT_N}  0.50   0.0     10.0    6.0   -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
    register_job   "yLam10"           "12_LambdaRap/yLam10"           ${DEFAULT_N}  0.50   0.0     10.0    10.0  -0.9    0.9     0.30  0.150  0.150  0.050  0.050  0
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
echo "  Multiply       : ${MULTIPLY}  (shards per job; total procs = jobs x multiply)"
echo "  NUMA nodes     : ${N_NUMA_NODES}  (numactl $([ "${NUMA_AVAILABLE}" -eq 1 ] && echo 'available -- pinning active' || echo 'NOT found -- no pinning'))"
echo "  Peak processes : $(( N_JOBS * MULTIPLY ))  (across ${N_NUMA_NODES} NUMA node(s))"
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
    IFS=':' read -r NAME SUBDIR N BZ PTMIN_LAM PTMAX_LAM RAPMAX ETAMIN ETAMAX T \
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
        --env EXE_BIN \
        --env GEN_MACRO \
        --env PLT_MACRO \
        --env PLOT_ONLY \
        --env MULTIPLY \
        --env NUMA_AVAILABLE \
        --env N_NUMA_NODES \
        --env HADD_EXE \
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

# =============================================================================
# CROSS-FAMILY COMPARISON
# Runs compareHelicityFamilies.cxx after all per-run generator + plotter jobs
# finish. Produces compare_families.root in BASE_DIR with within-family,
# cross-family, and brute-force comparison canvases.
#
# Missing generator ROOT files are silently skipped inside the macro, so
# this is safe to call even after a partial family run (--family N).
# The comparison always runs regardless of DISPATCH_STATUS so that partial
# results (from successfully completed jobs) are still captured.
# --> This also works if you use --plot to avoid regenerating!
# =============================================================================

CMP_MACRO="${SCRIPT_DIR}/compareHelicityFamilies.cxx"
CMP_ROOT="${BASE_DIR}/compare_families.root"
CMP_LOG="${LOG_DIR}/compareHelicityFamilies.log"

if [[ -f "${CMP_MACRO}" ]]; then
    echo ""
    echo "============================================================"
    echo "  Cross-family comparison"
    echo "  Macro   : ${CMP_MACRO}"
    echo "  Output  : ${CMP_ROOT}"
    echo "  Log     : ${CMP_LOG}"
    echo "============================================================"

    {
        echo "=== compareHelicityFamilies started at $(date '+%Y-%m-%d %H:%M:%S') ==="
        echo "Host    : $(hostname)"
        echo "BaseDir : ${BASE_DIR}"
        echo ""
    } > "${CMP_LOG}"

    T_CMP0=$(date +%s)
    if eval "${ROOT_EXE} -l -b -q '${CMP_MACRO}(\"${BASE_DIR}\")'" \
            >> "${CMP_LOG}" 2>&1; then
        T_CMP1=$(date +%s)
        echo "  [DONE] compareHelicityFamilies  ($(( T_CMP1 - T_CMP0 )) s)"
        echo ""
        echo "  Browse results:"
        echo "    root -l '${CMP_ROOT}'"
    else
        CMP_EXIT=$?
        echo "  [WARN] compareHelicityFamilies failed (exit ${CMP_EXIT})."
        echo "         See log: ${CMP_LOG}"
    fi
else
    echo ""
    echo "  [SKIP] compareHelicityFamilies.cxx not found at:"
    echo "         ${CMP_MACRO}"
    echo "         Place it alongside runHelicityToyModel.sh to enable"
    echo "         automatic cross-family comparison."
fi

# =============================================================================
# CLEANUP -- Removing the .exe as it is no longer needed
# =============================================================================
if [[ ${PLOT_ONLY} -eq 0 && -f "${EXE_BIN:-}" ]]; then
    rm -f "${EXE_BIN}"
    echo ""
    echo "  Cleaned up compiled executable."
fi

exit ${DISPATCH_STATUS}