#!/bin/bash
# InjectedPolO2MC -- validation step 1 ("smoke test", see README naming): o2-sim only (generator + Geant4) on a
# handful of events, then checkKine.C. TEST CODE, see ../README.md.
#
# Usage (needs O2, e.g. `o2envsim`):  smokeTestSim.sh [modelConfig] [nEvents] [seed]
# Defaults: configs/ring_leadP_0p30.cfg, 10 events, seed 42.
# Output: ${OUT_BASE}/tests/smoke/<config name>/  (ROOT files) and .../logs/ (all logs). Only a summary is printed.
set -uo pipefail
REPO_DIR="/home/users/cicerodm/PhD_codes"            # no trailing "/" needed (it is stripped anyway)
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # all outputs go here, never into the repository
source "${REPO_DIR%/}/InjectedPolO2MC/run/envcheck.sh"
injpol_require o2

MODELCFG=${1:-${INJPOL_ROOT}/generator/configs/ring_leadP_0p30.cfg}
NEV=${2:-10}
SEED=${3:-42}
NWORKERS=${NWORKERS:-4}

WORKDIR="${OUT_BASE}/tests/smoke/$(basename "${MODELCFG}" .cfg)"
mkdir -p "${WORKDIR}/logs" && cd "${WORKDIR}" || exit 1
cp "${MODELCFG}" model.cfg

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${WORKDIR}/model.cfg")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# --skipModules ZDC: the O2DPG convention; ZDC sits ~112 m from the IP and does not affect the central barrel.
echo "Running o2-sim (${NEV} events, ${NWORKERS} workers) in ${WORKDIR} ..."
T0=$(date +%s)
o2-sim -n "${NEV}" -g external --configFile injectedPol.ini -e TGeant4 --skipModules ZDC \
       --seed "${SEED}" -j "${NWORKERS}" -o o2sim > logs/o2sim.log 2>&1
SIMRC=$?
T1=$(date +%s)
injpol_move_logs "${WORKDIR}"   # o2sim_serverlog (generator), o2sim_workerlog* (Geant4), o2sim_mergerlog -> logs/
injpol_summarize_logs logs/o2sim.log logs/o2sim_serverlog logs/o2sim_workerlog0
if [ ${SIMRC} -ne 0 ] || [ ! -f o2sim_Kine.root ]; then
  echo "o2-sim FAILED (exit code ${SIMRC}). Please send the summary above (and the logs if needed)."
  exit 1
fi
NDONE=$(grep -a -c "EVENT FINISHED" logs/o2sim.log)
echo "[time] o2-sim: $((T1 - T0)) s wall for ${NDONE} events with ${NWORKERS} workers -> $(awk -v t=$((T1 - T0)) -v n="${NDONE}" 'BEGIN{printf "%.1f", (n>0)?t/n:0}') s/event"
echo "       (includes ~1 min of fixed initialisation: geometry, CCDB, Geant4 tables; meaningless per event for 10 events)"
injpol_disk_report "${WORKDIR}" "${NDONE}"

echo "Running checkKine.C (compiled) ..."
injpol_check_kine o2sim logs/checkKine.log
