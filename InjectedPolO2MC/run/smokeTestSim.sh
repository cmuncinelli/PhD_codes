#!/bin/bash
# InjectedPolO2MC -- validation step 1: o2-sim only (generator + Geant4), a handful of events, then checkKine.C.
# TEST CODE, see ../README.md.
#
# Usage (needs O2 only, e.g. `alienv enter O2Physics/latest`):
#   export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC
#   ${INJPOL_ROOT}/run/smokeTestSim.sh [modelConfig] [nEvents] [seed]
# Defaults: configs/ring_leadP_0p30.cfg, 10 events, seed 42. Runs in ./smoke_<config name>/.
# The o2-sim output (~20k lines) goes to o2sim.log only; a short summary is printed instead.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/envcheck.sh"
injpol_require o2

MODELCFG=${1:-${INJPOL_ROOT}/generator/configs/ring_leadP_0p30.cfg}
NEV=${2:-10}
SEED=${3:-42}
NWORKERS=${NWORKERS:-4}

WORKDIR=smoke_$(basename "${MODELCFG}" .cfg)
mkdir -p "${WORKDIR}" && cd "${WORKDIR}" || exit 1

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${MODELCFG}")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# --skipModules ZDC: the O2DPG convention; ZDC sits ~112 m from the IP and does not affect the central barrel.
echo "Running o2-sim (${NEV} events) in ${PWD} -> o2sim.log ..."
o2-sim -n "${NEV}" -g external --configFile injectedPol.ini -e TGeant4 --skipModules ZDC \
       --seed "${SEED}" -j "${NWORKERS}" -o o2sim > o2sim.log 2>&1
SIMRC=$?
injpol_summarize_log o2sim.log
if [ ${SIMRC} -ne 0 ] || [ ! -f o2sim_Kine.root ]; then
  echo "o2-sim FAILED (exit code ${SIMRC}). Please send the summary above (and o2sim.log if needed)."
  exit 1
fi

echo "Running checkKine.C -> checkKine.log ..."
root -l -b -q "${INJPOL_ROOT}/checks/checkKine.C(\"o2sim\")" > checkKine.log 2>&1
# The check results are the lines starting with [ OK ], [FAIL], [info] and the summary header
grep -a -E "^\[ OK \]|^\[FAIL\]|^\[info\]|^generator-level|^checkKine:|WARNING|  event " checkKine.log
echo "(full ROOT output: ${PWD}/checkKine.log)"
