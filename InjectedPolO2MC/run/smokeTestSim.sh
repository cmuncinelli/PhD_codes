#!/bin/bash
# InjectedPolO2MC -- validation step 1: o2-sim only (generator + Geant4), a handful of events, then checkKine.C.
# TEST CODE, see ../README.md.
#
# Usage (inside the O2sim / O2PDPSuite environment):
#   export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC
#   ${INJPOL_ROOT}/run/smokeTestSim.sh [modelConfig] [nEvents] [seed]
# Defaults: configs/ring_leadP_0p30.cfg, 10 events, seed 42. Runs in ./smoke_<config name>/.
set -euo pipefail

[ -z "${O2_ROOT:-}" ] && echo "ERROR: load the O2 environment first (O2_ROOT not set)" && exit 1
[ -z "${INJPOL_ROOT:-}" ] && echo "ERROR: export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC" && exit 1

MODELCFG=${1:-${INJPOL_ROOT}/generator/configs/ring_leadP_0p30.cfg}
NEV=${2:-10}
SEED=${3:-42}
NWORKERS=${NWORKERS:-4}

WORKDIR=smoke_$(basename "${MODELCFG}" .cfg)
mkdir -p "${WORKDIR}" && cd "${WORKDIR}"

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${MODELCFG}")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# --skipModules ZDC: the O2DPG convention (used by most O2DPG run scripts); ZDC sits ~112 m from the IP and does
# not affect the central barrel. No other module is skipped (material budget unchanged).
o2-sim -n "${NEV}" -g external --configFile injectedPol.ini -e TGeant4 --skipModules ZDC \
       --seed "${SEED}" -j "${NWORKERS}" -o o2sim 2>&1 | tee o2sim.log

root -l -b -q "${INJPOL_ROOT}/checks/checkKine.C(\"o2sim\")" 2>&1 | tee checkKine.log
