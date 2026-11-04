#!/bin/bash
# InjectedPolO2MC -- full chain (generator -> Geant4 -> digitization -> reconstruction -> AO2D) with O2DPG.
# TEST CODE, see ../README.md. Unanchored for now (see DESIGN.md, "Anchoring").
#
# Usage (inside the O2sim / O2PDPSuite environment, O2DPG loaded):
#   export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC
#   MODELCFG=${INJPOL_ROOT}/generator/configs/ring_leadP_0p30.cfg NSIGEVENTS=50 NTIMEFRAMES=1 \
#     ${INJPOL_ROOT}/run/runInjectedPol.sh
# Output: <WORKDIR>/tf*/AO2D.root (+ o2sim kinematics per TF, kept for truth-level checks).
#
# Every flag that departs from the o2dpg_sim_workflow.py defaults is listed and justified below.
set -euo pipefail

[ -z "${O2DPG_ROOT:-}" ] && echo "ERROR: load O2DPG first (O2DPG_ROOT not set)" && exit 1
[ -z "${O2_ROOT:-}" ] && echo "ERROR: load O2 first (O2_ROOT not set)" && exit 1
[ -z "${INJPOL_ROOT:-}" ] && echo "ERROR: export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC" && exit 1

MODELCFG=${MODELCFG:-${INJPOL_ROOT}/generator/configs/none.cfg}
NWORKERS=${NWORKERS:-8}
NSIGEVENTS=${NSIGEVENTS:-50}       # signal events per timeframe
NTIMEFRAMES=${NTIMEFRAMES:-1}
INTRATE=${INTRATE:-500000}         # Hz, used by digitization (pile-up); anchored productions take it from the run
SEED=${SEED:-1}                    # O2DPG derives one seed per TF from it; Pythia and the afterburner follow
SIMENGINE=${SIMENGINE:-TGeant4}
CPULIMIT=${CPULIMIT:-${NWORKERS}}
WORKDIR=${WORKDIR:-injpol_$(basename "${MODELCFG}" .cfg)_seed${SEED}}

mkdir -p "${WORKDIR}" && cd "${WORKDIR}"
cp "${MODELCFG}" model.cfg   # bookkeeping: keep the exact configuration next to the output

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${PWD}/model.cfg")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# Flags:
#   -gen external -ini ...        : our generator (Pythia8 + afterburner)
#   --skipModules ZDC             : O2DPG convention (most O2DPG run scripts); ZDC is ~112 m away, no effect on the barrel
#   -confKey Diamond.width[2]=6.  : 6 cm luminous-region sigma along z, as in the O2DPG run scripts
#   (NOT used: --readoutDets restrictions; anchored productions set them to the detectors read out in the anchor run)
${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -eCM 13600 -col pp -gen external -ini "${PWD}/injectedPol.ini" \
    -j "${NWORKERS}" -ns "${NSIGEVENTS}" -tf "${NTIMEFRAMES}" -interactionRate "${INTRATE}" \
    -seed "${SEED}" -e "${SIMENGINE}" --skipModules ZDC \
    -confKey "Diamond.width[2]=6."

${O2DPG_ROOT}/MC/bin/o2_dpg_workflow_runner.py -f workflow.json -tt aod --cpu-limit "${CPULIMIT}"

for d in tf*/; do
  ( cd "${d}" && root -l -b -q "${INJPOL_ROOT}/checks/checkKine.C(\"sgn\")" 2>&1 | tee checkKine.log )
done
for f in tf*/AO2D.root; do
  root -l -b -q "${INJPOL_ROOT}/checks/checkAODFlags.C(\"${f}\")" 2>&1 | tee "${f%/AO2D.root}/checkAODFlags.log"
done
