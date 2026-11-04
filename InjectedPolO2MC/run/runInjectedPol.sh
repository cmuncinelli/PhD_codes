#!/bin/bash
# InjectedPolO2MC -- full chain (generator -> Geant4 -> digitization -> reconstruction -> AO2D) with O2DPG.
# TEST CODE, see ../README.md. Unanchored for now (see DESIGN.md, "Anchoring").
#
# Usage (needs O2 AND O2DPG, e.g. `alienv enter O2Physics/latest O2DPG/latest-master-o2`):
#   MODELCFG=<model config> NSIGEVENTS=50 NTIMEFRAMES=1 runInjectedPol.sh
# Output: ${OUT_BASE}/full/<config name>_seed<SEED>/tf*/AO2D.root (+ o2sim kinematics per TF, kept for truth checks).
# Workflow output goes to workflow_creation.log / workflow_runner.log (+ per-task logs written by the runner).
#
# Every flag that departs from the o2dpg_sim_workflow.py defaults is listed and justified below.
set -uo pipefail
REPO_DIR="/home/users/cicerodm/PhD_codes"            # no trailing "/" needed (it is stripped anyway)
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # all outputs go here, never into the repository
source "${REPO_DIR%/}/InjectedPolO2MC/run/envcheck.sh"
injpol_require o2 o2dpg

MODELCFG=${MODELCFG:-${INJPOL_ROOT}/generator/configs/none.cfg}
NWORKERS=${NWORKERS:-8}
NSIGEVENTS=${NSIGEVENTS:-50}       # signal events per timeframe
NTIMEFRAMES=${NTIMEFRAMES:-1}
INTRATE=${INTRATE:-500000}         # Hz, used by digitization (pile-up); anchored productions take it from the run
SEED=${SEED:-1}                    # O2DPG derives one seed per TF from it; Pythia and the afterburner follow
SIMENGINE=${SIMENGINE:-TGeant4}
CPULIMIT=${CPULIMIT:-${NWORKERS}}
WORKDIR=${WORKDIR:-${OUT_BASE}/full/$(basename "${MODELCFG}" .cfg)_seed${SEED}}

mkdir -p "${WORKDIR}" && cd "${WORKDIR}" || exit 1
cp "${MODELCFG}" model.cfg   # bookkeeping: keep the exact configuration next to the output

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${WORKDIR}/model.cfg")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# Flags:
#   -gen external -ini ...        : our generator (Pythia8 + afterburner)
#   --skipModules ZDC             : O2DPG convention (most O2DPG run scripts); ZDC is ~112 m away, no effect on the barrel
#   -confKey Diamond.width[2]=6.  : 6 cm luminous-region sigma along z, as in the O2DPG run scripts
#   (NOT used: --readoutDets restrictions; anchored productions set them to the detectors read out in the anchor run)
echo "Creating workflow in ${WORKDIR} -> workflow_creation.log ..."
${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -eCM 13600 -col pp -gen external -ini "${WORKDIR}/injectedPol.ini" \
    -j "${NWORKERS}" -ns "${NSIGEVENTS}" -tf "${NTIMEFRAMES}" -interactionRate "${INTRATE}" \
    -seed "${SEED}" -e "${SIMENGINE}" --skipModules ZDC \
    -confKey "Diamond.width[2]=6." > workflow_creation.log 2>&1 || { tail -30 workflow_creation.log; echo "workflow creation FAILED"; exit 1; }

echo "Running workflow up to AOD -> workflow_runner.log (this takes a while) ..."
${O2DPG_ROOT}/MC/bin/o2_dpg_workflow_runner.py -f workflow.json -tt aod --cpu-limit "${CPULIMIT}" > workflow_runner.log 2>&1
RUNRC=$?
tail -15 workflow_runner.log
if [ ${RUNRC} -ne 0 ]; then
  echo "workflow runner FAILED (exit ${RUNRC}). The failing task and its log file are named above."
  echo "For a generator/transport problem, the o2-sim logs are in tf1/ (sgn_serverlog, sgn_workerlog0, ...)."
  exit 1
fi

for d in tf*/; do
  echo "== ${d}: o2-sim logs + checkKine.C"
  ( cd "${d}" && injpol_summarize_logs sgn_serverlog sgn_workerlog0 && injpol_check_kine sgn checkKine.log )
done
for f in tf*/AO2D.root; do
  echo "== ${f}: checkAODFlags.C"
  root -l -b -q "${INJPOL_ROOT}/checks/checkAODFlags.C(\"${f}\")" > "${f%/AO2D.root}/checkAODFlags.log" 2>&1
  grep -a -E "^\[ OK \]|^\[FAIL\]|^\[WARN\]|^\[info\]|generator-level|  DF_" "${f%/AO2D.root}/checkAODFlags.log"
done
