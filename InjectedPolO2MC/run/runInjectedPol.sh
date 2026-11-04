#!/bin/bash
# InjectedPolO2MC -- full chain (generator -> Geant4 -> digitization -> reconstruction -> AO2D) with O2DPG.
# TEST CODE, see ../README.md. Unanchored for now (see DESIGN.md, "Anchoring").
#
# Usage (needs O2 AND O2DPG, e.g. `o2envsim`):
#   MODELCFG=<model config> NSIGEVENTS=50 NTIMEFRAMES=1 runInjectedPol.sh
# Output: ${OUT_BASE}/full/<config name>_seed<SEED>/tf*/AO2D.root (+ o2sim kinematics per TF, kept for truth checks).
# Logs:   .../logs/ (our logs immediately; the O2DPG/o2-sim task logs are moved there once the run has SUCCEEDED,
#         so that a failed run can still be inspected and resumed in place).
# Progress from another terminal:  ${INJPOL_ROOT}/run/progress.py <output folder> --watch 30
#
# Every flag that departs from the o2dpg_sim_workflow.py defaults is listed and justified below.
set -uo pipefail
REPO_DIR="/home/users/cicerodm/PhD_codes"            # no trailing "/" needed (it is stripped anyway)
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # all outputs go here, never into the repository
source "${REPO_DIR%/}/InjectedPolO2MC/run/envcheck.sh"
injpol_require o2 o2dpg

MODELCFG=${MODELCFG:-${INJPOL_ROOT}/generator/configs/none.cfg}
NWORKERS=${NWORKERS:-8}
NSIGEVENTS=${NSIGEVENTS:-50}       # signal events per timeframe (O2DPG adds "orbits-early" pile-up collisions on top)
NTIMEFRAMES=${NTIMEFRAMES:-1}
INTRATE=${INTRATE:-500000}         # Hz, used by digitization (pile-up); anchored productions take it from the run
SEED=${SEED:-1}                    # O2DPG derives one seed per TF from it; Pythia and the afterburner follow
SIMENGINE=${SIMENGINE:-TGeant4}
CPULIMIT=${CPULIMIT:-${NWORKERS}}
WORKDIR=${WORKDIR:-${OUT_BASE}/full/$(basename "${MODELCFG}" .cfg)_seed${SEED}}

mkdir -p "${WORKDIR}/logs" && cd "${WORKDIR}" || exit 1
cp "${MODELCFG}" model.cfg   # bookkeeping: keep the exact configuration next to the output

cat > injectedPol.ini <<INI
[GeneratorExternal]
fileName = ${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C
funcName = generator_injectedPol("${WORKDIR}/model.cfg")

[GeneratorPythia8]
config = ${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg
INI

# Flags:
#   -gen external -ini ...          : our generator (Pythia8 + afterburner)
#   --event-gen-mode integrated     : REQUIRED. The O2DPG default ("separated") generates events in a separate o2-sim and
#                                     re-reads them with GeneratorFromO2Kine, which keeps only PRIMARY tracks and thus
#                                     drops our Lambda daughters (DESIGN.md section 4b). "integrated" runs the generator
#                                     inside the transport, as in the smoke test.
#   --skipModules ZDC               : O2DPG convention (most O2DPG run scripts); ZDC is ~112 m away, no effect on the barrel
#   -confKey Diamond.width[2]=6.    : 6 cm luminous-region sigma along z, as in the O2DPG run scripts
#   GenTPCLoopers.*                 : workaround for an O2 TPC-loopers bug in some O2 versions, chosen automatically
#                                     (INJPOL_LOOPERS=auto|default|noGeomProtection|veto; DESIGN.md section 4c)
#   (NOT used: --readoutDets restrictions; anchored productions set them to the detectors read out in the anchor run)
LOOPERSKEY=$(injpol_loopers_confkey)
CONFKEY="Diamond.width[2]=6.${LOOPERSKEY:+;${LOOPERSKEY}}"
echo "${CONFKEY}" > logs/confKey_used.txt

echo "Creating workflow in ${WORKDIR} -> logs/workflow_creation.log ..."
${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -eCM 13600 -col pp -gen external -ini "${WORKDIR}/injectedPol.ini" \
    --event-gen-mode integrated \
    -j "${NWORKERS}" -ns "${NSIGEVENTS}" -tf "${NTIMEFRAMES}" -interactionRate "${INTRATE}" \
    -seed "${SEED}" -e "${SIMENGINE}" --skipModules ZDC \
    -confKey "${CONFKEY}" > logs/workflow_creation.log 2>&1 || { tail -30 logs/workflow_creation.log; echo "workflow creation FAILED"; exit 1; }

echo "Running workflow up to AOD -> logs/workflow_runner.log (this takes a while) ..."
echo "  progress from another terminal:  ${INJPOL_ROOT}/run/progress.py ${WORKDIR} --watch 30"
${O2DPG_ROOT}/MC/bin/o2_dpg_workflow_runner.py -f workflow.json -tt aod --cpu-limit "${CPULIMIT}" > logs/workflow_runner.log 2>&1
RUNRC=$?
grep -a "Pipeline done" logs/workflow_runner.log || tail -15 logs/workflow_runner.log
if [ ${RUNRC} -ne 0 ]; then
  echo "workflow runner FAILED (exit ${RUNRC}). The failing task and its log file are named in logs/workflow_runner.log."
  echo "Task logs were left in place (tf1/*.log, tf1/sgn_serverlog, ...) so the run can be inspected or resumed."
  python3 "${INJPOL_ROOT}/run/progress.py" "${WORKDIR}"
  exit 1
fi

for d in tf*/; do
  d=${d%/}
  echo "== ${d}: o2-sim logs + checkKine.C"
  ( cd "${d}" && injpol_summarize_logs sgn_serverlog sgn_workerlog0 sgnsim_${d#tf}.log && injpol_check_kine sgn "${WORKDIR}/logs/${d}/checkKine.log" )
done
for f in tf*/AO2D.root; do
  d=$(dirname "${f}")
  echo "== ${f}: checkAODFlags.C"
  mkdir -p "logs/${d}"
  root -l -b -q "${INJPOL_ROOT}/checks/checkAODFlags.C(\"${f}\")" > "logs/${d}/checkAODFlags.log" 2>&1
  grep -a -E "^\[ OK \]|^\[FAIL\]|^\[WARN\]|^\[info\]|^\[HINT\]|generator-level|  DF_" "logs/${d}/checkAODFlags.log"
done

injpol_move_logs "${WORKDIR}"   # the run succeeded: gather all task logs under logs/
echo "== timing and disk report (also in logs/run_report.txt)"
python3 "${INJPOL_ROOT}/run/progress.py" "${WORKDIR}" --final | sed -n '/wall time per task kind/,$p'
