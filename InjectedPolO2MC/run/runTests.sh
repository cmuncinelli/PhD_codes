#!/bin/bash
# InjectedPolO2MC -- validation steps "unit" and 0:
#   - kernel unit tests (g++ only),
#   - stand-alone Pythia closure (needs pythia8-config, e.g. `alienv enter O2Physics/latest`).
# TEST CODE, see ../README.md. Binaries are built as *.exe in ${OUT_BASE}/bin, logs go to ${OUT_BASE}/tests.
#
# Usage: runTests.sh [modelConfig] [nAcceptedEvents] [seed]
#   defaults: configs/ring_leadP_0p30.cfg, 2000 events, seed 12345
set -uo pipefail
REPO_DIR="/home/users/cicerodm/PhD_codes"            # no trailing "/" needed (it is stripped anyway)
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # all outputs go here, never into the repository
source "${REPO_DIR%/}/InjectedPolO2MC/run/envcheck.sh"
injpol_require gxx

MODELCFG=${1:-${INJPOL_ROOT}/generator/configs/ring_leadP_0p30.cfg}
NEV=${2:-2000}
SEED=${3:-12345}
BIN="${OUT_BASE}/bin"; LOGS="${OUT_BASE}/tests"
mkdir -p "${BIN}" "${LOGS}"

echo "== kernel unit tests -> ${LOGS}/testKernel.log"
g++ -O2 -std=c++17 -I"${INJPOL_ROOT}/generator" "${INJPOL_ROOT}/tests/testKernel.cxx" -o "${BIN}/testKernel.exe" \
  && "${BIN}/testKernel.exe" > "${LOGS}/testKernel.log" 2>&1
grep -a -E "^\[FAIL\]|ALL PASSED|FAILED" "${LOGS}/testKernel.log"

if command -v pythia8-config >/dev/null 2>&1; then
  NAME=$(basename "${MODELCFG}" .cfg)
  echo "== stand-alone Pythia closure (${NAME}, ${NEV} events) -> ${LOGS}/standalonePythiaClosure_${NAME}.log"
  # -std=c++17 LAST: pythia8-config may add an older -std, and the last one wins
  g++ -O2 -I"${INJPOL_ROOT}/generator" "${INJPOL_ROOT}/tests/standalonePythiaClosure.cxx" -o "${BIN}/standalonePythiaClosure.exe" \
      $(pythia8-config --cxxflags --libs) -std=c++17 \
    && "${BIN}/standalonePythiaClosure.exe" "${INJPOL_ROOT}/generator/pythia8_pp_13600.cfg" "${MODELCFG}" "${NEV}" "${SEED}" \
         > "${LOGS}/standalonePythiaClosure_${NAME}.log" 2>&1
  sed -n '/=== summary ===/,$p' "${LOGS}/standalonePythiaClosure_${NAME}.log"
else
  echo "== pythia8-config not found: skipping the stand-alone Pythia closure (load e.g. O2Physics/latest)"
fi
