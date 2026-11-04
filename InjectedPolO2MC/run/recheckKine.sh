#!/bin/bash
# InjectedPolO2MC -- re-run ONLY checks/checkKine.C on existing o2-sim output (no new simulation).
# Usage: recheckKine.sh [directory] [kine prefix]
#   defaults: ${OUT_BASE}/smoke/ring_leadP_0p30, prefix "o2sim" (use "sgn" inside a tfN/ directory of a full chain)
set -uo pipefail
REPO_DIR="/home/users/cicerodm/PhD_codes"            # no trailing "/" needed (it is stripped anyway)
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # all outputs go here, never into the repository
source "${REPO_DIR%/}/InjectedPolO2MC/run/envcheck.sh"
injpol_require o2

DIR=${1:-${OUT_BASE}/smoke/ring_leadP_0p30}
PREFIX=${2:-o2sim}
cd "${DIR}" || { echo "ERROR: no directory ${DIR}"; exit 1; }
[ -f "${PREFIX}_Kine.root" ] || { echo "ERROR: ${DIR}/${PREFIX}_Kine.root not found"; exit 1; }
echo "Running checkKine.C (compiled) in ${DIR} ..."
injpol_check_kine "${PREFIX}" checkKine.log
