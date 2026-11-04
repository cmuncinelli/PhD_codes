# Sourced by the InjectedPolO2MC run scripts (do not execute directly).
# Expects REPO_DIR and OUT_BASE to be set by the calling script (hard-coded at its top).
#
# Environment needed:
#   - o2-sim, ROOT with the O2 libraries   -> O2 (also pulled in by O2Physics)
#   - o2dpg_sim_workflow.py & co.          -> O2DPG (scripts-only package: loading it ALONE gives no O2)
# Either of these works:
#   alienv enter O2Physics/latest O2DPG/latest-master-o2
#   alienv enter O2sim/latest              # official MC meta-package: O2Physics + O2DPG + QC + generators

# --- normalise the hard-coded paths (a trailing "/" is stripped, so "/a/b/" and "/a/b" are equivalent) ---
REPO_DIR="${REPO_DIR%/}"
OUT_BASE="${OUT_BASE%/}"
export INJPOL_ROOT="${REPO_DIR}/InjectedPolO2MC" # used by the generator macro (R__ADD_INCLUDE_PATH); never set it by hand

injpol_require() { # usage: injpol_require [o2] [o2dpg] [pythia] [gxx]
  local missing=0
  if [ ! -f "${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C" ]; then
    echo "ERROR: REPO_DIR='${REPO_DIR}' does not contain InjectedPolO2MC/ (edit REPO_DIR at the top of the script)"; missing=1
  fi
  for what in "$@"; do
    case "$what" in
      o2)     { [ -n "${O2_ROOT:-}" ] && command -v o2-sim >/dev/null 2>&1; } || { echo "ERROR: O2 is not loaded (O2_ROOT unset or o2-sim not in PATH)."; missing=1; } ;;
      o2dpg)  { [ -n "${O2DPG_ROOT:-}" ] && [ -x "${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py" ]; } || { echo "ERROR: O2DPG is not loaded (O2DPG_ROOT unset)."; missing=1; } ;;
      pythia) command -v pythia8-config >/dev/null 2>&1 || { echo "ERROR: pythia8-config not found (load e.g. O2Physics/latest)."; missing=1; } ;;
      gxx)    command -v g++ >/dev/null 2>&1 || { echo "ERROR: g++ not found."; missing=1; } ;;
    esac
  done
  if [ "$missing" -ne 0 ]; then
    echo "Load an environment with O2 and O2DPG, e.g.:"
    echo "  alienv enter O2Physics/latest O2DPG/latest-master-o2     (or: alienv enter O2sim/latest)"
    exit 1
  fi
  mkdir -p "${OUT_BASE}" || { echo "ERROR: cannot create OUT_BASE='${OUT_BASE}'"; exit 1; }
}

# Print only the lines worth reading from very long o2-sim / workflow logs.
# Known harmless messages (see README, "Known harmless messages") are filtered out of the error list.
injpol_summarize_logs() { # usage: injpol_summarize_logs log1 [log2 ...]
  local f
  echo "---- [InjectedPol] generator messages ----"
  for f in "$@"; do [ -f "$f" ] && grep -a -h "\[InjectedPol\]" "$f"; done | sed 's/^\[[0-9:]*\]//' | sort | uniq -c | sort -rn | head -20
  echo "---- errors / fatals / exceptions (harmless known messages filtered; first 30) ----"
  for f in "$@"; do
    [ -f "$f" ] || continue
    grep -a -n -i -E "\[fatal\]|\[error\]|error:|error in <|segmentation|exception|cannot open|failed to|no matching function" "$f" \
      | grep -a -v -E "TNetXNGFile::Open.*No servers are available|G4Exception-(START|END)|mat031|Beware of errors and use common sense" \
      | sed "s|^|$(basename "$f"):|" | head -30
  done
  for f in "$@"; do [ -f "$f" ] && echo "     (full log: ${PWD}/$f, $(wc -l < "$f") lines)"; done
}

# Run checks/checkKine.C COMPILED with ACLiC (the O2 macro convention). Interpreting it would make ROOT
# auto-parse the VMC dictionary (TMCProcess) recursively; compiling avoids that. ACLiC build products go
# to ${OUT_BASE}/aclic, never into the repository.
injpol_check_kine() { # usage: injpol_check_kine <kine prefix> <logfile>
  root -l -b -q \
    -e "gSystem->SetBuildDir(\"${OUT_BASE}/aclic\", kTRUE); gSystem->Load(\"libO2CommonUtils\"); gSystem->Load(\"libO2SimulationDataFormat\"); gSystem->Load(\"libO2Steer\");" \
    "${INJPOL_ROOT}/checks/checkKine.C+(\"$1\")" > "$2" 2>&1
  grep -a -E "^\[ OK \]|^\[FAIL\]|^\[info\]|^generator-level|^checkKine:|WARNING|  event " "$2" \
    || { echo "checkKine.C produced no results; first errors from $2:"; grep -a -n -i "error" "$2" | head -15; }
  echo "     (full ROOT output: ${PWD}/$2)"
}
