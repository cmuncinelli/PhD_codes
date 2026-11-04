# Sourced by the run scripts. Checks that the right aliBuild environment is loaded.
#
# What is needed:
#   - o2-sim, ROOT with the O2 libraries   -> provided by O2 (also pulled in by O2Physics)
#   - o2dpg_sim_workflow.py & co.          -> provided by O2DPG (a scripts-only package: loading it ALONE gives no O2)
# Either of these works:
#   alienv enter O2sim/latest                                 # official MC meta-package: O2Physics + O2DPG + QC + generators
#   alienv enter O2Physics/latest O2DPG/latest-master-o2      # what you already have built (versions must be compatible)
injpol_require() { # usage: injpol_require o2 [o2dpg]
  local missing=0
  for what in "$@"; do
    case "$what" in
      o2)
        if [ -z "${O2_ROOT:-}" ] || ! command -v o2-sim >/dev/null 2>&1; then
          echo "ERROR: O2 is not loaded (O2_ROOT unset or o2-sim not in PATH)."; missing=1
        fi ;;
      o2dpg)
        if [ -z "${O2DPG_ROOT:-}" ] || [ ! -x "${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py" ]; then
          echo "ERROR: O2DPG is not loaded (O2DPG_ROOT unset)."; missing=1
        fi ;;
    esac
  done
  if [ -z "${INJPOL_ROOT:-}" ] || [ ! -f "${INJPOL_ROOT}/generator/generator_pythia8_injectedPol.C" ]; then
    echo "ERROR: export INJPOL_ROOT=/path/to/PhD_codes/InjectedPolO2MC"; missing=1
  fi
  if [ "$missing" -ne 0 ]; then
    echo "Load an environment with both O2 and O2DPG, e.g.:"
    echo "  alienv enter O2Physics/latest O2DPG/latest-master-o2     (or: alienv enter O2sim/latest)"
    exit 1
  fi
}

# Print only the lines worth reading from a (very long) o2-sim / workflow log.
injpol_summarize_log() { # usage: injpol_summarize_log logfile
  echo "---- [InjectedPol] lines from $1 ----"
  grep -a "\[InjectedPol\]" "$1" | sort | uniq -c | sort -rn | head -40
  echo "---- errors / fatals / ROOT-macro problems in $1 (first 30) ----"
  grep -a -n -i -E "\[fatal\]|\[error\]|error:|segmentation|cannot open|Failed to|R__ADD_INCLUDE_PATH|no matching function" "$1" | head -30
  echo "---- (full log: $1, $(wc -l < "$1") lines) ----"
}
