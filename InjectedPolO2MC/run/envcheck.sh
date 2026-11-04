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
# Known harmless messages (see README, "Known harmless messages") are filtered out of the error list; the TPC-loopers
# failures (an O2 issue: the event is simply regenerated) are counted separately.
injpol_summarize_logs() { # usage: injpol_summarize_logs log1 [log2 ...]
  local f nloop=0 nout=0
  echo "---- [InjectedPol] generator messages ----"
  for f in "$@"; do [ -f "$f" ] && grep -a -h "\[InjectedPol\]" "$f"; done | sed 's/^\[[0-9:]*\]//' | sort | uniq -c | sort -rn | head -20
  echo "---- errors / fatals / exceptions (harmless known messages filtered; first 30) ----"
  for f in "$@"; do
    [ -f "$f" ] || continue
    nloop=$((nloop + $(grep -a -c "Failed to import loopers particles" "$f")))
    nout=$((nout + $(grep -a -c "GeneratorFromO2Kine: Ran out of events" "$f")))
    grep -a -n -i -E "\[fatal\]|\[error\]|error:|error in <|segmentation|exception|cannot open|failed to|no matching function" "$f" \
      | grep -a -v -E "TNetXNGFile::Open.*No servers are available|G4Exception-(START|END)|mat031|Beware of errors and use common sense|Failed to import loopers particles|ReadEvent failed in finalizeEvent|ReadEvent failed for generator ALICEo2" \
      | sed "s|^|$(basename "$f"):|" | head -30
  done
  [ "$nloop" -gt 0 ] && echo "[loopers] ${nloop} TPC-loopers import failure(s), each followed by a regenerated event: harmless with a recent O2 or with the
          workaround in logs/confKey_used.txt; on O2 2026-07-24..09-03 WITHOUT it they make the run hang later (DESIGN.md 4c)"
  for f in "$@"; do
    [ -f "$f" ] && awk '/Flat gas loopers: [0-9]+/{for(i=1;i<=NF;i++) if($i=="loopers:" && $(i+1)+0>5000000){print "[ALARM] runaway TPC loopers in " FILENAME ": " $0; exit}}' "$f"
  done
  [ "$nout" -gt 0 ] && echo "[WARNING] 'GeneratorFromO2Kine: Ran out of events' seen: the generator ran in a SEPARATE step (O2DPG --event-gen-mode separated), which DROPS the Lambda daughters. Use --event-gen-mode integrated."
  for f in "$@"; do [ -f "$f" ] && echo "     (full log: $(readlink -f "$f"), $(wc -l < "$f") lines)"; done
  return 0
}

# Move the o2-sim / O2DPG log files of a directory tree into <dir>/logs/<same relative path>.
# Only call it when a run is FINISHED: the O2DPG runner uses <task>.log_done to skip finished tasks on a re-run.
injpol_move_logs() { # usage: injpol_move_logs <dir>
  local dir="$1" f rel
  ( cd "${dir}" || exit 0
    find . -path ./logs -prune -o -type f \( -name "*.log" -o -name "*.log_done" -o -name "*.log_time" -o -name "*_serverlog" \
         -o -name "*_workerlog*" -o -name "*_mergerlog" \) -print | while read -r f; do
      rel="${f#./}"; mkdir -p "logs/$(dirname "${rel}")"; mv "${f}" "logs/${rel}"
    done )
}

# Disk usage of a directory, total and per event, with a warning for small samples.
injpol_disk_report() { # usage: injpol_disk_report <dir> <nEvents>
  local dir="$1" nev="$2" tot logs
  tot=$(du -sb "${dir}" | cut -f1); logs=$( [ -d "${dir}/logs" ] && du -sb "${dir}/logs" | cut -f1 || echo 0)
  awk -v t="${tot}" -v l="${logs}" -v n="${nev}" 'function h(x){u="B";if(x>1024){x/=1024;u="kB"};if(x>1024){x/=1024;u="MB"};if(x>1024){x/=1024;u="GB"};return sprintf("%.1f %s",x,u)}
    BEGIN{printf "[disk] %s total (logs %s)", h(t), h(l); if(n>0) printf ", %s per event over %d events", h(t/n), n; printf "\n";
          if(n>0 && n<500) printf "[disk] WARNING: %d events only -- per-event numbers are dominated by fixed per-run files (geometry, GRP,\n       CCDB objects, logs) and overestimate large productions.\n", n}'
}

# Run checks/checkKine.C COMPILED with ACLiC (the O2 macro convention). Interpreting it would make ROOT
# auto-parse the VMC dictionary (TMCProcess) recursively; compiling avoids that. ACLiC build products go
# to ${OUT_BASE}/aclic, never into the repository.
injpol_check_kine() { # usage: injpol_check_kine <kine prefix> <logfile> (the log's directory is created if needed)
  mkdir -p "$(dirname "$2")"
  root -l -b -q \
    -e "gSystem->SetBuildDir(\"${OUT_BASE}/aclic\", kTRUE); gSystem->Load(\"libO2CommonUtils\"); gSystem->Load(\"libO2SimulationDataFormat\"); gSystem->Load(\"libO2Steer\");" \
    "${INJPOL_ROOT}/checks/checkKine.C+(\"$1\")" > "$2" 2>&1
  grep -a -E "^\[ OK \]|^\[FAIL\]|^\[info\]|^generator-level|^checkKine:|WARNING|  event " "$2" \
    || { echo "checkKine.C produced no results; first errors from $2:"; grep -a -n -i "error" "$2" | head -15; }
  echo "     (full ROOT output: $(readlink -f "$2"))"
}

# --- TPC-loopers workaround, chosen from the installed O2 version (DESIGN.md section 4c) ---
# O2 versions between the geometrical protection (2026-07-24, 51d8507) and its fix (c575d91, first in daily-20260905-0000) treat "all candidate
# loopers outside the TPC" as a generation FAILURE. Each failure makes the primary server retry the event, which advances the
# loopers' internal event counter once more; at the end of the timeframe it then reads past the collision list
# (TPCLoopers.cxx: mInteractionTimeRecords[mCurrentEvent]) and asks for ~1e9 loopers -> the run hangs for hours.
# INJPOL_LOOPERS = auto (default) | default | noGeomProtection | veto
injpol_loopers_confkey() {
  local mode="${INJPOL_LOOPERS:-auto}" ver date inc="${O2_ROOT}/include/Generators"
  if [ "${mode}" = "auto" ]; then
    ver=$(basename "${O2_ROOT}")
    date=$(echo "${ver}" | grep -oE '20[0-9]{6}' | head -1)
    if [ -n "${date}" ]; then
      if   [ "${date}" -lt 20260724 ]; then mode=default           # no geometrical protection yet -> no failure path
      elif [ "${date}" -lt 20260905 ]; then mode=noGeomProtection  # affected versions (daily-20260904-0700 still lacks the fix)
      else mode=default; fi                                          # fixed: c575d91 is in daily-20260905-0000 and later
      echo "[loopers] O2 version '${ver}' (date ${date}) -> TPC loopers mode: ${mode}" >&2
    elif grep -q "ownsContextCollisions" "${inc}/TPCLoopers.h" 2>/dev/null; then
      mode=default
      echo "[loopers] O2 version '${ver}' not dated; headers contain the 2026-10-31 loopers fixes -> mode: ${mode}" >&2
    elif grep -q "geomProtection" "${inc}/TPCLoopersParam.h" 2>/dev/null; then
      mode=noGeomProtection
      echo "[loopers] O2 version '${ver}' not dated; may be affected by the loopers retry bug -> mode: ${mode}" >&2
    else
      mode=default
      echo "[loopers] O2 version '${ver}' not dated; no geometrical protection in its headers -> mode: ${mode}" >&2
    fi
  fi
  case "${mode}" in
    default)          echo "" ;;
    noGeomProtection) echo "GenTPCLoopers.geomProtection=false" ;;
    veto)             echo "GenTPCLoopers.loopersVeto=true" ;;
    *) echo "[loopers] unknown INJPOL_LOOPERS='${mode}'" >&2; echo "" ;;
  esac
}
