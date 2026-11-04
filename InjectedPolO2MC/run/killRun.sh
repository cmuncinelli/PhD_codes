#!/bin/bash
# InjectedPolO2MC -- stop ONE production cleanly, without touching anything else running on the machine
# (e.g. O2Physics analysis jobs): only processes of the current user whose working directory is inside the given
# output folder are selected (the O2DPG runner works in that folder; o2-sim and its workers in its tfN/ subfolders).
# Interactive shells (your terminals, `o2envsim` shells) are never selected, even if their current folder is in there.
#
# Usage: killRun.sh <output folder>            list the processes that would be stopped
#        killRun.sh <output folder> --kill     SIGTERM them, then SIGKILL whatever is left after 15 s
set -uo pipefail
DIR=$(readlink -f "${1:?usage: killRun.sh <output folder> [--kill]}")
[ -d "${DIR}" ] || { echo "ERROR: no directory ${DIR}"; exit 1; }

is_interactive_shell() { # a shell started without a script (e.g. your terminal, or an `o2envsim` shell) -> never killed
  local -a argv
  mapfile -d '' -t argv < "/proc/$1/cmdline" 2>/dev/null || return 1
  case "$(basename "${argv[0]#-}")" in bash|sh|zsh|tcsh|csh|ksh|fish) ;; *) return 1 ;; esac
  [ "${#argv[@]}" -le 1 ] && return 0
  local a; for a in "${argv[@]:1}"; do case "${a}" in -*) ;; *) [ -f "${a}" ] && return 1 ;; esac; done
  return 0
}

find_pids() {
  local p cwd
  for p in /proc/[0-9]*; do
    [ -O "${p}" ] || continue                     # only our own processes
    cwd=$(readlink "${p}/cwd" 2>/dev/null) || continue
    case "${cwd}/" in "${DIR}/"*) is_interactive_shell "${p#/proc/}" || echo "${p#/proc/}" ;; esac
  done
}

PIDS=$(find_pids | grep -v -x "$$")
if [ -z "${PIDS}" ]; then echo "No process of yours runs inside ${DIR}."; exit 0; fi
echo "Processes running inside ${DIR}:"
ps -o pid,etime,pcpu,rss,args -p $(echo ${PIDS} | tr ' ' ',') | cut -c1-200
if [ "${2:-}" != "--kill" ]; then echo "(dry run: add --kill to stop them)"; exit 0; fi

kill -TERM ${PIDS} 2>/dev/null
for i in $(seq 1 15); do sleep 1; LEFT=$(find_pids | grep -v -x "$$"); [ -z "${LEFT}" ] && break; done
LEFT=$(find_pids | grep -v -x "$$")
if [ -n "${LEFT}" ]; then echo "Still alive after 15 s, sending SIGKILL: ${LEFT}"; kill -KILL ${LEFT} 2>/dev/null; sleep 2; fi
LEFT=$(find_pids | grep -v -x "$$")
[ -z "${LEFT}" ] && echo "All processes inside ${DIR} stopped." || { echo "Could not stop: ${LEFT}"; exit 1; }
