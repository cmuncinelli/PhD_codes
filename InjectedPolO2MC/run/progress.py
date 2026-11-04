#!/usr/bin/env python3
"""InjectedPolO2MC -- progress, ETA, timing and disk-usage report for an O2DPG workflow directory.

TEST CODE, see ../README.md. Read-only: it only inspects files written by o2_dpg_workflow_runner.py and o2-sim.

Usage:
  progress.py <workdir>                 one-shot report (e.g. from another terminal while the chain runs)
  progress.py <workdir> --watch 30      refresh every 30 s until the workflow finishes (Ctrl-C to stop)
  progress.py <workdir> --final         full report, also written to <workdir>/logs/run_report.txt

What it reads (all written by O2DPG / O2 anyway):
  workflow.json                         task graph (names, dependencies, working directories)
  <tfN>/<task>.log_done                 marker left by the O2 task wrapper when a task succeeded
  <tfN>/<task>.log_time                 "#walltime <s>" written by the task wrapper (exact wall time per task)
  pipeline_action_<pid>.log             runner log; "Submitting task <name>" lines give task start times
  <tfN>/sgnsim_<N>.log                  o2-sim steering output; "EVENT FINISHED" lines give event-level progress and
                                        "Setting number of events to simulate to N" the expected total
  <tfN>/sgn_serverlog                   generator output; used to detect runaway TPC loopers (DESIGN.md section 4c)
Warnings: a running transport with no log activity for --stall-minutes (default 20) is flagged as possibly stalled.
After a successful run, runInjectedPol.sh moves these logs to <workdir>/logs/; both places are searched.
"""
import argparse
import datetime as dt
import glob
import json
import os
import re
import sys
import time

TARGET_DEFAULT = "aod"  # same target as `o2_dpg_workflow_runner.py -tt aod`
STALL_MINUTES = 20      # warn if a running transport task has written no log line for this long (--stall-minutes)


def find_file(workdir, rel):
    """Return the path of a task-related file, before (workdir/rel) or after (workdir/logs/rel) the log move."""
    for base in (workdir, os.path.join(workdir, "logs")):
        p = os.path.join(base, rel)
        if os.path.exists(p):
            return p
    return None


def load_tasks(workdir, targets):
    with open(os.path.join(workdir, "workflow.json")) as f:
        stages = json.load(f)["stages"]
    byname = {s["name"]: s for s in stages}
    # same selection as the runner: tasks whose name re.match-es a target, plus all their dependencies
    selected, todo = set(), [s["name"] for s in stages if any(re.match(t, s["name"]) for t in targets)]
    while todo:
        n = todo.pop()
        if n in selected or n not in byname:
            continue
        selected.add(n)
        todo.extend(byname[n].get("needs", []))
    order = [s["name"] for s in stages if s["name"] in selected]
    return [byname[n] for n in order]


def start_times(workdir):
    """task name -> datetime of the last 'Submitting task <name>' in the newest runner action log."""
    logs = sorted(glob.glob(os.path.join(workdir, "pipeline_action_*.log")) + glob.glob(os.path.join(workdir, "logs", "pipeline_action_*.log")),
                  key=os.path.getmtime)
    starts = {}
    if not logs:
        return starts, None
    first = None
    rx = re.compile(r"^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d),\d+ \w+ Submitting task (\S+)")
    with open(logs[-1], errors="replace") as f:
        for line in f:
            m = rx.match(line)
            if m:
                t = dt.datetime.strptime(m.group(1), "%Y-%m-%d %H:%M:%S")
                starts[m.group(2)] = t
                first = first or t
    return starts, first


def task_status(workdir, task):
    rel = os.path.join(task.get("cwd", "") or "", task["name"] + ".log")
    done = find_file(workdir, rel + "_done") is not None
    log = find_file(workdir, rel)
    wall = None
    tf = find_file(workdir, rel + "_time")
    if tf:
        with open(tf, errors="replace") as f:
            m = re.search(r"#walltime\s+([\d.]+)", f.read())
            if m:
                wall = float(m.group(1))
    state = "done" if done else ("running/failed" if log else "pending")
    return state, wall


def base_name(name):
    return re.sub(r"_\d+$", "", name)


def sim_events(workdir, task):
    """(finished, distributed, expected) events from an o2-sim steering log (sgnsim_N.log).
    expected comes from o2-sim's 'Setting number of events to simulate to N' (collision context), None if not found."""
    rel = os.path.join(task.get("cwd", "") or "", task["name"] + ".log")
    p = find_file(workdir, rel)
    if not p:
        return 0, 0, None
    fin, dist, expected = 0, 0, None
    rx = re.compile(r"Setting number of events to simulate to (\d+)")
    with open(p, errors="replace") as f:
        for line in f:
            if "EVENT FINISHED" in line:
                fin += 1
            elif "DISTRIBUTING EVENT" in line:
                dist += 1
            elif expected is None:
                m = rx.search(line)
                if m:
                    expected = int(m.group(1))
    return fin, dist, expected


def sim_health(workdir, task, stall_minutes):
    """Warnings for a RUNNING transport task: stalled logs, runaway TPC loopers (DESIGN.md section 4c)."""
    warns = []
    cwd = os.path.join(workdir, task.get("cwd", "") or "")
    files = [p for p in glob.glob(os.path.join(cwd, "sgn*")) if re.search(r"(\.log|serverlog|workerlog\d*|mergerlog)$", p)]
    if files:
        idle = (time.time() - max(os.path.getmtime(p) for p in files)) / 60.0
        if idle > stall_minutes:
            warns.append(f"   WARNING {task['name']}: no log activity for {idle:.0f} min (threshold {stall_minutes} min) -- possibly STALLED")
    server = os.path.join(cwd, "sgn_serverlog")
    if os.path.exists(server):
        rx = re.compile(r"Flat gas loopers: (\d+)")
        worst, neg = 0, False
        with open(server, errors="replace") as f:
            for line in f:
                m = rx.search(line)
                if m:
                    worst = max(worst, int(m.group(1)))
                if "Time constraint for loopers: -" in line:
                    neg = True
        if neg or worst > 5_000_000:
            warns.append(f"   ALARM {task['name']}: runaway TPC loopers (max {worst} requested{', negative time window' if neg else ''}).")
            warns.append("         This O2 bug never finishes (DESIGN.md section 4c): stop the run with run/killRun.sh and re-run")
            warns.append("         (runInjectedPol.sh applies the workaround automatically; or set INJPOL_LOOPERS=noGeomProtection).")
    return warns


def human_bytes(n):
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
        n /= 1024.0


def human_time(s):
    if s is None:
        return "n/a"
    s = int(round(s))
    return f"{s // 3600}h{(s % 3600) // 60:02d}m{s % 60:02d}s" if s >= 3600 else f"{s // 60}m{s % 60:02d}s"


def disk_report(workdir):
    cats = {"logs": 0, "AO2D": 0, "kinematics (Kine/MCHeader)": 0, "hits": 0, "digits": 0, "other ROOT (reco intermediates, geometry, ...)": 0, "other": 0}
    for root, _, files in os.walk(workdir):
        for fn in files:
            p = os.path.join(root, fn)
            try:
                size = os.lstat(p).st_size
            except OSError:
                continue
            if re.search(r"(\.log(_done|_time)?$|_serverlog$|_workerlog\d*$|_mergerlog$|\.txt$)", fn) or "/logs/" in p + "/":
                cats["logs"] += size
            elif fn.startswith("AO2D") and fn.endswith(".root"):
                cats["AO2D"] += size
            elif fn.endswith("_Kine.root") or fn.endswith("_MCHeader.root"):
                cats["kinematics (Kine/MCHeader)"] += size
            elif "Hits" in fn and fn.endswith(".root"):
                cats["hits"] += size
            elif re.search(r"dig", fn, re.I) and fn.endswith(".root"):
                cats["digits"] += size
            elif fn.endswith(".root"):
                cats["other ROOT (reco intermediates, geometry, ...)"] += size
            else:
                cats["other"] += size
    return cats


def report(workdir, targets, final=False):
    out = []
    tasks = load_tasks(workdir, targets)
    starts, t0 = start_times(workdir)
    now = dt.datetime.now()
    rows, done_walls = [], {}
    for t in tasks:
        state, wall = task_status(workdir, t)
        rows.append((t, state, wall))
        if state == "done" and wall is not None:
            done_walls.setdefault(base_name(t["name"]), []).append(wall)

    n_done = sum(1 for _, s, _ in rows if s == "done")
    n_tot = len(rows)
    running = [(t, s, w) for t, s, w in rows if s == "running/failed"]

    # expected wall time per task: average of finished tasks of the same kind (other TFs), else unknown
    def expected(t):
        v = done_walls.get(base_name(t["name"]))
        return sum(v) / len(v) if v else None

    known_done = sum(w for _, s, w in rows if s == "done" and w)
    rem, unknown, unknown_names = 0.0, 0, []
    for t, s, w in rows:
        if s == "done":
            continue
        e = expected(t)
        if e is None:
            unknown += 1
            unknown_names.append(t["name"])
        else:
            el = (now - starts[t["name"]]).total_seconds() if (s == "running/failed" and t["name"] in starts) else 0.0
            rem += max(e - el, 0.1 * e)
    frac_w = known_done / (known_done + rem) if (known_done + rem) > 0 and unknown == 0 else n_done / max(n_tot, 1)
    bar = "#" * int(40 * frac_w) + "." * (40 - int(40 * frac_w))
    elapsed = (now - t0).total_seconds() if t0 else None

    out.append(f"== InjectedPolO2MC workflow report: {workdir}")
    out.append(f"   targets {targets}, {n_tot} tasks: {n_done} done, {len(running)} running (or failed), {n_tot - n_done - len(running)} pending")
    out.append(f"   [{bar}] {100 * frac_w:5.1f}%   elapsed {human_time(elapsed)}")

    # event-level progress of the Geant4 stage (the dominant cost)
    health = []
    for t, s, w in rows:
        if base_name(t["name"]) == "sgnsim":
            fin, dist, expected = sim_events(workdir, t)
            total = expected if expected else dist
            line = f"   {t['name']}: {fin}/{expected if expected else '?'} events finished ({dist} distributed)"
            if s == "running/failed" and fin > 0 and t["name"] in starts:
                el = (now - starts[t["name"]]).total_seconds()
                line += f", {el / fin:.1f} s/event so far; transport of this TF ends in ~{human_time(el / fin * max(total - fin, 0))}"
                if not expected:
                    line += " (total unknown: counted from distributed events only)"
            elif w and fin > 0:
                line += f", wall {human_time(w)} = {w / fin:.2f} s/event (Geant4 + generator, {t.get('resources', {}).get('cpu', '?')} CPU booked)"
            out.append(line)
            if s == "running/failed":
                health += sim_health(workdir, t, STALL_MINUTES)
    out += health
    if rem > 0:
        extra = "" if unknown == 0 else f" + {unknown} task(s) with no timing reference yet ({', '.join(unknown_names[:4])}{', ...' if unknown > 4 else ''})"
        out.append(f"   ETA: <= ~{human_time(rem)}{extra}")
        out.append("        (sum of the wall times of finished tasks of the same kind in other timeframes; an UPPER bound, since")
        out.append("         timeframes run in parallel when CPUs allow)")
    elif n_done < n_tot:
        out.append("   ETA: not available yet -- it needs one finished timeframe (a single-TF run only has the per-event estimate of the transport above)")

    if final or n_done == n_tot:
        out.append("   -- wall time per task kind (sum over timeframes; tasks run in parallel, so the sum exceeds the elapsed time) --")
        agg = {}
        for t, s, w in rows:
            if w:
                agg[base_name(t["name"])] = agg.get(base_name(t["name"]), 0.0) + w
        for k, v in sorted(agg.items(), key=lambda kv: -kv[1])[:12]:
            out.append(f"      {k:<28s} {human_time(v):>10s}")

        nev = sum(sim_events(workdir, t)[0] for t, _, _ in rows if base_name(t["name"]) == "sgnsim")
        cats = disk_report(workdir)
        total = sum(cats.values())
        out.append(f"   -- disk usage of {workdir}: {human_bytes(total)} --")
        for k, v in cats.items():
            per = f"   ({human_bytes(v / nev)}/event)" if nev else ""
            out.append(f"      {k:<48s} {human_bytes(v):>10s}{per}")
        if nev:
            out.append(f"   simulated events (all timeframes): {nev}; total wall time {human_time(elapsed)} -> {elapsed / nev:.1f} s/event end-to-end" if elapsed else
                       f"   simulated events (all timeframes): {nev}")
            if nev < 500:
                out.append(f"   WARNING: only {nev} events. Per-event time and disk numbers are dominated by fixed per-timeframe costs")
                out.append("            (initialisation, geometry, CCDB objects, per-TF files) and OVERESTIMATE the cost of large productions.")
                out.append("            Use >= a few hundred events per TF (and several TFs) before planning long runs.")
        out.append("   (Number of events = collisions in the collision context: the requested -ns per TF plus 'orbits-early' pile-up")
        out.append("    collisions that O2DPG adds before the timeframe start, so it can exceed NSIGEVENTS x NTIMEFRAMES.)")

    text = "\n".join(out)
    if final:
        os.makedirs(os.path.join(workdir, "logs"), exist_ok=True)
        with open(os.path.join(workdir, "logs", "run_report.txt"), "w") as f:
            f.write(text + "\n")
    return text, n_done == n_tot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir")
    ap.add_argument("--targets", nargs="+", default=[TARGET_DEFAULT])
    ap.add_argument("--watch", type=int, default=0, help="refresh period in seconds (0 = one shot)")
    ap.add_argument("--final", action="store_true", help="full report, also written to <workdir>/logs/run_report.txt")
    ap.add_argument("--stall-minutes", type=int, default=STALL_MINUTES, help="warn when a running transport writes no log for this long")
    a = ap.parse_args()
    globals()["STALL_MINUTES"] = a.stall_minutes
    workdir = os.path.abspath(a.workdir.rstrip("/"))
    if not os.path.exists(os.path.join(workdir, "workflow.json")):
        sys.exit(f"no workflow.json in {workdir}")
    while True:
        text, finished = report(workdir, a.targets, a.final)
        if a.watch:
            print("\033[2J\033[H", end="")
        print(text, flush=True)
        if not a.watch or finished:
            break
        time.sleep(a.watch)


if __name__ == "__main__":
    main()
