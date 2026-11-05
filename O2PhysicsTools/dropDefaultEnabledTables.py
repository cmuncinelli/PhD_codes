#!/usr/bin/env python3
"""Remove `enabledTables` blocks that only restate the defaults (every entry "-1") from DPL configuration JSONs.

Why: service tasks such as mult-cent-table (Common/Tools/Multiplicity/MultModule.h) and propagation-service
(PWGLF/Utils/strangenessBuilderModule.h) read `enabledTables` as a LabeledArray and look every table up BY LABEL.
When O2Physics renames or adds a table, a JSON written by an older version no longer has that label, and the task
segfaults in Init() (o2::framework::Array2D<int>::operator()). Seen on 2026-11-05 with mult-cent-table: the JSON had
39 tables (with FV0AOuterMults), the code 42 (FITExtraMults, CentFT0MAnchorCols, CentFT0MAnchorBCs, CentFT0MOuterAs).

A block whose values are all "-1" ("enable if some task needs it", the default) carries no information: dropping it
lets the task use its own, always up-to-date, table list. Blocks with any explicit 0/1 are reported and kept: those
must be updated by hand (or regenerated from Hyperloop).

Usage:
  dropDefaultEnabledTables.py <json or folder> [...]          dry run: report what would change
  dropDefaultEnabledTables.py <json or folder> [...] --apply  rewrite, keeping a backup in <folder>/preEnabledTablesDrop/
"""
import argparse
import datetime
import json
import os
import shutil
import sys

BACKUP_DIRNAME = "preEnabledTablesDrop"


def all_default(block):
    vals = block.get("values")
    if not isinstance(vals, list) or not vals:
        return False
    flat = [v for row in vals for v in (row if isinstance(row, list) else [row])]
    return all(str(v).strip() == "-1" for v in flat)


def process(path, apply):
    with open(path) as f:
        try:
            cfg = json.load(f)
        except json.JSONDecodeError as e:
            print(f"  SKIP {path}: not valid JSON ({e})")
            return 0
    changed = 0
    for device, opts in cfg.items():
        if not isinstance(opts, dict) or "enabledTables" not in opts:
            continue
        block = opts["enabledTables"]
        n = len(block.get("labels_rows", [])) if isinstance(block, dict) else 0
        if isinstance(block, dict) and all_default(block):
            print(f"  {os.path.basename(path)}: {device}.enabledTables ({n} tables, all -1) -> DROP (task defaults apply)")
            del opts["enabledTables"]
            changed += 1
        else:
            print(f"  {os.path.basename(path)}: {device}.enabledTables ({n} tables) has explicit values -> KEEP; "
                  f"check its labels against the current O2Physics table list by hand")
    if changed and apply:
        bdir = os.path.join(os.path.dirname(os.path.abspath(path)), BACKUP_DIRNAME)
        os.makedirs(bdir, exist_ok=True)
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(path, os.path.join(bdir, f"{os.path.basename(path)}.{stamp}"))
        with open(path, "w") as f:
            json.dump(cfg, f, indent=4)
            f.write("\n")
    return changed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    files = []
    for p in a.paths:
        if os.path.isdir(p):
            files += sorted(os.path.join(p, f) for f in os.listdir(p) if f.endswith(".json"))
        elif p.endswith(".json"):
            files.append(p)
    if not files:
        sys.exit("no .json files given")
    total = sum(process(f, a.apply) for f in files)
    print(f"{total} block(s) {'dropped' if a.apply else 'would be dropped (dry run; add --apply)'} in {len(files)} file(s)")


if __name__ == "__main__":
    main()
