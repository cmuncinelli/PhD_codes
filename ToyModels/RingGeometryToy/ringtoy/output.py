"""
output.py -- Output location and provenance records.

Shared by run_ring_toy.py and the notebook exports, so that both write to the
same place with the same provenance. See README.md.
"""

import dataclasses
import getpass
import json
import math
import os
import subprocess
from datetime import datetime, timezone

import numpy as np

# Results accompany the other RingPol outputs. The driver's --output-dir and
# the output_root argument of the notebook exports override it.
DEFAULT_OUTPUT_DIR = "/home/users/cicerodm/RingPol/RingGeometryToy"

# ToyModels/RingGeometryToy/, whose git commit identifies the code of a run.
CODE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def git_commit():
    """Commit hash of the code, with -dirty for uncommitted changes, or 'unknown'."""
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=CODE_DIR, capture_output=True,
                                text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=CODE_DIR,
                                capture_output=True, text=True, check=True).stdout.strip()
        return commit + ("-dirty" if status else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def json_safe(value):
    """Convert dataclasses, numpy values and non-finite floats for strict JSON."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return json_safe(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        # Strict JSON has no NaN; null marks an undefined ratio (empty selection).
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def make_run_dir(label, output_root=None):
    """Create <output_root>/<label>_<UTC timestamp>/ and return its path."""
    root = DEFAULT_OUTPUT_DIR if output_root is None else output_root
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    run_dir = os.path.join(root, f"{label}_{stamp}")
    os.makedirs(run_dir, exist_ok=True)
    return run_dir


def write_provenance(run_dir, label, config, results, command):
    """Write provenance.json: code commit, user, time, command, configuration, results."""
    record = {
        "git_commit": git_commit(),
        "run_by": getpass.getuser(),
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "label": label,
        "config": config,
        "results": results,
    }
    with open(os.path.join(run_dir, "provenance.json"), "w", encoding="ascii") as handle:
        json.dump(json_safe(record), handle, indent=2)
