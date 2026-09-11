#!/usr/bin/env python3
"""
run_ring_toy.py -- Run a preset of the ring geometry toy and write its figures.

Usage (from ToyModels/RingGeometryToy/):
    python run_ring_toy.py --preset regimes
    python run_ring_toy.py --preset position --set shape.radius=1.5 --no-3d
    python run_ring_toy.py --list-presets

Needs numpy and matplotlib; plotly for the 3D scenes. Outputs go to
<output root>/<preset>_<UTC timestamp>/ with a provenance JSON, where the output
root is --output-dir, else /home/users/cicerodm/RingPol/RingGeometryToy/untrackedOutput.
See README.md for the physics and for what each preset studies.
"""

import argparse
import dataclasses
import getpass
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np

from ringtoy import plots2d
from ringtoy.scans import (aligned_counterpart, average_over_jet_azimuth, evaluate, evaluate_all,
                           scan_parameter, scan_position)
from ringtoy.observable import summary_from_sums
from ringtoy.scenarios import PRESETS, apply_override, get_preset, parse_override

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Default output root. Override with --output-dir when needed.
DEFAULT_OUTPUT_DIR = "/home/users/cicerodm/RingPol/RingGeometryToy/untrackedOutput"


def git_commit():
    """Commit hash with a -dirty suffix for uncommitted changes, or 'unknown'."""
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SCRIPT_DIR, capture_output=True,
                                text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=SCRIPT_DIR,
                                capture_output=True, text=True, check=True).stdout.strip()
        return commit + ("-dirty" if status else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _json_safe(value):
    """Convert numpy scalars/arrays and non-finite floats to JSON-compatible values."""
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _json_safe(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        # Strict JSON has no NaN; null marks an undefined ratio (empty selection).
        return float(value) if math.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def write_provenance(output_dir, preset, config, results_block, argv):
    record = {
        "git_commit": git_commit(),
        "run_by": getpass.getuser(),
        "run_timestamp": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(argv),
        "preset": preset,
        "config": config,
        "results": results_block,
    }
    with open(os.path.join(output_dir, "provenance.json"), "w", encoding="ascii") as handle:
        json.dump(_json_safe(record), handle, indent=2)


def run_regimes(config, output_dir, args):
    results = evaluate_all(config)
    stem = os.path.join(output_dir, "regimes")
    plots2d.plot_regime_overlays(results, stem)
    plots2d.plot_cos_gamma(results, f"{stem}_cos_gamma")
    plots2d.plot_jet_plane(results, f"{stem}_jet_plane")
    plots2d.plot_momentum_map(results, config.acceptance.eta_max, f"{stem}_momentum_map")
    if not args.no_3d:
        from ringtoy import viz3d
        viz3d.scene_figure(config, results, tube_radius=args.tube_radius).write_html(
            f"{stem}_scene.html", include_plotlyjs=True)
        strengths = np.linspace(0.0, max(config.flow.strength, 1e-9) * 2.0, 9)
        viz3d.strength_slider_figure(config, strengths, evaluate, regime="R3",
                                     tube_radius=args.tube_radius).write_html(
            f"{stem}_R3_strength_slider.html", include_plotlyjs=True)
    return {regime: res.summary for regime, res in results.items()}


def run_position(config, output_dir, args):
    distances = config.scan.values
    results = scan_position(config, distances)
    block = {}
    for regime, per_distance in results.items():
        for variable in ("delta_phi", "delta_theta"):
            curves = []
            for distance, res in zip(distances, per_distance):
                label = f"L = {distance:+.1f} fm" + (" (inward)" if distance < 0 else "")
                shade = 0.15 + 0.7 * abs(distance) / max(1e-9, max(abs(d) for d in distances))
                # Outward rings in blue, inward in red, the centred ring in black;
                # lighter shades are closer to the vertex.
                if distance == 0:
                    colour = "#000000"
                else:
                    colour = blend_towards_white("#cc0000" if distance < 0 else "#0033cc", shade)
                style = dict(color=colour, linestyle="--" if distance < 0 else "-")
                curves.append((label, getattr(res, variable), style))
            plots2d.plot_binned_overlay(curves, variable,
                                        os.path.join(output_dir, f"position_{regime}_{variable}"),
                                        title=plots2d.REGIME_TITLE[regime])
        block[regime] = {f"{d:+.2f}": res.summary for d, res in zip(distances, per_distance)}
    return block


def blend_towards_white(hex_color, shade):
    """Blend a colour towards white; shade in (0, 1], 1 is the full colour."""
    rgb = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)
    blended = 255.0 - shade * (255.0 - rgb)
    return "#" + "".join(f"{int(round(c)):02x}" for c in blended)


def run_misaligned(config, output_dir, args):
    misaligned = average_over_jet_azimuth(config, config.scan.n_jet_azimuths)
    aligned = average_over_jet_azimuth(aligned_counterpart(config), config.scan.n_jet_azimuths)
    block = {}
    for regime in misaligned:
        for index, variable in ((1, "delta_phi"), (2, "delta_theta")):
            curves = [("aligned (ring centred on the jet)", aligned[regime][index],
                       dict(color="#0033cc", linestyle="-")),
                      ("misaligned (ring centre fixed in the lab)", misaligned[regime][index],
                       dict(color="#cc0000", linestyle="--"))]
            plots2d.plot_binned_overlay(curves, variable,
                                        os.path.join(output_dir, f"misaligned_{regime}_{variable}"),
                                        title=plots2d.REGIME_TITLE[regime])
        s_mis = summary_from_sums(misaligned[regime][0])
        s_ali = summary_from_sums(aligned[regime][0])
        block[regime] = {"misaligned": s_mis, "aligned": s_ali,
                         "signal_ratio_misaligned_over_aligned":
                             s_mis.signal_per_total_weight / s_ali.signal_per_total_weight
                             if s_ali.signal_per_total_weight != 0.0 else None}
    return block


def run_parameter(config, output_dir, args):
    scan = config.scan
    sums = scan_parameter(config, scan.parameter, scan.values)
    xlabel = {"jet_eta": r"$\eta_{\mathrm{Jet}}$",
              "delta": r"$\delta = \phi_{\mathrm{Jet}} - \Psi_2$ [rad]",
              "flow.strength": "shape-map strength $s$"}.get(scan.parameter, scan.parameter)
    plots2d.plot_scan(np.asarray(scan.values), sums, xlabel,
                      os.path.join(output_dir, f"scan_{scan.parameter.replace('.', '_')}"))
    return {regime: [summary_from_sums(s) for s in sums_list]
            for regime, sums_list in sums.items()}


RUNNERS = {"none": run_regimes, "position": run_position,
           "misaligned": run_misaligned, "parameter": run_parameter}


def main(argv):
    parser = argparse.ArgumentParser(description="Ring geometry toy.")
    parser.add_argument("--preset", default="regimes", help="preset name (see --list-presets)")
    parser.add_argument("--set", dest="overrides", action="append", default=[],
                        metavar="KEY=VALUE", help="override a configuration field, e.g. flow.v2=0.1")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                        help="output root")
    parser.add_argument("--no-3d", action="store_true", help="skip the plotly scenes")
    parser.add_argument("--tube-radius", type=float, default=None,
                        help="draw a cosmetic torus of this radius [fm] in the 3D scenes")
    parser.add_argument("--list-presets", action="store_true")
    args = parser.parse_args(argv[1:])

    if args.list_presets:
        for name in sorted(PRESETS):
            print(name)
        return 0

    config = get_preset(args.preset)
    for text in args.overrides:
        key, value = parse_override(text)
        config = apply_override(config, key, value)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    output_dir = os.path.join(args.output_dir, f"{args.preset}_{stamp}")
    os.makedirs(output_dir, exist_ok=True)

    plots2d.apply_style()
    results_block = RUNNERS[config.scan.kind](config, output_dir, args)
    # The 3D scene of a scan preset shows its base configuration.
    if config.scan.kind != "none" and not args.no_3d:
        from ringtoy import viz3d
        viz3d.scene_figure(config, evaluate_all(config), tube_radius=args.tube_radius).write_html(
            os.path.join(output_dir, "base_scene.html"), include_plotlyjs=True)

    write_provenance(output_dir, args.preset, config, results_block, argv)
    print(f"outputs written to {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))