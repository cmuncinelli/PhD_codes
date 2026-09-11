#!/usr/bin/env python3
"""
run_ring_toy.py -- Run a preset or an exploration of the ring geometry toy.

Usage (from ToyModels/RingGeometryToy/):
    python run_ring_toy.py --preset regimes
    python run_ring_toy.py --preset position --set shape.radius=1.5 --no-3d
    python run_ring_toy.py --explore morph
    python run_ring_toy.py --list

Needs numpy and matplotlib; plotly for the 3D scenes and movies; scipy for the
ellipticity check. Outputs go to <output root>/<label>_<UTC timestamp>/ with a
provenance JSON; the output root is --output-dir, else the path set in
ringtoy/output.py. The interactive explorer lives in RingGeometryToy.ipynb.
See README.md for the physics and for what each preset studies.
"""

import argparse
import os
import sys

import numpy as np

from ringtoy import plots2d
from ringtoy.observable import summary_from_sums
from ringtoy.output import DEFAULT_OUTPUT_DIR, make_run_dir, write_provenance
from ringtoy.scans import (aligned_counterpart, average_over_jet_azimuth, evaluate_all, scan_parameter,
                           scan_position)
from ringtoy.scenarios import PRESETS, apply_override, get_preset, parse_override


def blend_towards_white(hex_color, shade):
    """Blend a colour towards white; shade in (0, 1], 1 is the full colour."""
    rgb = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)
    blended = 255.0 - shade * (255.0 - rgb)
    return "#" + "".join(f"{int(round(c)):02x}" for c in blended)


def run_regimes(config, run_dir, args):
    results = evaluate_all(config)
    stem = os.path.join(run_dir, "regimes")
    plots2d.plot_regime_overlays(results, stem)
    plots2d.plot_cos_gamma(results, f"{stem}_cos_gamma")
    plots2d.plot_jet_plane(results, f"{stem}_jet_plane")
    plots2d.plot_momentum_map(results, config.acceptance.eta_max, f"{stem}_momentum_map")
    if not args.no_3d:
        from ringtoy.views import SweepSpec
        from ringtoy.viz3d import FigureOptions, animated_figure
        strengths = np.linspace(0.0, max(config.flow.strength, 0.5) * 2.0, 21)
        sweep = SweepSpec(parameter="flow.strength", values=tuple(strengths))
        animated_figure(config, "flow.strength", strengths, FigureOptions(tube_radius=args.tube_radius),
                        sweep=sweep).write_html(f"{stem}_strength_movie.html", include_plotlyjs=True)
    return {regime: res.summary for regime, res in results.items()}


def run_position(config, run_dir, args):
    distances = config.scan.values
    results = scan_position(config, distances)
    block = {}
    reach = max(1e-9, max(abs(d) for d in distances))
    for regime, per_distance in results.items():
        for variable in ("delta_phi", "delta_theta"):
            curves = []
            for distance, res in zip(distances, per_distance):
                label = f"L = {distance:+.1f} fm" + (" (inward)" if distance < 0 else "")
                # Outward rings in blue, inward in red, the centred ring in black;
                # lighter shades are closer to the vertex.
                shade = 0.15 + 0.7 * abs(distance) / reach
                if distance == 0:
                    colour = "#000000"
                else:
                    colour = blend_towards_white("#cc0000" if distance < 0 else "#0033cc", shade)
                curves.append((label, getattr(res, variable), dict(color=colour,
                                                                   linestyle="--" if distance < 0 else "-")))
            plots2d.plot_binned_overlay(curves, variable, os.path.join(run_dir, f"position_{regime}_{variable}"),
                                        title=plots2d.REGIME_TITLE[regime])
        block[regime] = {f"{d:+.2f}": res.summary for d, res in zip(distances, per_distance)}
    return block


def run_misaligned(config, run_dir, args):
    misaligned = average_over_jet_azimuth(config, config.scan.n_jet_azimuths)
    aligned = average_over_jet_azimuth(aligned_counterpart(config), config.scan.n_jet_azimuths)
    block = {}
    for regime in misaligned:
        for index, variable in ((1, "delta_phi"), (2, "delta_theta")):
            curves = [("aligned (ring centred on the jet)", aligned[regime][index],
                       dict(color="#0033cc", linestyle="-")),
                      ("misaligned (ring centre fixed in the lab)", misaligned[regime][index],
                       dict(color="#cc0000", linestyle="--"))]
            plots2d.plot_binned_overlay(curves, variable, os.path.join(run_dir, f"misaligned_{regime}_{variable}"),
                                        title=plots2d.REGIME_TITLE[regime])
        s_mis = summary_from_sums(misaligned[regime][0])
        s_ali = summary_from_sums(aligned[regime][0])
        block[regime] = {"misaligned": s_mis, "aligned": s_ali,
                         "signal_ratio_misaligned_over_aligned":
                             s_mis.signal_per_total_weight / s_ali.signal_per_total_weight
                             if s_ali.signal_per_total_weight != 0.0 else None}
    return block


def run_parameter(config, run_dir, args):
    scan = config.scan
    sums = scan_parameter(config, scan.parameter, scan.values)
    xlabel = {"jet_eta": r"$\eta_{\mathrm{Jet}}$",
              "delta": r"$\delta = \phi_{\mathrm{Jet}} - \Psi_2$ [rad]",
              "flow.strength": "shape-map strength $s$"}.get(scan.parameter, scan.parameter)
    plots2d.plot_scan(np.asarray(scan.values), sums, xlabel,
                      os.path.join(run_dir, f"scan_{scan.parameter.replace('.', '_')}"))
    return {regime: [summary_from_sums(s) for s in sums_list] for regime, sums_list in sums.items()}


RUNNERS = {"none": run_regimes, "position": run_position,
           "misaligned": run_misaligned, "parameter": run_parameter}


def main(argv):
    parser = argparse.ArgumentParser(description="Ring geometry toy.")
    parser.add_argument("--preset", default=None, help="preset name (see --list)")
    parser.add_argument("--explore", default=None, help="exploration name (see --list)")
    parser.add_argument("--set", dest="overrides", action="append", default=[],
                        metavar="KEY=VALUE", help="override a preset field, e.g. flow.v2=0.1")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR,
                        help=f"output root (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--no-3d", action="store_true", help="skip the plotly scenes")
    parser.add_argument("--tube-radius", type=float, default=0.0,
                        help="cosmetic torus radius [fm] in the 3D scenes (0: none)")
    parser.add_argument("--list", action="store_true", help="list presets and explorations")
    args = parser.parse_args(argv[1:])

    if args.list:
        from ringtoy.explorations import EXPLORATIONS
        print("presets:      " + " ".join(sorted(PRESETS)))
        print("explorations: " + " ".join(sorted(EXPLORATIONS)))
        return 0

    if args.explore is not None:
        # Explorations use their own base presets and write their own run directories.
        from ringtoy.explorations import EXPLORATIONS
        if args.explore not in EXPLORATIONS:
            parser.error(f"unknown exploration '{args.explore}'")
        outcome = EXPLORATIONS[args.explore](output_root=args.output_dir)
        print(f"outputs written to {outcome[-1]}")
        return 0

    preset = args.preset if args.preset is not None else "regimes"
    config = get_preset(preset)
    for text in args.overrides:
        key, value = parse_override(text)
        config = apply_override(config, key, value)

    run_dir = make_run_dir(preset, args.output_dir)
    plots2d.apply_style()
    results = RUNNERS[config.scan.kind](config, run_dir, args)
    if config.scan.kind != "none" and not args.no_3d:
        # A scan preset also gets the static scene of its base configuration.
        from ringtoy.views import compute_view
        from ringtoy.viz3d import FigureOptions, make_figure
        make_figure(compute_view(config), FigureOptions(tube_radius=args.tube_radius)).write_html(
            os.path.join(run_dir, "base_scene.html"), include_plotlyjs=True)

    write_provenance(run_dir, preset, config, results, " ".join(argv))
    print(f"outputs written to {run_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
