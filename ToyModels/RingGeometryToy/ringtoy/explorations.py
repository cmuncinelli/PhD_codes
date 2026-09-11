"""
explorations.py -- Movies and checks that answer the geometric questions of the toy.

Part of the ring geometry toy (README.md lists the questions). Every function
writes into its own run directory with provenance, and returns what it drew
so that the notebook can display it.
"""

import dataclasses

import numpy as np

from . import plots2d
from .geometry import RingPlacement, RingShape
from .observable import ellipse_ring_average, summary_from_sums
from .output import make_run_dir, write_provenance
from .scans import evaluate, scan_parameter
from .scenarios import get_preset
from .views import SweepSpec
from .viz3d import FigureOptions, animated_figure


def _movie(label, config, parameter, values, options, output_root, command):
    sweep = SweepSpec(parameter=parameter, values=tuple(float(v) for v in values))
    fig = animated_figure(config, parameter, values, options, sweep=sweep)
    run_dir = make_run_dir(label, output_root)
    fig.write_html(f"{run_dir}/{label}.html", include_plotlyjs=True, auto_play=False)
    write_provenance(run_dir, label, config, {"parameter": parameter, "values": list(values)}, command)
    return fig, run_dir


def inward_outward_morph(config=None, distances=np.linspace(-4.0, 4.0, 33), options=FigureOptions(),
                         output_root=None):
    """Ring distance along the jet from -4 fm (inward) to +4 fm (outward).

    Expected: the cone on the momentum sphere flips from -t_hat to +t_hat and
    the Delta phi peaks migrate from +-pi towards 0.
    """
    config = config if config is not None else get_preset("position")
    return _movie("morph_inward_outward", config, "placement.distance_along_jet", distances, options,
                  output_root, "explorations.inward_outward_morph")


def misalignment_movie(config=None, n_frames=48, options=FigureOptions(), output_root=None):
    """Ring centre fixed in the lab while the jet azimuth turns once around."""
    config = config if config is not None else get_preset("misaligned")
    values = np.linspace(0.0, 2.0 * np.pi, n_frames, endpoint=False)
    return _movie("movie_misalignment", config, "jet_phi", values, options, output_root,
                  "explorations.misalignment_movie")


def eta_jet_movie(config=None, values=np.linspace(-1.5, 1.5, 31), options=FigureOptions(), output_root=None):
    """Jet pseudorapidity sweep: the ring cone crosses the eta window."""
    config = config if config is not None else get_preset("eta_jet")
    return _movie("movie_eta_jet", config, "jet_eta", values, options, output_root,
                  "explorations.eta_jet_movie")


def evolution_movie(config=None, strengths=None, options=FigureOptions(), output_root=None):
    """Shape-map strength with ring and almond evolving together.

    The almond starts elongated out of plane, becomes round at
    s = (b - a) / (a - b r), and ends elongated in plane.
    """
    config = config if config is not None else get_preset("strength")
    config = dataclasses.replace(config, almond=dataclasses.replace(config.almond, evolve_with_flow=True))
    if strengths is None:
        strengths = np.linspace(0.0, 2.5, 26)
    return _movie("movie_evolution", config, "flow.strength", strengths, options, output_root,
                  "explorations.evolution_movie")


def delta_test(config=None, values=np.linspace(0.0, np.pi, 25), output_root=None):
    """Does the ring response depend on the jet angle to the event plane?

    Returns {regime: relative variation (max - min) / |mean|} of the
    integrated ring average; a small number means the idea of selecting jets
    relative to Psi2 can be dropped.
    """
    config = config if config is not None else get_preset("delta")
    sums = scan_parameter(config, "delta", values)
    variation = {}
    for regime, sums_list in sums.items():
        averages = np.array([summary_from_sums(s).ring_average for s in sums_list])
        mean = np.mean(averages)
        variation[regime] = float((np.max(averages) - np.min(averages)) / abs(mean)) if mean != 0.0 else float("nan")
    run_dir = make_run_dir("delta_test", output_root)
    plots2d.apply_style()
    plots2d.plot_scan(np.asarray(values), sums, r"$\delta = \phi_{\mathrm{Jet}} - \Psi_2$ [rad]",
                      f"{run_dir}/delta_test", title="relative variation: " + ", ".join(
                          f"{k} {v:.3f}" for k, v in variation.items()))
    write_provenance(run_dir, "delta_test", config, {"relative_variation": variation}, "explorations.delta_test")
    return variation, run_dir


def ellipticity_check(ratios=np.linspace(0.1, 1.0, 19), major=1.7, output_root=None):
    """Toy ring average of a coaxial ellipse versus (b/a) K(k) / E(k)."""
    base = get_preset("regimes")
    toy, analytic = [], []
    for ratio in ratios:
        config = dataclasses.replace(
            base, shape=RingShape(kind="ellipse", radius=major, radius_minor=major * ratio,
                                  ellipse_angle=0.37, n_points=2048),
            placement=RingPlacement(distance_along_jet=1.3))
        toy.append(summary_from_sums(evaluate(config, "R0").sums).ring_average)
        analytic.append(float(ellipse_ring_average(major, major * ratio)))
    run_dir = make_run_dir("ellipticity_check", output_root)
    plots2d.apply_style()
    plots2d.plot_ellipticity_check(ratios, toy, analytic, f"{run_dir}/ellipticity_check")
    max_difference = float(np.max(np.abs(np.asarray(toy) - np.asarray(analytic))))
    write_provenance(run_dir, "ellipticity_check", base,
                     {"ratios": list(ratios), "toy": toy, "analytic": analytic,
                      "max_abs_difference": max_difference}, "explorations.ellipticity_check")
    return max_difference, run_dir


EXPLORATIONS = {
    "morph": inward_outward_morph,
    "misalignment": misalignment_movie,
    "eta_jet_movie": eta_jet_movie,
    "evolution": evolution_movie,
    "delta_test": delta_test,
    "ellipticity": ellipticity_check,
}
