"""
views.py -- Everything a scene shows, computed once from a configuration.

Part of the ring geometry toy (README.md describes the scenes). This module
has no plotting dependency: the interactive explorer, the static scenes and
the movies all draw the same View, so the physics lives in one place.
"""

import dataclasses
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from . import almond as almond_model
from .observable import summary_from_sums
from .regimes import REGIMES
from .scans import configure_parameter, evaluate, evaluate_all, get_parameter
from .scenarios import apply_override

# Contours drawn for the almond, in units of its rms widths.
ALMOND_SIGMAS = (1.0, 2.0)


@dataclass(frozen=True)
class SweepSpec:
    """Integrated results versus one parameter, shown next to a scene.

    parameter : dotted configuration key or "delta" (see scans.configure_parameter).
    values    : swept values.
    n_points  : ring resolution used for the sweep. Integrated results
                converge spectrally in the number of points, so a coarser
                ring keeps sweeps fast without visible loss.
    """
    parameter: str = "jet_phi"
    values: tuple = tuple(np.linspace(-np.pi, np.pi, 41))
    n_points: int = 240


@dataclass(frozen=True)
class View:
    """Results of one configuration, ready to draw.

    results           : {regime: RegimeResult}.
    almond_contours   : ((n_sigma, (N, 3) points), ...); empty when hidden.
    almond_eccentricity : eccentricity of the drawn almond (NaN when hidden).
    sweep             : SweepSpec or None.
    sweep_summaries   : {regime: (RingSummary, ...)} or None.
    sweep_value       : current value of the swept parameter, or None.
    """
    config: object
    results: dict
    almond_contours: tuple
    almond_eccentricity: float
    sweep: object
    sweep_summaries: dict
    sweep_value: float


@lru_cache(maxsize=64)
def _sweep_summaries(neutral_config, parameter, values, regimes):
    # Cached on a configuration whose swept parameter is pinned to a fixed
    # value, so animating that parameter reuses one sweep instead of
    # recomputing it for every frame.
    summaries = {regime: [] for regime in regimes}
    for value in values:
        point = configure_parameter(neutral_config, parameter, value)
        for regime in regimes:
            summaries[regime].append(summary_from_sums(evaluate(point, regime).sums))
    return {regime: tuple(items) for regime, items in summaries.items()}


def sweep_summaries(config, sweep, regimes=REGIMES):
    """Integrated summaries per regime along a sweep (cached)."""
    values = tuple(float(v) for v in sweep.values)
    neutral = configure_parameter(config, sweep.parameter, values[0])
    neutral = apply_override(neutral, "shape.n_points", min(config.shape.n_points, sweep.n_points))
    return _sweep_summaries(neutral, sweep.parameter, values, tuple(regimes))


def compute_view(config, regimes=REGIMES, sweep=None):
    """Evaluate a configuration and everything drawn next to it."""
    results = evaluate_all(config, regimes)
    if config.almond.show:
        contours = tuple((n_sigma, almond_model.contour(config.flow, config.almond, n_sigma))
                         for n_sigma in ALMOND_SIGMAS)
        eccentricity = almond_model.current_eccentricity(config.flow, config.almond)
    else:
        contours = ()
        eccentricity = float("nan")
    if sweep is not None:
        summaries = sweep_summaries(config, sweep, regimes)
        value = float(get_parameter(config, sweep.parameter))
    else:
        summaries = None
        value = None
    return View(config=config, results=results, almond_contours=contours,
                almond_eccentricity=eccentricity, sweep=sweep,
                sweep_summaries=summaries, sweep_value=value)


def mean_cos_gamma(result):
    """Weighted mean of cos(gamma) over accepted elements (NaN if none)."""
    m = result.measurement
    sel = m.accepted & np.isfinite(m.cos_gamma)
    total = np.sum(m.weight[sel])
    return float(np.sum(m.weight[sel] * m.cos_gamma[sel]) / total) if total > 0.0 else float("nan")


def readouts(view):
    """Per-regime numbers quoted next to a scene."""
    rows = {}
    for regime, result in view.results.items():
        summary = result.summary
        rows[regime] = dict(ring_average=summary.ring_average,
                            signal_per_total_weight=summary.signal_per_total_weight,
                            acceptance_fraction=summary.acceptance_fraction,
                            mean_cos_gamma=mean_cos_gamma(result))
    return rows


def highlight_mask(result, bin_index):
    """Accepted elements that fall in one Delta phi bin (all False for None)."""
    m = result.measurement
    if bin_index is None:
        return np.zeros(len(m.r), dtype=bool)
    edges = result.delta_phi.edges
    return m.accepted & (m.delta_phi >= edges[bin_index]) & (m.delta_phi < edges[bin_index + 1])


def with_parameter(config, parameter, value):
    """Configuration with a parameter set, for movies and sweeps."""
    return configure_parameter(config, parameter, float(value))


def replace_almond(config, **changes):
    return dataclasses.replace(config, almond=dataclasses.replace(config.almond, **changes))
