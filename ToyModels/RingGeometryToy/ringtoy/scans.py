"""
scans.py -- Evaluation of the regimes and of the parameter scans.

Part of the ring geometry toy (README.md describes the scans).
"""

import dataclasses
from dataclasses import dataclass

import numpy as np

from .geometry import RingPlacement, build_ring, make_jet_frame
from .observable import bin_measurement, integrate, measure, summary_from_sums
from .regimes import REGIMES, build_state
from .scenarios import apply_override


@dataclass(frozen=True)
class RegimeResult:
    """Everything computed for one regime of one configuration."""
    state: object
    measurement: object
    sums: object
    delta_phi: object
    delta_theta: object

    @property
    def summary(self):
        return summary_from_sums(self.sums)


def delta_phi_edges(binning):
    return np.linspace(-np.pi, np.pi, binning.n_delta_phi + 1)


def delta_theta_edges(binning):
    return np.linspace(0.0, np.pi, binning.n_delta_theta + 1)


def evaluate(config, regime):
    """Build the ring, apply one regime, and measure it."""
    frame = make_jet_frame(config.jet_eta, config.jet_phi)
    curve = build_ring(config.shape, config.placement, frame)
    state = build_state(curve, config.flow, config.polarization, regime)
    measurement = measure(state, config.acceptance)
    return RegimeResult(
        state=state, measurement=measurement, sums=integrate(measurement),
        delta_phi=bin_measurement(measurement.delta_phi, measurement, delta_phi_edges(config.binning)),
        delta_theta=bin_measurement(measurement.delta_theta, measurement, delta_theta_edges(config.binning)))


def evaluate_all(config, regimes=REGIMES):
    return {regime: evaluate(config, regime) for regime in regimes}


def configure_parameter(config, parameter, value):
    """Set a scanned parameter; "delta" places the jet at psi2 + value."""
    if parameter == "delta":
        return dataclasses.replace(config, jet_phi=config.flow.psi2 + value)
    return apply_override(config, parameter, value)


def scan_parameter(config, parameter, values, regimes=REGIMES):
    """Integrated sums per regime for each scanned value."""
    sums = {regime: [] for regime in regimes}
    for value in values:
        point = configure_parameter(config, parameter, value)
        for regime in regimes:
            sums[regime].append(evaluate(point, regime).sums)
    return sums


def scan_position(config, distances, regimes=REGIMES):
    """Full results per regime for each distance of the ring along the jet."""
    results = {regime: [] for regime in regimes}
    for distance in distances:
        point = dataclasses.replace(config, placement=RingPlacement(distance_along_jet=distance))
        for regime in regimes:
            results[regime].append(evaluate(point, regime))
    return results


def average_over_jet_azimuth(config, n_azimuths, regimes=REGIMES):
    """Add results over jet azimuths uniform in [0, 2 pi).

    Returns {regime: (sums, delta_phi, delta_theta)}. Sums and histograms are
    added before any ratio is formed, so each jet orientation contributes in
    proportion to its ring weight.
    """
    totals = {}
    for k in range(n_azimuths):
        point = dataclasses.replace(config, jet_phi=2.0 * np.pi * k / n_azimuths)
        for regime in regimes:
            result = evaluate(point, regime)
            if regime not in totals:
                totals[regime] = (result.sums, result.delta_phi, result.delta_theta)
            else:
                sums, dphi, dtheta = totals[regime]
                totals[regime] = (sums + result.sums, dphi + result.delta_phi,
                                  dtheta + result.delta_theta)
    return totals


def aligned_counterpart(config):
    """Same ring distance from the vertex, but centred on the jet axis."""
    if config.placement.center_lab is None:
        return config
    distance = float(np.linalg.norm(config.placement.center_lab))
    return dataclasses.replace(config, placement=RingPlacement(distance_along_jet=distance))
