"""Consistency and symmetry tests of the four regimes and of the bookkeeping."""

import numpy as np
import pytest

from ringtoy.flow import EllipticFlow
from ringtoy.geometry import RingPlacement, RingShape, build_ring, make_jet_frame
from ringtoy.observable import Acceptance, bin_measurement, measure, summarize
from ringtoy.regimes import REGIMES, PolarizationModel, build_state

SHAPE = RingShape(kind="fourier", radius=1.1, fourier_modes=((2, 0.15, 0.3), (3, 0.07, -0.5)),
                  tilt=0.2, tilt_axis_angle=0.9, n_points=541)
PLACEMENT = RingPlacement(distance_along_jet=1.4, offset_e1=0.6, offset_e2=-0.3)
FLOW = EllipticFlow(v2=0.18, psi2=0.35, strength=0.7, h_z=1.3)
POLARIZATION = PolarizationModel(magnitude=0.74)


def _state(regime, flow, eta=0.21, phi=0.8):
    return build_state(build_ring(SHAPE, PLACEMENT, make_jet_frame(eta, phi)), flow, POLARIZATION, regime)


@pytest.mark.parametrize("regime, flow", [
    ("R1", EllipticFlow(v2=0.18, psi2=0.35, strength=0.0, h_z=1.3)),
    ("R2", EllipticFlow(v2=0.0, psi2=0.35, strength=0.7, h_z=1.0)),
    ("R3", EllipticFlow(v2=0.0, psi2=0.35, strength=0.0, h_z=1.0)),
])
def test_regimes_reduce_to_reference(regime, flow):
    reference = _state("R0", FLOW)
    reduced = _state(regime, flow)
    np.testing.assert_allclose(reduced.x, reference.x, atol=1e-13)
    np.testing.assert_allclose(reduced.p_hat, reference.p_hat, atol=1e-13)
    np.testing.assert_allclose(reduced.P, reference.P, atol=1e-13)


@pytest.mark.parametrize("regime", REGIMES)
def test_rotation_covariance(regime):
    # Rotating the jet and the event plane together about the beam must leave
    # every relative quantity unchanged.
    alpha = 1.1
    rotated_flow = EllipticFlow(v2=FLOW.v2, psi2=FLOW.psi2 + alpha, strength=FLOW.strength, h_z=FLOW.h_z)
    acceptance = Acceptance(eta_max=0.9)
    original = measure(_state(regime, FLOW), acceptance)
    rotated = measure(_state(regime, rotated_flow, phi=0.8 + alpha), acceptance)
    np.testing.assert_allclose(rotated.r, original.r, atol=1e-12)
    np.testing.assert_allclose(rotated.delta_theta, original.delta_theta, atol=1e-12)
    np.testing.assert_allclose(np.cos(rotated.delta_phi), np.cos(original.delta_phi), atol=1e-12)
    np.testing.assert_allclose(np.sin(rotated.delta_phi), np.sin(original.delta_phi), atol=1e-12)
    s_original, s_rotated = summarize(original), summarize(rotated)
    np.testing.assert_allclose(s_rotated.ring_average, s_original.ring_average, rtol=1e-10)
    np.testing.assert_allclose(s_rotated.acceptance_fraction, s_original.acceptance_fraction, rtol=1e-10)


def test_signal_density_integrates_to_total_signal():
    measurement = measure(_state("R3", FLOW), Acceptance(eta_max=0.7))
    binned = bin_measurement(measurement.delta_phi, measurement, np.linspace(-np.pi, np.pi, 38))
    summary = summarize(measurement)
    np.testing.assert_allclose(np.sum(binned.signal_density() * binned.widths),
                               summary.signal_per_total_weight, rtol=1e-12)


def test_stretch_scaling_follows_line_stretching():
    scaled = build_state(build_ring(SHAPE, PLACEMENT, make_jet_frame(0.21, 0.8)), FLOW,
                         PolarizationModel(magnitude=0.5, stretch_scaling=True,
                                           weight_mode="deformed_arclength"), "R1")
    # With deformed arc-length weights, weight / dlam is the deformed speed and
    # the undeformed speed is recovered from an R0 state.
    reference = _state("R0", FLOW)
    stretch = scaled.weight / reference.weight
    np.testing.assert_allclose(np.linalg.norm(scaled.P, axis=1), 0.5 * stretch, rtol=1e-12)
