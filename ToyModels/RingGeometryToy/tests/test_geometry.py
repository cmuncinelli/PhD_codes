"""Geometry and sign-convention tests for the ring geometry toy."""

import numpy as np
import pytest

from ringtoy.flow import EllipticFlow
from ringtoy.geometry import RingPlacement, RingShape, build_ring, make_jet_frame
from ringtoy.observable import measure, ring_directions
from ringtoy.regimes import PolarizationModel, build_state

# Generic, non-special values: symmetric choices can hide sign or factor errors.
ETA_JET = 0.37
PHI_JET = -1.2


def test_jet_frame_is_right_handed_orthonormal():
    frame = make_jet_frame(ETA_JET, PHI_JET)
    basis = np.stack([frame.e1, frame.e2, frame.t_hat])
    np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=1e-14)
    np.testing.assert_allclose(np.cross(frame.e1, frame.e2), frame.t_hat, atol=1e-14)


def test_sign_convention_right_handed_ring():
    # Jet along +x, point at +y: n_hat must be +z.
    n_hat, _, _, valid = ring_directions(np.array([1.0, 0.0, 0.0]), np.array([[0.0, 1.0, 0.0]]))
    assert valid[0]
    np.testing.assert_allclose(n_hat[0], [0.0, 0.0, 1.0], atol=1e-15)

    # The default circulation is right-handed about the jet, so r = +f.
    frame = make_jet_frame(0.0, 0.0)
    curve = build_ring(RingShape(kind="circle", radius=0.5, n_points=64),
                       RingPlacement(distance_along_jet=1.0), frame)
    state = build_state(curve, EllipticFlow(), PolarizationModel(magnitude=0.8), "R0")
    np.testing.assert_allclose(measure(state).r, 0.8, atol=1e-14)


@pytest.mark.parametrize("distance", [2.3, -1.7, 0.0])
def test_coaxial_circle_gives_full_ring_value(distance):
    frame = make_jet_frame(ETA_JET, PHI_JET)
    curve = build_ring(RingShape(kind="circle", radius=1.3, n_points=257),
                       RingPlacement(distance_along_jet=distance), frame)
    state = build_state(curve, EllipticFlow(), PolarizationModel(magnitude=0.62), "R0")
    measurement = measure(state)
    assert measurement.n_undefined == 0
    np.testing.assert_allclose(measurement.r, 0.62, atol=1e-12)
    # Material weights of a circle add up to its circumference.
    np.testing.assert_allclose(np.sum(state.weight), 2.0 * np.pi * 1.3, rtol=1e-12)


def test_ellipse_ring_value_is_rho_dpsi_ds():
    a, b = 1.9, 0.8
    frame = make_jet_frame(ETA_JET, PHI_JET)
    curve = build_ring(RingShape(kind="ellipse", radius=a, radius_minor=b,
                                 ellipse_angle=0.3, n_points=301),
                       RingPlacement(distance_along_jet=1.1), frame)
    state = build_state(curve, EllipticFlow(), PolarizationModel(magnitude=1.0), "R0")
    lam = curve.lam
    rho = np.sqrt(a**2 * np.cos(lam)**2 + b**2 * np.sin(lam)**2)
    ds_dlam = np.sqrt(a**2 * np.sin(lam)**2 + b**2 * np.cos(lam)**2)
    # rho^2 dpsi/dlam = u v' - v u' = a b, invariant under the in-plane rotation.
    expected = a * b / (rho * ds_dlam)
    np.testing.assert_allclose(measure(state).r, expected, atol=1e-12)


def test_reversed_circulation_flips_ring_value():
    frame = make_jet_frame(ETA_JET, PHI_JET)
    placement = RingPlacement(distance_along_jet=0.9, offset_e1=0.2)
    shape = RingShape(kind="fourier", radius=1.0, fourier_modes=((2, 0.2, 0.4),), n_points=199)
    reversed_shape = RingShape(kind="fourier", radius=1.0, fourier_modes=((2, 0.2, 0.4),),
                               circulation=-1, n_points=199)
    r_plus = measure(build_state(build_ring(shape, placement, frame),
                                 EllipticFlow(), PolarizationModel(), "R0")).r
    r_minus = measure(build_state(build_ring(reversed_shape, placement, frame),
                                  EllipticFlow(), PolarizationModel(), "R0")).r
    np.testing.assert_allclose(r_minus, -r_plus, atol=1e-14)


def test_fourier_ring_rejects_non_positive_radius():
    frame = make_jet_frame(ETA_JET, PHI_JET)
    with pytest.raises(ValueError):
        build_ring(RingShape(kind="fourier", radius=1.0, fourier_modes=((3, 1.2, 0.0),)),
                   RingPlacement(), frame)
