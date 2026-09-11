"""Tests of the elliptic flow: v2 calibration and affine tangent transport."""

import numpy as np
import pytest

from ringtoy.flow import (EllipticFlow, anisotropy_ratio, apply_linear_map, hubble_matrix,
                          momentum_directions, shape_map_matrix, v2_from_ratio)
from ringtoy.geometry import RingPlacement, RingShape, build_ring, make_jet_frame


def test_v2_calibration_is_exact():
    v2, psi2 = 0.137, 0.41
    n = 4096
    phi = 2.0 * np.pi * np.arange(n) / n
    # Emitters uniform in azimuth about the vertex, with a longitudinal
    # component to check that it does not leak into the azimuthal map.
    x = np.stack([np.cos(phi), np.sin(phi), np.full(n, 0.35)], axis=1)
    p_hat, valid = momentum_directions(hubble_matrix(EllipticFlow(v2=v2, psi2=psi2, h_z=1.6)), x)
    assert np.all(valid)
    phi_p = np.arctan2(p_hat[:, 1], p_hat[:, 0])
    np.testing.assert_allclose(np.mean(np.cos(2.0 * (phi_p - psi2))), v2, atol=1e-12)
    np.testing.assert_allclose(v2_from_ratio(anisotropy_ratio(v2)), v2, atol=1e-15)


def test_isotropic_flow_is_identity():
    np.testing.assert_allclose(hubble_matrix(EllipticFlow(v2=0.0, psi2=0.73, h_z=1.0)),
                               np.eye(3), atol=1e-15)


def test_affine_tangent_transport_matches_finite_differences():
    frame = make_jet_frame(0.23, 0.9)
    curve = build_ring(RingShape(kind="fourier", radius=1.2,
                                 fourier_modes=((2, 0.15, 0.3), (3, 0.07, -0.5)),
                                 tilt=0.2, tilt_axis_angle=0.9, n_points=4000),
                       RingPlacement(distance_along_jet=1.4, offset_e1=0.6), frame)
    deformation = shape_map_matrix(EllipticFlow(v2=0.21, psi2=-0.7, strength=0.9, h_z=1.4))
    x_deformed = apply_linear_map(deformation, curve.x)
    analytic = apply_linear_map(deformation, curve.dx_dlam)
    # Centred periodic differences have error dlam^2/6 * |x'''| ~ 4e-7 * |x'''|,
    # and |x'''| reaches ~10 fm for these modes; atol sits just above that
    # bound, while an error in the transport rule would be O(1).
    numeric = (np.roll(x_deformed, -1, axis=0) - np.roll(x_deformed, 1, axis=0)) / (2.0 * curve.dlam)
    np.testing.assert_allclose(numeric, analytic, rtol=1e-4, atol=2e-5)


def test_shape_map_rejects_folding():
    with pytest.raises(ValueError):
        shape_map_matrix(EllipticFlow(v2=0.3, strength=-1.5))
