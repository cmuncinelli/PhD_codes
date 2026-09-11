"""Tests of the almond drawn for context."""

import numpy as np
import pytest

from ringtoy.almond import AlmondModel, contour, current_eccentricity, eccentricity, rms_widths
from ringtoy.flow import EllipticFlow, anisotropy_ratio


def test_toy_relation_gives_flow_anisotropy():
    # kappa2 = 1: the Gaussian pressure-gradient ratio a^2/b^2 is the flow ratio r.
    v2 = 0.173
    a, b = rms_widths(v2, AlmondModel(kappa2=1.0, rms_radius=2.1))
    np.testing.assert_allclose(a**2 / b**2, anisotropy_ratio(v2), rtol=1e-14)
    np.testing.assert_allclose((b**2 - a**2) / (b**2 + a**2), v2, rtol=1e-14)
    np.testing.assert_allclose(np.sqrt((a**2 + b**2) / 2.0), 2.1, rtol=1e-14)


def test_kappa2_scales_eccentricity_and_limits_it():
    np.testing.assert_allclose(eccentricity(0.07, 0.25), 0.28, rtol=1e-14)
    with pytest.raises(ValueError):
        eccentricity(0.3, 0.25)


def test_contour_is_rotated_ellipse():
    flow = EllipticFlow(v2=0.21, psi2=0.63)
    almond = AlmondModel(kappa2=0.8, rms_radius=1.7)
    a, b = rms_widths(flow.v2, almond)
    points = contour(flow, almond, 2.0)
    c, s = np.cos(flow.psi2), np.sin(flow.psi2)
    u = c * points[:, 0] + s * points[:, 1]
    v = -s * points[:, 0] + c * points[:, 1]
    np.testing.assert_allclose((u / (2.0 * a))**2 + (v / (2.0 * b))**2, 1.0, rtol=1e-12)
    np.testing.assert_allclose(points[:, 2], 0.0, atol=1e-15)


def test_evolving_almond_becomes_round_at_inversion_strength():
    v2 = 0.19
    almond = AlmondModel(kappa2=0.6, rms_radius=1.4, evolve_with_flow=True)
    a, b = rms_widths(v2, almond)
    r = anisotropy_ratio(v2)
    inversion = (b - a) / (a - b * r)
    assert inversion > 0.0
    np.testing.assert_allclose(current_eccentricity(EllipticFlow(v2=v2, strength=inversion), almond), 0.0,
                               atol=1e-14)
    assert current_eccentricity(EllipticFlow(v2=v2, strength=2.0 * inversion), almond) < 0.0
