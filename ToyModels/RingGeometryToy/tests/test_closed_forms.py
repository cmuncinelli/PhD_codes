"""Tests of the closed-form results of the Ring Design Discussion note."""

import numpy as np
import pytest

from ringtoy.geometry import direction_from_eta_phi
from ringtoy.observable import (beam_referenced_basis, cos_sin_c_closed_form,
                                elliptic_average_cos_c, ring_directions)


def test_cos_c_sin_c_closed_form_matches_vector_algebra():
    rng = np.random.default_rng(20260911)
    n = 200
    eta_jet = rng.uniform(-0.9, 0.9, n)
    eta_lambda = rng.uniform(-0.9, 0.9, n)
    phi_jet = rng.uniform(-np.pi, np.pi, n)
    phi_lambda = rng.uniform(-np.pi, np.pi, n)
    for i in range(n):
        t_hat = direction_from_eta_phi(eta_jet[i], phi_jet[i])
        p_hat = direction_from_eta_phi(eta_lambda[i], phi_lambda[i])[None, :]
        n_hat, _, _, valid = ring_directions(t_hat, p_hat)
        assert valid[0]
        phi_hat, theta_hat = beam_referenced_basis(p_hat)
        cos_c, sin_c = cos_sin_c_closed_form(eta_jet[i], eta_lambda[i], phi_lambda[i] - phi_jet[i])
        np.testing.assert_allclose(np.dot(n_hat[0], phi_hat[0]), cos_c, atol=1e-12)
        np.testing.assert_allclose(np.dot(n_hat[0], theta_hat[0]), sin_c, atol=1e-12)


@pytest.mark.parametrize("eta_jet", [0.05, 0.2, 0.5, 0.9])
def test_elliptic_average(eta_jet):
    integrate = pytest.importorskip("scipy.integrate")
    integrand = lambda dphi: cos_sin_c_closed_form(eta_jet, 0.0, dphi)[0]
    numeric = integrate.quad(integrand, 0.0, 2.0 * np.pi, points=[np.pi], limit=400)[0] / (2.0 * np.pi)
    np.testing.assert_allclose(elliptic_average_cos_c(eta_jet), numeric, rtol=1e-7)


def _quadrupole_leakage(eta_jet, eta_lambda, delta):
    # <n_z sin 2(u + delta)> over a uniform azimuth grid u = phi_lambda - phi_jet.
    n = 4096
    u = 2.0 * np.pi * np.arange(n) / n
    t_hat = direction_from_eta_phi(eta_jet, 0.0)
    p_hat = direction_from_eta_phi(np.full(n, eta_lambda), u)
    n_hat, _, _, valid = ring_directions(t_hat, p_hat)
    assert np.all(valid)
    return np.mean(n_hat[:, 2] * np.sin(2.0 * (u + delta)))


def test_quadrupole_leakage_selection_rule():
    # No leakage when either pseudorapidity vanishes, for any orientation.
    for delta in (0.0, 0.37, np.pi / 4):
        assert abs(_quadrupole_leakage(0.0, 0.4, delta)) < 1e-12
        assert abs(_quadrupole_leakage(0.3, 0.0, delta)) < 1e-12
    # Otherwise the leakage scales as cos 2 delta.
    base = _quadrupole_leakage(0.3, 0.4, 0.0)
    assert abs(base) > 1e-3
    np.testing.assert_allclose(_quadrupole_leakage(0.3, 0.4, np.pi / 8) / base,
                               np.cos(np.pi / 4), rtol=1e-9)
    assert abs(_quadrupole_leakage(0.3, 0.4, np.pi / 4)) < 1e-12
