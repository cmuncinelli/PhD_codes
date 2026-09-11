"""
flow.py -- Elliptic ("anisotropic Hubble") flow: shape map and momentum map.

Part of the ring geometry toy (README.md holds the derivations).

The flow is the linear velocity field u(x) = H x about the vertex, with

    H = R_z(psi2) diag(1, r, h_z) R_z(-psi2),   r = (1 - v2) / (1 + v2).

It drives two independent operations:
  * shape map    : x' = (1 + s H) x, an affine map that turns a transverse
                   circle about the vertex into an ellipse elongated along
                   the event plane; tangents transform with the same matrix.
  * momentum map : p_hat = normalize(H x), emission along the local flow.
                   For emitters uniform in azimuth about the vertex this
                   gives exactly <cos 2(phi_p - psi2)> = v2.
"""

from dataclasses import dataclass

import numpy as np

from .geometry import DEGENERACY_EPS


@dataclass(frozen=True)
class EllipticFlow:
    """Parameters of the elliptic flow.

    v2       : momentum-space anisotropy produced by the momentum map for
               emitters uniform in azimuth about the vertex, |v2| < 1.
    psi2     : event-plane angle [rad].
    strength : s >= 0 in the shape map x' = (1 + s H) x (dimensionless:
               expansion rate times elapsed time, in units of the in-plane rate).
    h_z      : longitudinal expansion rate relative to the in-plane rate, > 0.
    """
    v2: float = 0.0
    psi2: float = 0.0
    strength: float = 0.0
    h_z: float = 1.0


def anisotropy_ratio(v2):
    """Out-of-plane to in-plane rate ratio r that produces a given v2."""
    if not -1.0 < v2 < 1.0:
        raise ValueError(f"v2 must satisfy |v2| < 1, got {v2}")
    return (1.0 - v2) / (1.0 + v2)


def v2_from_ratio(ratio):
    """Exact inverse of anisotropy_ratio (see README.md for the derivation)."""
    return (1.0 - ratio) / (1.0 + ratio)


def hubble_matrix(flow):
    """Symmetric velocity-gradient matrix H of the elliptic flow."""
    if not flow.h_z > 0.0:
        raise ValueError(f"h_z must be positive, got {flow.h_z}")
    ratio = anisotropy_ratio(flow.v2)
    c = np.cos(flow.psi2)
    s = np.sin(flow.psi2)
    rotation = np.array([[c, -s, 0.0],
                         [s, c, 0.0],
                         [0.0, 0.0, 1.0]])
    return rotation @ np.diag([1.0, ratio, flow.h_z]) @ rotation.T


def shape_map_matrix(flow):
    """Deformation gradient F = 1 + s H of the shape map."""
    ratio = anisotropy_ratio(flow.v2)
    # H has eigenvalues (1, r, h_z), so F has (1 + s, 1 + s r, 1 + s h_z).
    # A non-positive eigenvalue would fold the ring through the vertex, which
    # no expansion can do.
    eigenvalues = 1.0 + flow.strength * np.array([1.0, ratio, flow.h_z])
    if np.any(eigenvalues <= 0.0):
        raise ValueError("shape map is not orientation preserving for this strength")
    return np.eye(3) + flow.strength * hubble_matrix(flow)


def apply_linear_map(matrix, vectors):
    """Apply a 3x3 matrix to (N, 3) row vectors."""
    return vectors @ matrix.T


def momentum_directions(matrix, x):
    """Unit momentum directions normalize(matrix x), with a validity mask.

    Points mapped to (numerically) zero have no direction; they are returned
    as NaN and flagged invalid instead of being given an arbitrary direction.
    """
    p = apply_linear_map(matrix, x)
    norm = np.linalg.norm(p, axis=1)
    valid = norm > DEGENERACY_EPS
    p_hat = np.full_like(p, np.nan)
    p_hat[valid] = p[valid] / norm[valid, None]
    return p_hat, valid
