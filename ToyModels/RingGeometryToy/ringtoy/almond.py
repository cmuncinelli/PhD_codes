"""
almond.py -- QGP-like almond drawn for context, consistent with the flow v2.

Part of the ring geometry toy (README.md holds the derivation). The almond
does not enter any observable.

A Gaussian almond with rms widths a (along the event plane) and b
(perpendicular) accelerates the fluid as (x/a^2, y/b^2), so the flow
anisotropy ratio of flow.py is r = a^2/b^2. With r = (1 - v2)/(1 + v2) this
gives an eccentricity eps2 = (b^2 - a^2)/(b^2 + a^2) = v2 exactly, i.e. a
response coefficient kappa2 = 1. Realistic hydrodynamics has kappa2 ~ 0.2-0.3,
so kappa2 is a parameter: eps2 = v2 / kappa2.
"""

from dataclasses import dataclass

import numpy as np

from .flow import anisotropy_ratio, apply_linear_map, shape_map_matrix


@dataclass(frozen=True)
class AlmondModel:
    """Parameters of the almond.

    show             : draw the almond in scenes.
    kappa2           : response coefficient, eps2 = v2 / kappa2 (1 = toy relation).
    rms_radius       : sqrt((a^2 + b^2) / 2) [fm].
    evolve_with_flow : carry the almond with the same shape map as the ring.
    n_points         : points per contour.
    """
    show: bool = True
    kappa2: float = 1.0
    rms_radius: float = 1.8
    evolve_with_flow: bool = False
    n_points: int = 181


def eccentricity(v2, kappa2):
    """Initial eccentricity eps2 = v2 / kappa2."""
    if not kappa2 > 0.0:
        raise ValueError(f"kappa2 must be positive, got {kappa2}")
    eps2 = v2 / kappa2
    if not -1.0 < eps2 < 1.0:
        raise ValueError(f"|v2 / kappa2| must be below 1, got {eps2}")
    return eps2


def rms_widths(v2, almond):
    """rms widths (a, b) along and perpendicular to the event plane [fm]."""
    eps2 = eccentricity(v2, almond.kappa2)
    # a^2 + b^2 = 2 R^2 and (b^2 - a^2) / (b^2 + a^2) = eps2.
    a = almond.rms_radius * np.sqrt(1.0 - eps2)
    b = almond.rms_radius * np.sqrt(1.0 + eps2)
    return a, b


def current_widths(flow, almond):
    """rms widths after the shape map when the almond evolves with the flow."""
    a, b = rms_widths(flow.v2, almond)
    if almond.evolve_with_flow:
        # The shape map scales the event-plane axes by 1 + s and 1 + s r.
        a *= 1.0 + flow.strength
        b *= 1.0 + flow.strength * anisotropy_ratio(flow.v2)
    return a, b


def current_eccentricity(flow, almond):
    """Eccentricity of the drawn almond (it changes sign as the flow inverts it)."""
    a, b = current_widths(flow, almond)
    return (b**2 - a**2) / (b**2 + a**2)


def contour(flow, almond, n_sigma):
    """Lab-frame (N, 3) points of the n_sigma contour at z = 0 [fm]."""
    a, b = rms_widths(flow.v2, almond)
    t = np.linspace(0.0, 2.0 * np.pi, almond.n_points)
    c, s = np.cos(flow.psi2), np.sin(flow.psi2)
    u = n_sigma * a * np.cos(t)
    v = n_sigma * b * np.sin(t)
    points = np.stack([c * u - s * v, s * u + c * v, np.zeros_like(t)], axis=1)
    if almond.evolve_with_flow:
        points = apply_linear_map(shape_map_matrix(flow), points)
    return points
