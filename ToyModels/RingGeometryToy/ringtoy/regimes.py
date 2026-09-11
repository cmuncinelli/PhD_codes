"""
regimes.py -- Polarized ring states for the four distortion regimes.

Part of the ring geometry toy (README.md holds the regime definitions).

    R0 : ring undeformed, momenta along the position direction
    R1 : ring deformed by the shape map, momenta along the new positions
    R2 : ring undeformed, momenta deflected by the momentum map
    R3 : ring deformed and momenta deflected from the deformed positions
"""

from dataclasses import dataclass

import numpy as np

from .flow import apply_linear_map, hubble_matrix, momentum_directions, shape_map_matrix
from .geometry import JetFrame

REGIMES = ("R0", "R1", "R2", "R3")


@dataclass(frozen=True)
class PolarizationModel:
    """How the ring carries polarization.

    magnitude       : f in [0, 1]; each ring element has |P| = f unless
                      stretch_scaling is on.
    stretch_scaling : scale |P| by the local line stretching |F tau| / |tau|
                      (vortex stretching). Off by default.
                      TODO: decide whether stretch scaling belongs in the default model
    weight_mode     : "material" keeps the weight of each element fixed under
                      the shape map (number of emitters conserved);
                      "deformed_arclength" re-weights by the deformed length.
                      TODO: revisit ring rarefaction (density of emitters along a stretched ring)
    """
    magnitude: float = 1.0
    stretch_scaling: bool = False
    weight_mode: str = "material"


@dataclass(frozen=True)
class RingState:
    """A polarized ring after applying one regime.

    x               : (N, 3) positions after the regime's shape map [fm].
    tangent         : (N, 3) unit tangents along the circulation.
    P               : (N, 3) polarization vectors.
    p_hat           : (N, 3) unit momentum directions (NaN where undefined).
    weight          : (N,) element weights.
    valid           : (N,) True where the momentum direction is defined.
    center          : (3,) centre of the undeformed ring [fm].
    center_deformed : (3,) image of the centre under the regime's shape map [fm].
    """
    regime: str
    x: np.ndarray
    tangent: np.ndarray
    P: np.ndarray
    p_hat: np.ndarray
    weight: np.ndarray
    valid: np.ndarray
    center: np.ndarray
    center_deformed: np.ndarray
    frame: JetFrame


def build_state(curve, flow, polarization, regime):
    """Apply one regime to a material curve and attach its polarization."""
    if regime not in REGIMES:
        raise ValueError(f"unknown regime '{regime}', expected one of {REGIMES}")
    if not 0.0 <= polarization.magnitude <= 1.0:
        raise ValueError(f"polarization magnitude must lie in [0, 1], got {polarization.magnitude}")
    if polarization.weight_mode not in ("material", "deformed_arclength"):
        raise ValueError(f"unknown weight_mode '{polarization.weight_mode}'")

    deform = regime in ("R1", "R3")
    deflect = regime in ("R2", "R3")

    if deform:
        deformation = shape_map_matrix(flow)
        x = apply_linear_map(deformation, curve.x)
        dx_dlam = apply_linear_map(deformation, curve.dx_dlam)
        center_deformed = deformation @ curve.center
    else:
        x = curve.x
        dx_dlam = curve.dx_dlam
        center_deformed = curve.center.copy()

    # R3 deflects from the deformed positions: the ring is first carried by
    # the flow and then emits along the flow at its new location.
    emission_matrix = hubble_matrix(flow) if deflect else np.eye(3)
    p_hat, valid = momentum_directions(emission_matrix, x)

    speed_undeformed = np.linalg.norm(curve.dx_dlam, axis=1)
    speed = np.linalg.norm(dx_dlam, axis=1)
    tangent = dx_dlam / speed[:, None]

    if polarization.stretch_scaling:
        magnitude = polarization.magnitude * speed / speed_undeformed
    else:
        magnitude = np.full(len(speed), polarization.magnitude)
    P = magnitude[:, None] * tangent

    if polarization.weight_mode == "material":
        weight = speed_undeformed * curve.dlam
    else:
        weight = speed * curve.dlam

    return RingState(regime=regime, x=x, tangent=tangent, P=P, p_hat=p_hat,
                     weight=weight, valid=valid, center=curve.center.copy(),
                     center_deformed=center_deformed, frame=curve.frame)
