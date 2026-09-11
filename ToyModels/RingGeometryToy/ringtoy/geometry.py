"""
geometry.py -- Jet frame, ring-curve families, and ring placement.

Part of the ring geometry toy (README.md holds conventions and physics).
All vectors are lab-frame Cartesian float64 arrays; lengths are in fm.
"""

from dataclasses import dataclass

import numpy as np

# Floor used to detect degenerate directions (zero-length projections and
# cross products). Rings are fm-sized, so any genuine geometric quantity is
# many orders of magnitude above it, while double round-off stays below it.
DEGENERACY_EPS = 1.0e-12

X_HAT = np.array([1.0, 0.0, 0.0])
Z_HAT = np.array([0.0, 0.0, 1.0])


def direction_from_eta_phi(eta, phi):
    """Unit vector(s) with pseudorapidity eta and azimuth phi (broadcasts)."""
    eta = np.asarray(eta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    # cos(theta) = tanh(eta) and sin(theta) = 1/cosh(eta) are used directly:
    # going through theta = 2*atan(exp(-eta)) loses precision at large |eta|.
    sin_theta = 1.0 / np.cosh(eta)
    components = np.broadcast_arrays(sin_theta * np.cos(phi),
                                     sin_theta * np.sin(phi),
                                     np.tanh(eta))
    return np.stack(components, axis=-1)


@dataclass(frozen=True)
class JetFrame:
    """Right-handed basis (e1, e2, t_hat) attached to the jet axis.

    e2 is the projection of the beam axis onto the plane transverse to the
    jet ("up" in jet-plane plots) and e1 = e2 x t_hat ("right"). With this
    choice a jet pointing at the viewer sees right-handed circulation as
    counterclockwise in the (e1, e2) plane.
    """
    eta: float
    phi: float
    t_hat: np.ndarray
    e1: np.ndarray
    e2: np.ndarray


def make_jet_frame(eta, phi):
    """Build the jet frame for a jet with pseudorapidity eta and azimuth phi."""
    t_hat = direction_from_eta_phi(eta, phi)
    up = Z_HAT - np.dot(Z_HAT, t_hat) * t_hat
    norm_up = np.linalg.norm(up)
    if norm_up < DEGENERACY_EPS:
        # Jet along the beam: the beam has no transverse projection, so the
        # "up" axis falls back to the projection of x_hat.
        up = X_HAT - np.dot(X_HAT, t_hat) * t_hat
        norm_up = np.linalg.norm(up)
    e2 = up / norm_up
    e1 = np.cross(e2, t_hat)
    return JetFrame(eta=float(eta), phi=float(phi), t_hat=t_hat, e1=e1, e2=e2)


@dataclass(frozen=True)
class RingShape:
    """Shape of the ring core curve, defined in the jet frame.

    kind            : "circle", "ellipse" or "fourier".
    radius          : circle radius, ellipse semi-axis along e1 (before
                      ellipse_angle), or mean radius of the Fourier ring [fm].
    radius_minor    : ellipse semi-axis along e2 (before ellipse_angle) [fm].
    ellipse_angle   : in-plane rotation of the ellipse axes [rad].
    fourier_modes   : ((m, amplitude, phase), ...), relative modulation
                      rho = radius * (1 + sum amplitude*cos(m*(lambda - phase))).
    tilt            : rotation of the ring plane out of the plane transverse
                      to the jet [rad].
    tilt_axis_angle : in-plane direction of the tilt axis, from e1 [rad].
    circulation     : +1 for right-handed circulation about t_hat, -1 reversed.
    n_points        : number of material points on the curve.
    """
    kind: str = "circle"
    radius: float = 1.0
    radius_minor: float = 1.0
    ellipse_angle: float = 0.0
    fourier_modes: tuple = ()
    tilt: float = 0.0
    tilt_axis_angle: float = 0.0
    circulation: int = +1
    n_points: int = 720


@dataclass(frozen=True)
class RingPlacement:
    """Position of the ring centre.

    distance_along_jet : L [fm]; L > 0 places the ring outward along the jet,
                         L < 0 inward (behind the vertex).
    offset_e1/offset_e2: transverse offsets of the centre from the jet line [fm].
    center_lab         : explicit lab-frame centre [fm]; when given it
                         overrides the three fields above (used for jet--flow
                         misalignment, where the centre is fixed in the lab
                         while the jet direction varies).
    """
    distance_along_jet: float = 0.0
    offset_e1: float = 0.0
    offset_e2: float = 0.0
    center_lab: tuple = None


@dataclass(frozen=True)
class MaterialCurve:
    """Material points of the undeformed ring and their analytic tangents.

    lam     : curve parameter on a uniform periodic grid [rad].
    x       : (N, 3) lab positions [fm].
    dx_dlam : (N, 3) derivative of x with respect to lam [fm/rad], pointing
              along the circulation.
    dlam    : grid spacing of lam [rad].
    center  : (3,) ring centre [fm].
    frame   : jet frame used to build the ring.
    """
    lam: np.ndarray
    x: np.ndarray
    dx_dlam: np.ndarray
    dlam: float
    center: np.ndarray
    frame: JetFrame


def ring_center(placement, frame):
    """Lab-frame centre of the ring for a given placement and jet frame."""
    if placement.center_lab is not None:
        return np.asarray(placement.center_lab, dtype=np.float64).reshape(3)
    return (placement.distance_along_jet * frame.t_hat
            + placement.offset_e1 * frame.e1
            + placement.offset_e2 * frame.e2)


def _rotation_matrix(axis, angle):
    """Rodrigues rotation matrix about a (not necessarily unit) axis."""
    k = np.asarray(axis, dtype=np.float64)
    k = k / np.linalg.norm(k)
    k_cross = np.array([[0.0, -k[2], k[1]],
                        [k[2], 0.0, -k[0]],
                        [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(angle) * k_cross + (1.0 - np.cos(angle)) * (k_cross @ k_cross)


def _require_positive(value, name):
    if not value > 0.0:
        raise ValueError(f"{name} must be positive, got {value}")


def build_ring(shape, placement, frame):
    """Material points and tangents of the ring core curve in the lab frame."""
    n = int(shape.n_points)
    if n < 8:
        raise ValueError(f"n_points must be at least 8, got {n}")
    if shape.circulation not in (+1, -1):
        raise ValueError(f"circulation must be +1 or -1, got {shape.circulation}")

    # Uniform grid in the curve parameter. The curves are smooth and periodic,
    # so the trapezoid rule on this grid is spectrally accurate for every ring
    # average; adaptive quadrature would add cost without adding accuracy.
    dlam = 2.0 * np.pi / n
    lam = dlam * np.arange(n)
    cos_l = np.cos(lam)
    sin_l = np.sin(lam)

    if shape.kind == "circle":
        _require_positive(shape.radius, "radius")
        u = shape.radius * cos_l
        v = shape.radius * sin_l
        du = -shape.radius * sin_l
        dv = shape.radius * cos_l
    elif shape.kind == "ellipse":
        _require_positive(shape.radius, "radius")
        _require_positive(shape.radius_minor, "radius_minor")
        u0 = shape.radius * cos_l
        v0 = shape.radius_minor * sin_l
        du0 = -shape.radius * sin_l
        dv0 = shape.radius_minor * cos_l
        c_a = np.cos(shape.ellipse_angle)
        s_a = np.sin(shape.ellipse_angle)
        u = c_a * u0 - s_a * v0
        v = s_a * u0 + c_a * v0
        du = c_a * du0 - s_a * dv0
        dv = s_a * du0 + c_a * dv0
    elif shape.kind == "fourier":
        _require_positive(shape.radius, "radius")
        modulation = np.ones(n)
        dmodulation = np.zeros(n)
        for m, amplitude, phase in shape.fourier_modes:
            # Non-integer m would break periodicity and silently open the curve.
            if not float(m).is_integer():
                raise ValueError(f"Fourier mode numbers must be integers, got {m}")
            argument = m * (lam - phase)
            modulation += amplitude * np.cos(argument)
            dmodulation -= amplitude * m * np.sin(argument)
        rho = shape.radius * modulation
        if np.min(rho) <= 0.0:
            raise ValueError("Fourier modulation makes the radius non-positive")
        drho = shape.radius * dmodulation
        u = rho * cos_l
        v = rho * sin_l
        du = drho * cos_l - rho * sin_l
        dv = drho * sin_l + rho * cos_l
    else:
        raise ValueError(f"unknown ring kind '{shape.kind}'")

    zeros = np.zeros(n)
    local = np.stack([u, v, zeros], axis=1)
    dlocal = np.stack([du, dv, zeros], axis=1)

    if shape.tilt != 0.0:
        axis = np.array([np.cos(shape.tilt_axis_angle), np.sin(shape.tilt_axis_angle), 0.0])
        rotation = _rotation_matrix(axis, shape.tilt)
        local = local @ rotation.T
        dlocal = dlocal @ rotation.T

    # Reversing the circulation reverses only the traversal direction; the
    # material points themselves stay where they are.
    dlocal = shape.circulation * dlocal

    # Rows of the basis are (e1, e2, t_hat): local (u, v, w) --> u*e1 + v*e2 + w*t_hat.
    basis = np.stack([frame.e1, frame.e2, frame.t_hat], axis=0)
    center = ring_center(placement, frame)
    x = center + local @ basis
    dx_dlam = dlocal @ basis

    return MaterialCurve(lam=lam, x=x, dx_dlam=dx_dlam, dlam=dlam, center=center, frame=frame)
