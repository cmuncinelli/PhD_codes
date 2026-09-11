"""
observable.py -- Ring observable, angular variables, and closed-form references.

Part of the ring geometry toy (README.md holds definitions and derivations).

For each ring element the observable uses only its momentum direction and
its polarization vector: n_hat = t_hat x p_hat / |t_hat x p_hat| and
r = P . n_hat. No hyperons or decays are simulated, so every average below is
an exact weighted mean over ring elements.
"""

from dataclasses import dataclass

import numpy as np

from .geometry import DEGENERACY_EPS, Z_HAT


def wrap_to_pi(angle):
    """Wrap angles into [-pi, pi)."""
    return np.mod(angle + np.pi, 2.0 * np.pi) - np.pi


def pseudorapidity(p_hat):
    """Pseudorapidity of unit direction(s): eta = atanh(cos theta)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.arctanh(p_hat[..., 2])


def azimuth(p_hat):
    """Azimuth of direction(s) in [-pi, pi]."""
    return np.arctan2(p_hat[..., 1], p_hat[..., 0])


def ring_directions(t_hat, p_hat):
    """Ring basis about the jet for each momentum direction.

    Returns (n_hat, theta_hat, delta_theta, valid):
      n_hat       : t_hat x p_hat normalized (ring direction).
      theta_hat   : n_hat x p_hat, the direction of increasing opening angle;
                    (p_hat, theta_hat, n_hat) is right-handed.
      delta_theta : opening angle between jet and momentum, Delta theta_{Jet,Lambda}.
      valid       : False where p_hat is (anti)collinear with the jet.
    """
    cross = np.cross(t_hat, p_hat)
    sin_angle = np.linalg.norm(cross, axis=-1)
    cos_angle = p_hat @ t_hat
    # atan2 keeps full precision near 0 and pi, where arccos does not.
    delta_theta = np.arctan2(sin_angle, cos_angle)
    valid = sin_angle > DEGENERACY_EPS
    n_hat = np.full_like(cross, np.nan)
    n_hat[valid] = cross[valid] / sin_angle[valid, None]
    theta_hat = np.cross(n_hat, p_hat)
    return n_hat, theta_hat, delta_theta, valid


@dataclass(frozen=True)
class Acceptance:
    """Kinematic window; eta stands in for rapidity (directions only)."""
    eta_max: float = None


@dataclass(frozen=True)
class RingMeasurement:
    """Per-element quantities of the ring observable.

    r, theta_component, helicity_component : P projected on n_hat,
        theta_hat and p_hat.
    cos_gamma   : cosine of the angle between P and n_hat (NaN where P = 0).
    delta_phi   : phi_p - phi_jet wrapped into [-pi, pi).
    delta_theta : opening angle Delta theta_{Jet,Lambda}.
    accepted    : element has a defined ring direction and passes the window.
    """
    r: np.ndarray
    theta_component: np.ndarray
    helicity_component: np.ndarray
    cos_gamma: np.ndarray
    delta_phi: np.ndarray
    delta_theta: np.ndarray
    eta: np.ndarray
    weight: np.ndarray
    accepted: np.ndarray
    n_undefined: int


def measure(state, acceptance=Acceptance()):
    """Evaluate the ring observable element by element."""
    frame = state.frame
    n_hat, theta_hat, delta_theta, valid_ring = ring_directions(frame.t_hat, state.p_hat)
    defined = state.valid & valid_ring

    r = np.einsum("ij,ij->i", state.P, n_hat)
    theta_component = np.einsum("ij,ij->i", state.P, theta_hat)
    helicity_component = np.einsum("ij,ij->i", state.P, state.p_hat)

    p_norm = np.linalg.norm(state.P, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cos_gamma = np.where(p_norm > 0.0, r / p_norm, np.nan)

    eta = pseudorapidity(state.p_hat)
    delta_phi = wrap_to_pi(azimuth(state.p_hat) - frame.phi)

    accepted = defined.copy()
    if acceptance.eta_max is not None:
        with np.errstate(invalid="ignore"):
            accepted &= np.abs(eta) < acceptance.eta_max

    return RingMeasurement(r=r, theta_component=theta_component,
                           helicity_component=helicity_component,
                           cos_gamma=cos_gamma, delta_phi=delta_phi,
                           delta_theta=delta_theta, eta=eta,
                           weight=state.weight, accepted=accepted,
                           n_undefined=int(np.count_nonzero(~defined)))


@dataclass(frozen=True)
class IntegratedSums:
    """Additive sums behind a RingSummary.

    Kept separate from the summary so that results from several rings or jet
    orientations can be added before any ratio is formed.
    """
    sum_w: float
    sum_wr: float
    sum_w_theta: float
    sum_w_helicity: float
    total_weight: float
    n_undefined: int

    def __add__(self, other):
        return IntegratedSums(sum_w=self.sum_w + other.sum_w,
                              sum_wr=self.sum_wr + other.sum_wr,
                              sum_w_theta=self.sum_w_theta + other.sum_w_theta,
                              sum_w_helicity=self.sum_w_helicity + other.sum_w_helicity,
                              total_weight=self.total_weight + other.total_weight,
                              n_undefined=self.n_undefined + other.n_undefined)


@dataclass(frozen=True)
class RingSummary:
    """Integrated results.

    ring_average            : <P . n_hat> over accepted elements.
    theta_component         : <P . theta_hat> over accepted elements.
    helicity_component      : <P . p_hat> over accepted elements.
    acceptance_fraction     : accepted weight / total weight.
    signal_per_total_weight : sum(w r) over accepted / total weight; the value
                              a flat unpolarized background would dilute to,
                              up to its normalization.
    n_undefined             : elements without a defined ring direction.
    """
    ring_average: float
    theta_component: float
    helicity_component: float
    acceptance_fraction: float
    signal_per_total_weight: float
    n_undefined: int


def integrate(measurement):
    """Additive integrated sums over the accepted elements of a measurement."""
    sel = measurement.accepted
    w = measurement.weight[sel]
    return IntegratedSums(sum_w=float(np.sum(w)),
                          sum_wr=float(np.sum(w * measurement.r[sel])),
                          sum_w_theta=float(np.sum(w * measurement.theta_component[sel])),
                          sum_w_helicity=float(np.sum(w * measurement.helicity_component[sel])),
                          total_weight=float(np.sum(measurement.weight)),
                          n_undefined=measurement.n_undefined)


def summary_from_sums(sums):
    """Form the ratios of a RingSummary from additive sums."""
    if sums.sum_w > 0.0:
        ring_average = sums.sum_wr / sums.sum_w
        theta_component = sums.sum_w_theta / sums.sum_w
        helicity_component = sums.sum_w_helicity / sums.sum_w
    else:
        ring_average = theta_component = helicity_component = np.nan
    return RingSummary(ring_average=float(ring_average),
                       theta_component=float(theta_component),
                       helicity_component=float(helicity_component),
                       acceptance_fraction=float(sums.sum_w / sums.total_weight),
                       signal_per_total_weight=float(sums.sum_wr / sums.total_weight),
                       n_undefined=sums.n_undefined)


def summarize(measurement):
    """Integrate a measurement into a RingSummary."""
    return summary_from_sums(integrate(measurement))


@dataclass(frozen=True)
class BinnedRing:
    """Per-bin sums of weights and weighted ring values.

    Sums, not ratios, are stored so that histograms can be added before any
    ratio is formed.
    """
    edges: np.ndarray
    sum_w: np.ndarray
    sum_wr: np.ndarray
    total_weight: float

    @property
    def centers(self):
        return 0.5 * (self.edges[1:] + self.edges[:-1])

    @property
    def widths(self):
        return np.diff(self.edges)

    def __add__(self, other):
        if not np.array_equal(self.edges, other.edges):
            raise ValueError("cannot add BinnedRing objects with different edges")
        return BinnedRing(edges=self.edges, sum_w=self.sum_w + other.sum_w,
                          sum_wr=self.sum_wr + other.sum_wr,
                          total_weight=self.total_weight + other.total_weight)

    def ring_average(self):
        """Mean of P . n_hat over the ring elements in each bin (NaN if empty)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.sum_w > 0.0, self.sum_wr / self.sum_w, np.nan)

    def signal_density(self):
        """sum(w r) per unit of the binned variable, per total ring weight.

        With a background of unpolarized emitters that is flat in the binned
        variable and dominates the yield, the measured ring value in each bin
        is proportional to this density; it is the counterpart of the
        R(Delta phi) curves of the model papers.
        """
        return self.sum_wr / (self.total_weight * self.widths)


def bin_measurement(values, measurement, edges):
    """Histogram accepted elements of a measurement in a per-element variable."""
    sel = measurement.accepted
    edges = np.asarray(edges, dtype=np.float64)
    sum_w, _ = np.histogram(values[sel], bins=edges, weights=measurement.weight[sel])
    sum_wr, _ = np.histogram(values[sel], bins=edges,
                             weights=measurement.weight[sel] * measurement.r[sel])
    return BinnedRing(edges=edges, sum_w=sum_w, sum_wr=sum_wr,
                      total_weight=float(np.sum(measurement.weight)))


# ---------------------------------------------------------------------------
# Closed-form references from the Ring Design Discussion note (used by tests)
# ---------------------------------------------------------------------------

def beam_referenced_basis(p_hat):
    """Beam-referenced (phi_hat, theta_hat) at direction(s) p_hat."""
    cross = np.cross(Z_HAT, p_hat)
    phi_hat = cross / np.linalg.norm(cross, axis=-1)[..., None]
    theta_hat = np.cross(phi_hat, p_hat)
    return phi_hat, theta_hat


def cos_sin_c_closed_form(eta_jet, eta_lambda, delta_phi):
    """cos C = n_hat . phi_hat and sin C = n_hat . theta_hat in closed form."""
    big_n = np.sinh(eta_jet) - np.sinh(eta_lambda) * np.cos(delta_phi)
    big_m = np.cosh(eta_lambda) * np.sin(delta_phi)
    norm = np.hypot(big_n, big_m)
    return big_n / norm, -big_m / norm


def elliptic_average_cos_c(eta_jet):
    """Uniform-Delta-phi average of cos C at eta_lambda = 0.

    (2/pi) tanh(eta_jet) K(k), k = sech(eta_jet). scipy's ellipk takes the
    parameter m = k^2, not the modulus k.
    """
    from scipy.special import ellipk
    if eta_jet == 0.0:
        # The formula is 0 * infinity at the origin; the limit is 0.
        return 0.0
    return 2.0 / np.pi * np.tanh(eta_jet) * ellipk(1.0 / np.cosh(eta_jet) ** 2)
