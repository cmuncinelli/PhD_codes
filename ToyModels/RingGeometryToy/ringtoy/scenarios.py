"""
scenarios.py -- Toy configuration, named presets, and command-line overrides.

Part of the ring geometry toy (README.md describes what each preset studies).
"""

import ast
import dataclasses
from dataclasses import dataclass, field

import numpy as np

from .almond import AlmondModel
from .flow import EllipticFlow
from .geometry import RingPlacement, RingShape
from .observable import Acceptance
from .regimes import PolarizationModel


@dataclass(frozen=True)
class Binning:
    """Number of bins of the differential observables."""
    n_delta_phi: int = 72
    n_delta_theta: int = 60


@dataclass(frozen=True)
class ScanSpec:
    """Parameter scan attached to a preset.

    kind           : "none", "parameter", "position" or "misaligned".
    parameter      : dotted configuration key for kind="parameter", or the
                     special name "delta" (jet azimuth relative to psi2).
    values         : scanned values.
    n_jet_azimuths : jet orientations averaged over for kind="misaligned".
    """
    kind: str = "none"
    parameter: str = ""
    values: tuple = ()
    n_jet_azimuths: int = 180


@dataclass(frozen=True)
class ToyConfig:
    """Complete configuration of one toy run."""
    jet_eta: float = 0.0
    jet_phi: float = 0.0
    shape: RingShape = field(default_factory=RingShape)
    placement: RingPlacement = field(default_factory=RingPlacement)
    flow: EllipticFlow = field(default_factory=EllipticFlow)
    polarization: PolarizationModel = field(default_factory=PolarizationModel)
    acceptance: Acceptance = field(default_factory=Acceptance)
    binning: Binning = field(default_factory=Binning)
    scan: ScanSpec = field(default_factory=ScanSpec)
    almond: AlmondModel = field(default_factory=AlmondModel)


def _base():
    # Generic, non-special defaults: psi2 is deliberately not aligned with the
    # jet, so that R1, R2 and R3 are all distinct from R0.
    return ToyConfig(
        jet_eta=0.0, jet_phi=0.0,
        shape=RingShape(kind="circle", radius=1.0, n_points=720),
        placement=RingPlacement(distance_along_jet=2.0),
        flow=EllipticFlow(v2=0.2, psi2=0.3, strength=1.0, h_z=1.0),
        polarization=PolarizationModel(magnitude=1.0),
        acceptance=Acceptance(eta_max=None),
    )


def _preset_regimes():
    return _base()


def _preset_position():
    # Analogue of the insertion-position scan of the model papers: positive
    # distances are outward along the jet, negative ones inward.
    return dataclasses.replace(
        _base(), acceptance=Acceptance(eta_max=0.5),
        scan=ScanSpec(kind="position", values=(-4.0, -2.0, 0.0, 2.0, 4.0)))


def _preset_misaligned():
    # Ring centre fixed in the lab at x = 2 fm while the jet azimuth is
    # scanned uniformly, as in the jet-alignment study of the model papers.
    return dataclasses.replace(
        _base(), placement=RingPlacement(center_lab=(2.0, 0.0, 0.0)),
        flow=EllipticFlow(v2=0.2, psi2=0.0, strength=1.0, h_z=1.0),
        acceptance=Acceptance(eta_max=0.5),
        scan=ScanSpec(kind="misaligned", n_jet_azimuths=180))


def _preset_eta_jet():
    return dataclasses.replace(
        _base(), acceptance=Acceptance(eta_max=0.5),
        scan=ScanSpec(kind="parameter", parameter="jet_eta",
                      values=tuple(np.linspace(-1.5, 1.5, 31))))


def _preset_delta():
    # Elliptic flow has period pi in the jet angle relative to psi2.
    return dataclasses.replace(
        _base(), jet_eta=0.3, acceptance=Acceptance(eta_max=0.8),
        scan=ScanSpec(kind="parameter", parameter="delta",
                      values=tuple(np.linspace(0.0, np.pi, 25))))


def _preset_strength():
    return dataclasses.replace(
        _base(),
        scan=ScanSpec(kind="parameter", parameter="flow.strength",
                      values=tuple(np.linspace(0.0, 2.0, 21))))


PRESETS = {
    "regimes": _preset_regimes,
    "position": _preset_position,
    "misaligned": _preset_misaligned,
    "eta_jet": _preset_eta_jet,
    "delta": _preset_delta,
    "strength": _preset_strength,
}


def get_preset(name):
    """Configuration of a named preset."""
    if name not in PRESETS:
        raise ValueError(f"unknown preset '{name}', expected one of {sorted(PRESETS)}")
    return PRESETS[name]()


def apply_override(config, dotted_key, value):
    """Return a copy of a (nested, frozen) configuration with one field replaced."""
    head, _, rest = dotted_key.partition(".")
    if not dataclasses.is_dataclass(config) or head not in {f.name for f in dataclasses.fields(config)}:
        raise KeyError(f"unknown configuration key '{dotted_key}'")
    if rest:
        value = apply_override(getattr(config, head), rest, value)
    return dataclasses.replace(config, **{head: value})


def parse_override(text):
    """Parse 'key=value'; the value is a Python literal, or a bare string."""
    key, sep, raw = text.partition("=")
    if not sep:
        raise ValueError(f"override '{text}' is not of the form key=value")
    try:
        value = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        # Bare words such as kind=ellipse are taken as strings.
        value = raw
    return key.strip(), value
