# HelicityToyModel

> **Work in progress -- use with caution.**
> This code is still being tested and reviewed.
> Results should be treated as preliminary and cross-checked carefully
> before drawing any physics conclusions.

---

## Purpose

This folder contains a standalone Monte Carlo toy model designed to
isolate and diagnose two distinct detector-induced fake-polarization
effects that can appear in Lambda (and AntiLambda) hyperon measurements
at ALICE, specifically in the context of the **ring observable** analysis
described in PhysRevC.109.014905.

This follows Joseph Richard Adams' dissertation from 2021, titled
*"A Measurement of Lambda-Hyperon Spin Polarization in Au+Au Collisions
at sqrt(sNN) = 3 GeV with STAR"*. We thank Prof. Mike Lisa for the
suggestions and guidance that anticipated many of the problems of
fake-polarization signals due to efficiency effects.

The toy model is deliberately simple: it has no GEANT4 material budget,
no detector response, no smearing, and no beam-induced backgrounds.
Its value is precisely that simplicity: we aim to isolate each individual
effect in a controlled way that a full GEANT4 simulation cannot, because
in a full simulation all effects are active simultaneously and cannot
be switched on or off independently.

---

## Physics background

Crucially, both effects described below arise from correlations between
decay kinematics and detector- or analysis-induced efficiency variations.
They do not reflect any physical polarization of the Lambda, but instead
represent biases introduced by the measurement process.

### Effect 1: the helicity efficiency effect (forward-backward asymmetry)

When a minimum transverse-momentum threshold is applied to daughter
tracks, the two decay configurations are not treated equally.

- **Proton emitted forward** in the Lambda rest frame: the proton
  carries most of the lab-frame momentum (large pT in the lab), and the
  pion is soft. The soft pion is harder to reconstruct, so this
  configuration is **partially suppressed**.
- **Proton emitted backward** in the Lambda rest frame: the pion carries
  most of the lab-frame momentum and is well-reconstructed; the proton is
  softer. This configuration is **less affected** by the pT threshold.

The net result is a fake asymmetry in `cos(theta*)`, where `theta*` is
the angle between the proton emission direction in the Lambda rest frame
and the Lambda lab-frame momentum direction. This is the
**forward-backward asymmetry**, also known in the literature as the
**helicity efficiency effect**.

This effect primarily induces a distortion along the Lambda momentum
direction (a "longitudinal" component of the reconstructed polarization).
It does **not** create a fake ring observable by itself in a symmetric
setup, because the ring observable measures a *different* projection
(the left-right direction, see below). However, in the presence of
detector or acceptance asymmetries, or when combined with additional
cuts such as DCA selections, it can contribute indirectly to distortions
of the full decay angular distribution.

### Effect 2: the DCA-cut left-right asymmetry (phi* modulation)

When a minimum distance-of-closest-approach (DCA) cut is applied to
daughter tracks, the two daughters -- which carry opposite charges --
are bent in opposite directions by the solenoidal magnetic field.
Depending on the azimuthal angle `phi*` of the proton emission direction
around the Lambda axis, one of the daughters may be bent **toward** the
primary vertex (failing the DCA cut) while the other is bent away.

This creates, to leading order, a `sin(phi*)` modulation in the azimuthal
decay angle distribution -- a **left-right asymmetry** -- where
`phi* = 0` is the direction defined by `p_Lambda x z_hat` (the beam
direction). Crucially:

- The sign of this asymmetry is determined by the magnetic field
  direction **and** by the sign of the Lambda pseudorapidity.
- Flipping the magnetic field flips the asymmetry sign.
- Lambdas at `eta > 0` and Lambdas at `eta < 0` experience the asymmetry
  with **opposite signs**.

This effect originates from a combination of track curvature, DCA-based
selection, and reconstruction biases (e.g. vertex finding based on DCA
minimization), all of which depend on the relative geometry of the decay
and the detector.

Because the ring observable is constructed as a projection onto exactly
the `p_Lambda x z_hat` direction (when `z_hat` is used as the jet proxy),
this left-right asymmetry directly produces a **fake ring signal**.
Moreover, since the sign flips with eta, the fake signal is
**antisymmetric in eta** -- which is precisely the pattern observed in the
O-O data.

The traditional correction for this effect in ALICE global polarization
measurements is a Monte Carlo efficiency correction. However, if the
required correction is 100x larger than the physical signal, the
systematic uncertainty from that correction may dominate the measurement.
This toy model is intended to help design a **data-driven correction
strategy** that minimizes dependence on Monte Carlo efficiency modelling.

---

## Files in this folder

### `helicityEfficiencyToyModel.cxx`

**Run this file to generate data.**

The core Monte Carlo generator. Produces Lambdas from a Boltzmann mT
spectrum, propagates them to their decay vertex using the proper lifetime,
decays them isotropically via `TGenPhaseSpace` (no polarization), and
evaluates four cut scenarios:

| Scenario     | Applied cuts                          |
|--------------|---------------------------------------|
| `NoCuts`     | None (bug-check: must be flat)        |
| `pTCutOnly`  | Minimum pT on each daughter           |
| `DCACutOnly` | Minimum DCA_xy to PV on each daughter |
| `BothCuts`   | Both pT and DCA cuts simultaneously   |

Each scenario produces histograms separately for `eta > 0`, `eta < 0`,
and the combined sample. The key outputs are the `(cos theta*, phi*)` 2D
maps, the ring observable proxy distributions and profiles, and the
daughter kinematic distributions.

Output: one `.root` file per run, containing a structured directory tree.

**Run with ROOT:**

```bash
# Default (1M Lambdas, all default parameters):
root -l -b -q 'helicityEfficiencyToyModel.cxx'

# Custom parameters (see preamble in the file for full signature):
root -l -b -q 'helicityEfficiencyToyModel.cxx(2000000,"output.root",0.5,...)'
```

---

### `plotHelicityEfficiency.cxx`

**Run this file after the generator.**

Reads the ROOT file produced by `helicityEfficiencyToyModel.cxx` and
generates 10 multi-panel diagnostic figure PDFs. The most important
figures are:

- **Fig 1**: 2D `(cos theta*, phi*)` maps for each cut scenario -- the
  primary diagnostic for identifying which effect is active.
- **Fig 3**: `phi*` projections with a `sin(phi*)` fit to the DCA-cut
  distribution -- directly measures the left-right modulation amplitude.
- **Fig 4**: `<R_proxy>` vs Lambda pseudorapidity -- the **antisymmetric
  eta dependence** is the clearest signature of the DCA-cut fake signal.
- **Fig 10**: `(eta > 0) - (eta < 0)` difference plots -- isolates the
  antisymmetric fake component that appears as a non-zero ring observable
  in the data.

**Run with ROOT:**

```bash
root -l -b -q 'plotHelicityEfficiency.cxx("output.root","plots_dir")'
```

---

### `runHelicityToyModel.sh`

**The coordinator script -- start here for a full parameter scan.**

Runs the generator and plotter across seven families of parameter
variations (see below) and organises all output into a structured
directory tree under a configurable base directory. Supports parallel
execution (configurable via `MAX_PARALLEL`) and produces one log file
per run.

**Quick start:**

```bash
chmod +x runHelicityToyModel.sh

# Run all families (may take 30-60 minutes with MAX_PARALLEL=4):
./runHelicityToyModel.sh

# Run only one specific family:
./runHelicityToyModel.sh --family 1

# Preview all commands without executing:
./runHelicityToyModel.sh --dry-run
```

The scan families and their physics motivation are documented in detail
inside the script header. A brief summary:

| Family | Varied parameter        | Key question                                       |
|--------|-------------------------|----------------------------------------------------|
| 0      | Baseline                | Establish reference point                          |
| 1      | Magnetic field sign     | Does DCA asymmetry flip with B?                    |
| 2      | Daughter pT threshold   | How does helicity effect grow with cut?            |
| 3      | DCA threshold           | How does left-right asymmetry grow with cut?       |
| 4      | Lambda kinematic window | How large is the effect in the ring-analysis window? |
| 5      | Thermal temperature     | How does the spectrum shape affect the fake signal? |
| 6      | Statistics scaling      | Is 1/sqrt(N) scaling confirmed?                    |

**Default output directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/`

**Default log directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/logs/`

Both can be changed at the top of the script.

---

### `README.md`

This file.

---

## What this toy model cannot do

This is a **kinematic and geometric toy**, not a full detector simulation.
It deliberately excludes:

- Material budget and multiple scattering (which smear DCA distributions)
- TPC space-charge distortions (relevant in central Pb-Pb)
- ITS hit requirements and their eta/phi efficiency maps
- Momentum resolution smearing
- Secondary interactions and photon conversions
- Pile-up and beam-induced backgrounds
- Any polarization physics (the decay is always isotropic)
- Continuous efficiency variations with pT, eta, and detector occupancy
  (the model effectively applies hard selection thresholds)

The toy model gives **qualitative and mechanistic** insight: it tells
you which effect is responsible for which asymmetry, with which sign,
and how it scales with the cut parameters. It does **not** give
quantitative correction factors for the data. Those require a full
GEANT4-level simulation with the ALICE detector geometry and the
same reconstruction algorithm used on data.

---

## Known limitations and open issues

- [ ] The DCA calculation uses an ideal analytical helix with a perfectly
      uniform axial magnetic field. Real ALICE tracks are affected by
      space-charge distortions and imperfect field maps.
- [ ] The pT spectrum uses a simple Boltzmann mT exponential. A measured
      Lambda pT spectrum (e.g., from ALICE Pb-Pb or O-O) could be used
      as a more realistic input via a TH1-based sampler.
- [ ] Only `Lambda -> p + pi-` is simulated. `AntiLambda -> pbar + pi+`
      is charge-conjugate and should behave identically by symmetry,
      but an explicit AntiLambda simulation would confirm this.
- [ ] The ring observable proxy uses the beam direction (+z) as the jet
      axis. A more realistic proxy would sample jet directions from a
      measured or simulated jet spectrum with realistic eta/phi coverage.
- [ ] No feed-down Lambdas from Xi or Sigma0 decays are simulated.
      Feed-down adds displaced secondary Lambdas whose decay geometry
      is qualitatively different (larger decay radius, different DCA
      distributions).
- [ ] The plotting macro has not been tested on small-statistics runs
      (N < 50k) where some histogram bins may be empty and the sin-fit
      in Fig 3 may fail to converge.
- [ ] The code has not been compiled and tested with ROOT 6.28+; only
      static analysis (brace/parenthesis balance, ASCII cleanliness)
      has been performed. **The code must be run and its output
      carefully inspected before any physics conclusions are drawn.**