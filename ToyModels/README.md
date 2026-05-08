# ToyModels

> **Work in progress -- use with caution.**
> This code is still being tested and reviewed. Particularly the helicityEfficiencyToyModel and the AEE numerical integration parts!
> Results should be treated as preliminary and cross-checked carefully before drawing any physics conclusions.

---

## Purpose

This folder contains a collection of standalone toy models and numerical integration scripts designed to isolate and diagnose distinct detector-induced fake-polarization effects that can appear in Lambda (and AntiLambda) hyperon measurements at ALICE. These are specifically tailored to explore the context of the **ring observable** analysis described in PhysRevC.109.014905.

The models herein are deliberately simple: they have no GEANT4 material budget, no detailed detector response, no momentum smearing, and no beam-induced backgrounds. Their value lies precisely in this simplicity: we aim to isolate each individual kinematic or geometric effect in a controlled way that a full GEANT4 simulation cannot, because in a full simulation all effects are active simultaneously and cannot be easily switched on or off independently.

The physical intuition is based on Joseph Richard Adams' dissertation from 2021, titled *"A Measurement of Lambda-Hyperon Spin Polarization in Au+Au Collisions at sqrt(sNN) = 3 GeV with STAR"*. For that, we thank the author and his advisor, Prof. Mike Lisa.

---

## Physics background

Crucially, the effects described below arise from correlations between decay kinematics and detector- or analysis-induced efficiency variations. They do not reflect any physical polarization of the Lambda, but instead represent biases introduced by the measurement process.

### Effect 1: The helicity efficiency effect (forward-backward asymmetry)

When a minimum transverse-momentum threshold is applied to daughter tracks, the two decay configurations are not treated equally.

- **Proton emitted forward** in the Lambda rest frame: the proton carries most of the lab-frame momentum (large pT in the lab), and the pion is soft. The soft pion is harder to reconstruct, so this configuration is **partially suppressed**.
- **Proton emitted backward** in the Lambda rest frame: the pion carries most of the lab-frame momentum and is well-reconstructed; the proton is softer. This configuration is **less affected** by the pT threshold.

The net result is a fake asymmetry in `cos(theta*)`, where `theta*` is the angle between the proton emission direction in the Lambda rest frame and the Lambda lab-frame momentum direction. This is the **forward-backward asymmetry**, also known in the literature as the **helicity efficiency effect** (HEE).

This effect primarily induces a distortion along the Lambda momentum direction. It should **not** create a fake ring observable by itself, as the observable kills all polarization components colinear to the Lambda momentum, but in the presence of detector or acceptance asymmetries, it is possible that HEE contributes indirectly to distortions of the full decay angular distribution. Still needs to be checked!

### Effect 2: The DCA-cut left-right asymmetry

When a minimum distance-of-closest-approach (DCA) cut is applied to daughter tracks, the two daughters -- which carry opposite charges -- are bent in opposite directions by the solenoidal magnetic field. Depending on the azimuthal angle `phi*` of the proton emission direction around the Lambda axis, the DCA between daughters may be larger or smaller, creating a prioritized geometry during reconstruction and thus a fake polarization signal.

*Note on the Toy Model framework:* This Toy Model explicitly lacks a simulated DCA between the daughters themselves (i.e., the daughters originate exactly from the same spatial point). Therefore, any left-right asymmetry observed here confirms that the effect is driven predominantly by the **DCA to the Primary Vertex (PV)**, rather than the DCA between the daughters!

### Effect 3: The Azimuthal Emission Efficiency (AEE)

The Azimuthal Emission Efficiency (AEE) creates an artificial structure in the angular distributions due to the non-uniform acceptance of daughters in the laboratory frame, strongly correlating the Lambda kinematics with the jet axis. The DCA between daughters may be larger or smaller, creating a prioritized geometry during reconstruction and thus a fake polarization signal. If not driven by DCA between the daughters, it can be driven by DCA
of the daughters to the PV. This is one of the main things that the Toy Model
is testing for!

Because the ring observable is constructed as a projection onto the `p_Lambda x t_hat` direction (where t_hat is the jet direction), this left-right asymmetry directly produces a **fake ring signal**.

When evaluating the Ring Observable, the AEE induces a spurious component that can be mathematically formulated as:
`R^{AEE}_Lambda(theta, theta_Lambda, Delta_phi) = [ -sin(theta) cos(theta_Lambda) cos(Delta_phi) + cos(theta) sin(theta_Lambda) ] / sqrt( 1 - (sin(theta) sin(theta_Lambda) cos(Delta_phi) + cos(theta) cos(theta_Lambda))^2 )`

Numerical integration reveals that kinematic biases (such as the distributions of the jet correlation `Delta_phi` and the Lambda pseudorapidity `eta_Lambda`) produce visible distortions matching experimental data. However, full integration over `dOmega_Jet` and `dOmega*_Lambda` theoretically cancels this AEE-induced effect out, leading to a pure vortex-polarization observable.

---

## Files in this folder

### `numericalAEEIntegration.cxx`

**Run this script to evaluate AEE effects.**

This script performs the numerical integration of the Ring Observable component induced by the AEE effect. It compares an idealized flat phase space baseline against a realistic parameterization derived from experimental data weights (accounting for jet correlation `Delta_phi` and `eta_Lambda` distributions). 

The script evaluates this observable as a function of the trigger (jet) pseudorapidity (`eta_t`) and calculates the marginal distribution as a function of `eta_Lambda` to assess detector-induced kinematic biases. It serves as a proof-of-concept that full phase space integration cancels the spurious AEE effect.

**Run with ROOT:**
```
root -l numericalAEEIntegration.cxx
```

### `helicityEfficiencyToyModel.cxx`

**Run this file to generate data for Effects 1 & 2**.

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

Each run produces a ROOT file containing two parallel histogram directories to isolate geometric acceptance boundaries:
- `WithoutEtaGate/` -- legacy, no daughter eta requirement (kept for reference)
- `WithEtaGate/` -- physically consistent set (BOTH daughters in acceptance)

Each family contains the four scenarios mentioned above, further split into `EtaPos`, `EtaNeg`, and `All` sub-directories. The key outputs are the `(cos theta*, phi*)` 2D maps, the ring observable proxy distributions and profiles, and the daughter kinematic distributions.

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

Runs the generator and plotter across ten families of parameter
variations specifically designed to isolate one effect by silencing the other.
It organizes all output into a structured directory tree under a configurable
base directory. Supports parallel execution (default `MAX_PARALLEL=60`) and
produces one log file per run.

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

| Family | Varied parameter        | Key question                                                                 |
|--------|-------------------------|------------------------------------------------------------------------------|
| 0      | Baseline                | Establish flat reference point with no cuts.                                 |
| 1      | Asymmetric DCA cuts     | [AEE] Does the sign depend on which daughter's DCA dominates?                |
| 2      | Symmetric DCA cuts      | [AEE] How does magnitude grow as a clean function of cut strength?           |
| 3      | Magnetic field          | [AEE] Does asymmetry flip with B field polarity and scale with strength?     |
| 4      | Lambda pT minimum       | How does the combined fake signal depend on the Lambda pT regime?            |
| 5      | Daughter pT cuts        | [HEE] How does HEE isolate and grow with symmetric/asymmetric pT cuts?       |
| 6      | Eta acceptance window   | How does fake signal magnitude change with detector acceptance?              |
| 7      | Temperature scan        | How does the Lambda pT spectrum shape (Boltzmann T) affect the fake signal?  |
| 8      | Ring kinematic windows  | Where is the fake signal largest across successive pT windows?               |
| 9      | Realistic ALICE cuts    | What is the combined HEE+AEE estimate under experimental conditions?         |

**Default output directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/`
**Default log directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/logs/`

Both can be changed at the top of the script. The entire run generates approximately 600 MB - 1.4 GB of data between compressed ROOT files and PDFs.

---

## What these toy models cannot do

This is a **kinematic and geometric tool**, not a full detector simulation.
It deliberately excludes:

- Material budget and multiple scattering (which smear DCA distributions)
- TPC space-charge distortions (relevant in central Pb-Pb or O-O)
- ITS hit requirements and their eta/phi efficiency maps
- Momentum resolution smearing and secondary interactions
- Secondary interactions and photon conversions
- Pile-up and beam-induced backgrounds
- Any real polarization physics (decays are isotropic or flat baselines).
- Continuous efficiency variations with pT, eta, and detector occupancy
  (the `helicityEfficiencyToyModel.cxx` model only applies hard selection thresholds)

The toy models give **qualitative and mechanistic** insights ONLY: they tell
you which effect is responsible for which asymmetry, with which sign,
and how it scales with the cut parameters. They **do not** give
quantitative correction factors for the data. Those require a full
GEANT4-level simulation with the ALICE detector geometry and the
same reconstruction algorithm used on data.

---

## Known limitations and open issues for `helicityEfficiencyToyModel.cxx`

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