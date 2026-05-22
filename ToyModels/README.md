# ToyModels

> **Work in progress -- use with caution.**
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
evaluates four cut scenarios. It mostly utilizes `ROOT::Math::PxPyPzEVector` and `ROOT::Math::XYZVector` for optimized kinematic calculations (even though TLorentzVector was kept to use the -- very -- convenient TGenPhaseSpace class), and features a native `main()` entrypoint allowing it to be compiled directly into a standalone C++ binary for performance gains.

| Scenario     | Applied cuts                          |
|--------------|---------------------------------------|
| `NoCuts`     | None (bug-check: must be flat)        |
| `pTCutOnly`  | Minimum pT on each daughter           |
| `DCACutOnly` | Minimum DCA_xy to PV on each daughter |
| `BothCuts`   | Both pT and DCA cuts simultaneously   |

Each run produces a ROOT file containing two parallel histogram directories to isolate geometric acceptance boundaries:
- `WithoutEtaGate/` -- legacy, no daughter eta requirement (kept for reference and for comparison on what type of biases are introduced by the eta gate)
- `WithEtaGate/` -- physically consistent set (BOTH daughters in acceptance)

Each family contains the four scenarios mentioned above, further split into
`EtaPos`, `EtaNeg`, and `All` sub-directories.

#### Key output histograms

The main outputs per scenario directory are:

**Emission angle distributions:**
- `h2d_cosTheta_phi` -- 2D `(cos theta*, phi*)` map; primary diagnostic for
  identifying which effect is active.
- `h1d_cosTheta`, `h1d_phi` -- 1D projections.

**Ring proxy along the beam axis (z-hat proxy):**
- `h1d_ringProxy` -- per-Lambda distribution of the ring observable computed
  with the beam direction (+z) as the jet axis.
- `pRingProxyVsEta` -- TProfile of `<R_proxy>` vs Lambda pseudorapidity.
- `pRingProxyVsPt` -- TProfile of `<R_proxy>` vs Lambda pT.
- `pRingProxyIntegrated` -- single-bin TProfile giving the integrated mean.
- `hRingProxy_Kahan` -- 1-bin TH1D storing the high-precision integrated mean and its standard error, computed via a Kahan compensated accumulator to prevent floating-point precision degradation over massive event loops.

**Ring proxy along a random jet direction:**
- `h1d_ringProxyJet` -- per-Lambda distribution of the ring observable using
  a randomly drawn jet direction (reshuffled every `nLambdasSinceJetShuffle`
  Lambdas to mimic a realistic jet-Lambda association).
- `pRingProxyJetVsEta` -- TProfile of `<R_proxyJet>` vs Lambda pseudorapidity.
- `pRingProxyJetVsPt` -- TProfile of `<R_proxyJet>` vs Lambda pT.
- `pRingProxyJet` -- single-bin TProfile giving the integrated mean.
- `pRingProxyJet_JetEtaPos` -- TProfile tracking integrated `<R_proxyJet>` for configurations where jet pseudorapidity `eta_jet >= 0` across all Lambda eta.
- `pRingProxyJet_JetEtaNeg` -- TProfile tracking integrated `<R_proxyJet>` for configurations where jet pseudorapidity `eta_jet < 0` across all Lambda eta.
- `hRingProxyJet_Kahan` -- 1-bin TH1D storing the Kahan-compensated integrated mean.
- `hRingProxyJet_JetEtaPos_Kahan`, `hRingProxyJet_JetEtaNeg_Kahan` -- Kahan-summed equivalents of the jet eta-split profiles.

**Error estimation for the integrated jet ring proxy:**

Three parallel estimators are computed and stored to allow cross-checking of
statistical uncertainty, motivated by the observation that kinematic cuts
(min-pT, DCA) create preferentially sampled regions of phase space that
introduce correlations between Lambdas and can cause the naive TProfile SEM to
underestimate the true uncertainty:

1. **TProfile SEM / Kahan SEM** (`pRingProxyJet`, `hRingProxyJet_Kahan`): the standard SEM, which assumes all Lambda entries are statistically independent. Tends to underestimate the error, specially once cuts are active because DCA and min pT cuts select preferential regions of the Lambda phase space. The Kahan TH1D accumulator uses the exact same sample variance logic but cleanly prevents precision loss during summation.
2. **Per-event mean accumulators** (`hEventMeanRingProxyJet`,
   `pEventMeanRingProxyJetIntegrated`): Lambdas that share the same randomly
   drawn jet direction are grouped into an implicit *event* (jet group). At each
   jet reshuffle, the mean `R_e = evtSumRpj / evtCntRpj` is computed and filled
   (with a weight equal to the number of Lambdas that passed the cuts in that
   event) into `hEventMeanRingProxyJet`. The final uncertainty is
   `StdDev(R_e) / sqrt(N_events)`, which correctly accounts for intra-event
   correlations arising from all Lambdas in an event sharing the same jet
   direction. See `FlushEventMean()` and `FlushEventMeansFamily()`.
3. **Data-chunking estimator** (`hChunkMeansRingProxyJet`,
   `pChunkMeansRingProxyJetIntegrated`): the full run is partitioned into
   `kChunks = 200` sequential chunks, with jet groups assigned round-robin
   (chunk index = `evtIdxRpj % kChunks`). Crucially, the chunk index advances
   for every jet-group event regardless of whether any Lambda survived the cuts,
   keeping the assignment consistent across all twelve `ScenarioHistos`. After
   the main loop, `FinalizeChunks()` computes the per-chunk mean
   `mu_k = chunkSumRpj[k] / chunkCntRpj[k]` for each non-empty chunk and fills
   it into `hChunkMeansRingProxyJet`. The final uncertainty is then
   `StdDev(mu_k) / sqrt(K_filled)`, i.e. `hChunkMeansRingProxyJet->GetMeanError()`.
   The motivation for this third level is to capture the phase-space covariance
   introduced by kinematic cuts across events, which the per-event estimator
   (limited to ~4 Lambdas per event) cannot fully resolve.

**Polarization vector field:**
- `pPstarX_vsPxPy`, `pPstarY_vsPxPy`, `pPstarZ_vsPxPy` -- TProfile2Ds of
  `<p*_x>`, `<p*_y>`, `<p*_z>` in bins of `(px_Lambda, py_Lambda)` over
  `[-3, 3] GeV/c` in a 40×40 grid. Together they form a 2D vector field in
  Lambda transverse momentum space. The expected ring-like structure from DCA
  cuts appears as arrows tangent to circles of constant pT.
- `pPstarX_vsPhiLam`, `pPstarY_vsPhiLam`, `pPstarZ_vsPhiLam` -- 1D TProfiles
  of the same three components vs the Lambda azimuthal angle `phi_Lambda`, in
  32 bins over `[-pi, pi]`. These directly show the sinusoidal modulation around
  the ring without requiring any vector-field rendering.
- `pPstarX_vsPzPx`, `pPstarZ_vsPzPx`, `pPstarY_vsPzPx` -- TProfile2Ds displaying the proton rest-frame direction vector field mapped in the longitudinal-transverse ZX momentum plane against `(pz_lam, px_lam)`. Horizontal components display `<p*_z>`, vertical components display `<p*_x>`, and the out-of-plane background color map charts `<p*_y>`.
- `pPstarX_vsEtaLam`, `pPstarY_vsEtaLam`, `pPstarZ_vsEtaLam` -- 1D compact TProfile projections tracking the individual `<p*>` components against Lambda pseudorapidity `eta_lam` to reveal longitudinal structures driven by kinematics.

**Daughter kinematics:**
- `h1d_pT_proton`, `h1d_pT_pion` -- daughter pT distributions.
- `h1d_DCA_proton`, `h1d_DCA_pion` -- daughter DCA_xy to PV distributions.
- `h1d_decayR` -- transverse decay radius.

**Run natively (compiled - recommended):**
```bash
# Compile natively with optimizations for best performance
g++ -O3 -std=c++20 -march=native -flto -pipe helicityEfficiencyToyModel.cxx -o helicityToyModel $(root-config --cflags --libs)

# Run the compiled binary (an example -- if you really want to use all the options I heavily recommend you stick to the shell coordinator script! It is much more organized!)
./helicityToyModel 2000000 "output.root" 0.5

```

**Run via ROOT CLING (interpreted, or at least JIT-compiled):**

```bash
# Default (all default parameters):
root -l -b -q 'helicityEfficiencyToyModel.cxx'

# Custom parameters (see preamble in the file for full signature):
root -l -b -q 'helicityEfficiencyToyModel.cxx(2000000,"output.root",0.5,...)'

```

---

### `plotHelicityEfficiency.cxx`

**Run this file after the generator.**

Reads the ROOT file produced by `helicityEfficiencyToyModel.cxx` and generates
diagnostic figures organized into output subdirectories within the plots
ROOT file. The output structure per family is:

```
<FamilyDir>/
  EmissionAngles/    -- Figs 1–3
  RingZHat/          -- Figs 4–6
  Daughters/         -- Figs 7–8
  RingJet/           -- Figs 10–15 (including Fig 12b split by jet eta)
  PolarizationVectorField/  -- Figs 16–19

```

Fig 9 (eta asymmetry) is saved directly in the family root directory.

#### Figure descriptions

**Figs 1–3 (EmissionAngles):** one canvas per eta selection (`EtaPos`, `EtaNeg`, `All`).

- **Fig 1**: 2D `(cos theta*, phi*)` maps for each cut scenario. The primary
diagnostic for identifying which effect is active and the distortion it
introduces in the full decay angular distribution.
- **Fig 2**: `cos(theta*)` projections for all four cut scenarios overlaid.
- **Fig 3**: `phi*` projections with a `sin(phi*)` fit to the DCA-cut
distribution. Directly measures the left-right modulation amplitude.

**Figs 4–6 (RingZHat):** ring observable using the beam direction (+z) as
the jet axis.

- **Fig 4**: `<R_proxy>` vs Lambda pseudorapidity -- the **antisymmetric eta
dependence** is the clearest signature of the DCA-cut fake signal.
- **Fig 5**: `<R_proxy>` vs Lambda pT for all cut scenarios.
- **Fig 6**: Integrated `<R_proxy>` bar chart per cut scenario and eta
selection, using TProfile SEMs. Overlaid with mathematically "lossless" (or at least much less loss-prone) Kahan sum values to evaluate if standard TProfile floating-point swamping occurred.

**Figs 7–8 (Daughters):**

- **Fig 7**: Daughter pT distributions for each cut scenario.
- **Fig 8**: Daughter DCA_xy-to-PV distributions for each cut scenario.

**Fig 9** (saved at the family directory level):

- **Fig 9**: `(eta > 0) - (eta < 0)` difference plots of the 2D
`(cos theta*, phi*)` distributions. Isolates the antisymmetric fake component
that would appear as a non-zero ring observable in the data.

**Figs 10–15 (RingJet):** ring observable using a randomly drawn jet direction.

- **Fig 10**: `<R_proxyJet>` vs Lambda pseudorapidity for all four cut scenarios.
- **Fig 11**: `<R_proxyJet>` vs Lambda pT for all four cut scenarios.
- **Fig 12** (two variants):
  - TProfile SEM variant: integrated `<R_proxyJet>` bar chart per cut scenario,
    one panel per eta selection, using the standard TProfile error.
  - Event-mean variant (`IntegratedRingJetEvt`): same layout, but uncertainty
    taken from `hEventMeanRingProxyJet->GetMeanError()`, which accounts for
    intra-event correlations.


- **Fig 12b**: Integrated `<R_proxyJet>` bar chart layout split into three distinct panels based on jet pseudorapidity selections (`eta_jet >= 0`, `eta_jet < 0`, and All `eta_jet`) to directly confirm spatial jet-axis dependence.
- **Fig 13**: Left panel -- `R_proxyJet` per-Lambda distribution for all
  scenarios, normalised and overlaid. Right panel -- `<R_proxyJet>` vs
  Lambda pseudorapidity, all scenarios.
- **Fig 14**: Three-way error comparison per scenario. For each of the four
  cut scenarios (and three eta selections) the figure overlays three
  estimators of the integrated mean: TProfile SEM (filled marker, shifted
  left), event-mean sigma (open marker, centred), and chunking sigma
  (5-pointed star, shifted right). Kahan comparisons are also tracked here. A ratio > 1 (event or chunk vs TProfile)
  means the TProfile was optimistic.
- **Fig 15**: Distribution of per-event (per jet-group) means
  `R_e` for all four cut scenarios. Expected to be much narrower than the
  per-Lambda distribution by ~ `1/sqrt(N_per_event)` and approximately
  Gaussian. Its standard deviation / sqrt(N_events) is the event-mean
  uncertainty plotted in Figs 12 and 14.

**Figs 16–19 (PolarizationVectorField):**

- **Fig 16** (main + supplemental): Spurious `<p*>` polarization vector field
  in the `(px_Lambda, py_Lambda)` plane. The colormap shows `<p*_z>` (scaled to
  percent); overlaid arrows show the mean transverse direction `<p*_T>` in each
  tile, scaled via a robust percentile-based normalization. Arrow length is
  normalized to the **95th-percentile tile magnitude** (`scalePercentile = 0.95`
  by default): tiles whose averaged transverse magnitude exceeds that reference
  are drawn at full length (capped); all others are drawn proportionally. Tiles
  with fewer than `minEntries = 50` occupied bins, or whose averaged magnitude
  is below 5% of the reference, are suppressed. An annotation box in the upper-
  left corner reports the percentile reference magnitude (mean and propagated
  error) as a percentage. The supplemental canvas shows the same field for an
  intermediate set of cuts.
- **Fig 17**: `<p*_x>`, `<p*_y>`, `<p*_z>` as a function of Lambda azimuthal
  angle `phi_Lambda`, one panel per component and one curve per cut scenario.
  Directly shows the sinusoidal modulation that drives the ring-like structure
  in Fig 16.
- **Fig 18**: Main proton rest-frame polarization vector field rendered in the longitudinal-transverse ZX momentum plane `(pz_lam, px_lam)` for the "No Cuts" and "Both Cuts" conditions. The background diverging palette encodes the out-of-plane `<p*_y>` component, while overlaid block-averaged arrows map the mean 2D proton direction (`<p*_z>`, `<p*_x>`) to diagnose longitudinal structures.
- **Fig 18s**: Supplemental ZX-plane polarization vector field displaying the behavior of the two intermediate cut configurations (`pTCutOnly` and `DCACutOnly`).
- **Fig 19**: 1D projections showing `<p*_x>`, `<p*_y>`, and `<p*_z>` components vs. Lambda pseudorapidity `eta_lam` across all four cut scenarios, tracking beam-direction modulations or kinematic biases introduced by selection rules.

**Run with ROOT:**

```bash
root -l -b -q 'plotHelicityEfficiency.cxx("output.root","plots_dir")'

```

---

## `runHelicityToyModel.sh`

**The coordinator script -- start here for a full parameter scan.**

Runs the generator and plotter across eleven families (0–10) of parameter
variations designed to isolate one effect by silencing the other. It organizes
all output into a structured directory tree under a configurable base directory.

Before dispatching the generator jobs, the script automatically pre-compiles `helicityEfficiencyToyModel.cxx` using `g++` with "aggressive" (but still, should by IEEE-standard safe) optimizations (`-O3`, `-march=native`, `-flto`, etc.). The parallel workers then execute this compiled binary rather than interpreting the script via ROOT, yielding significant speed improvements. The plotting step remains I/O bound and is still executed via ROOT CLING.

Supports parallel execution (default `MAX_PARALLEL=64`) and produces one log
file per run.

#### Command-line options

| Flag | Alias | Description |
|------|-------|-------------|
| `--family N [M ...]` | `-f` | Run only the specified family(ies). Multiple integers accepted. |
| `--jobs N` | `-j` | Override the maximum number of concurrent ROOT processes. |
| `--dry-run` | `-n` | Print all commands without executing. |
| `--list` | `-l` | Print the full job table and exit without running. |
| `--plot` | `-p` | Skip the generator; run the plotter only. Useful when only `plotHelicityEfficiency.cxx` has changed. |
| `--help` | `-h` | Print usage and exit. |

**Quick start:**

```bash
chmod +x runHelicityToyModel.sh

# Run all families (may take ~15 min with MAX_PARALLEL=64 on a 192-core machine):
./runHelicityToyModel.sh

# Run only one specific family:
./runHelicityToyModel.sh --family 1

# Run families 2 and 3 simultaneously:
./runHelicityToyModel.sh --family 2 3

# Override parallelism:
./runHelicityToyModel.sh --jobs 8

# Preview all commands without executing:
./runHelicityToyModel.sh --dry-run

# List all registered jobs and their parameters, then exit:
./runHelicityToyModel.sh --list

# Skip data generation and re-run only the plotter:
./runHelicityToyModel.sh --plot

```

When `--plot` is set, the generator step is skipped for every job; the plotter
always runs regardless (it is fast). This is exported as `PLOT_ONLY` so that
all subshells dispatched via GNU parallel or the bash background pool see the
same flag.

#### Scan families

| Family | Subdir | Varied parameter | Key question |
|--------|--------|-----------------|--------------|
| 0 | `0_Baseline/` | Baseline | Flat reference point with no cuts. |
| 1 | `1_AsymDCA/` | Asymmetric DCA cuts | Does the AEE sign depend on which daughter's DCA dominates? |
| 2 | `2_SymDCA/` | Symmetric DCA cuts | How does AEE magnitude grow as a function of cut strength? |
| 3 | `3_BField/` | Magnetic field (DCA cuts only, pT = 0) | Does the asymmetry flip with B-field polarity and scale with strength? |
| 4 | `4_LamPtMin/` | Lambda pT minimum | How does the combined fake signal depend on the Lambda pT regime? |
| 5 | `5_DaughterPt/` | Daughter pT cuts | How does HEE isolate and grow with symmetric/asymmetric pT cuts? |
| 6 | `6_EtaMax/` | Eta acceptance window | How does the fake signal magnitude change with detector acceptance? |
| 7 | `7_Temperature/` | Boltzmann temperature | How does the Lambda pT spectrum shape affect the fake signal? |
| 8 | `8_KinWindow/` | Ring kinematic windows | Where is the fake signal largest across successive pT windows? |
| 9 | `9_RealisticAlice/` | Realistic ALICE cuts | Combined HEE+AEE estimate under experimental conditions. |
| 10 | `10_BField/` | Magnetic field (with pT and DCA cuts) | Combined AEE+HEE field dependence with both cuts simultaneously active; sign flip and strength scan. |

Families 1, 2, 3 are **AEE probes** (DCA cuts active, pT cuts = 0). Family 10
re-runs the field scan from Family 3 but with pT cuts enabled, allowing the
combined AEE+HEE dependence on field polarity and strength to be assessed.
Families 5 and 6 are **HEE probes** (pT cuts active, DCA cuts = 0). Families
4, 7, 8 and 9 use standard-ish ALICE cuts to study kinematic context and
combined effects.

#### Parallelism

All jobs from all selected families are collected into a single queue and
dispatched simultaneously up to `MAX_PARALLEL` concurrent processes. Two
dispatch strategies are tried in order:

1. **GNU parallel** (`parallel --jobs N`): preferred; handles job tracking,
per-job logging via `--joblog`, and clean failure reporting.
2. **Pure bash background pool**: fallback if GNU parallel is not installed;
uses a lightweight semaphore.

A job log file is written to `LOG_DIR/<run_name>.log` for each run (generator
and plotter combined). If GNU parallel is used, a `parallel_joblog.tsv` timing
file is also produced in `LOG_DIR/`.

**Default output directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/`  
**Default log directory:** `/home/users/cicerodm/RingPol/HelicityToyModel/logs/`

Both can be changed at the top of the script. The entire scan generates
approximately 600 MB of ROOT files (~59 runs × ~4–8 MB each; eyeballed).

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
      axis for the z-hat proxy family, and a randomly drawn jet direction
      for the jet proxy family. A more realistic proxy would sample jet
      directions from a measured or simulated jet spectrum with realistic
      eta/phi coverage.
- [ ] The chunking error estimator (`kChunks = 200`) is designed to capture
      phase-space covariance introduced by kinematic cuts, but the optimal
      number of chunks and their relationship to the true covariance
      structure has not been formally derived. Comparisons with the
      TProfile SEM and per-event estimator (Fig 14) are the current
      validation tool.