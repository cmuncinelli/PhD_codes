# InjectedPolO2MC -- injected-polarization full-detector MC (O2 / O2DPG)

> [!WARNING]
> **THIS IS A PERSONAL TEST / DESIGN SANDBOX. NOT AN OFFICIAL ALICE PRODUCTION SETUP.**
>
> Nothing in this directory has been reviewed by the ALICE Simulation/DPG experts or by the O2DPG
> maintainers. It changes the MC truth (Λ decays are done by a custom afterburner and not by Geant4).
> Its output must **never** be used as, or confused with, an official MC production.
> If this proves useful, the generator should be proposed upstream to O2DPG and reviewed there.

## Goal

Produce full-detector (Pythia8 → our afterburner → Geant4 → digitization → reconstruction → AO2D) pp
MC where the Λ/Λ̄ hyperons carry a **known, injected polarization**. We then analyse the AO2Ds with the
same O2Physics workflow used on real data (as if they came from Hyperloop) to run **closure tests** of
the ring-observable measurement and its correction framework:

- **no polarization** → every estimator (R_z, STAR-like signal extraction, correction framework) must recover zero;
- **ring-like polarization** P = P₀ t̂, t̂ = (p̂_jet × p_Λ)/|p̂_jet × p_Λ| → compare injected vs reconstructed magnitude;
- **contaminants** (global, longitudinal, helicity-like polarization...) → study how they leak into the ring observable.

See "Closure-test for the framework" in the paper roadmap, section 4.

## Status

| Step | Status |
|---|---|
| Archaeology of O2 / O2DPG / O2Physics | done, see [`DESIGN.md`](DESIGN.md) |
| Design of the afterburner | first draft, see [`DESIGN.md`](DESIGN.md) |
| Generator macro (`generator_pythia8_injectedPol.C`) | not started |
| Smoke test (o2-sim only, Kine inspection) | not started |
| Small full chain → AO2D, flag checks | not started |
| Zero-polarization baseline vs standard Geant4 decays | not started |
| Injected-polarization productions + closure | not started |

## Planned layout

```
InjectedPolO2MC/
├── README.md                          <- this file
├── DESIGN.md                          <- archaeology findings + afterburner design
├── generator/
│   ├── generator_pythia8_injectedPol.C   <- GeneratorPythia8 subclass (the afterburner)
│   ├── injectedPol.ini                   <- [GeneratorExternal] / [GeneratorPythia8] config
│   └── pythia8_pp_13600.cfg              <- Pythia card (Λ left undecayed)
├── run/
│   ├── runInjectedPol.sh                 <- o2dpg_sim_workflow.py + o2_dpg_workflow_runner.py
│   └── batch/                            <- batch helpers (inspired by gianniliveraro/ALICE_PhotonReconstruction)
└── checks/
    ├── checkKine.C                       <- inspect o2sim_Kine.root (bookkeeping smoke test)
    └── checkAODFlags.C                   <- inspect McParticles flags in AO2D
```

## Pipeline (where each repo fits)

1. **O2DPG** (unmodified, from the O2sim/O2PDPSuite release in use) builds and runs the workflow:
   `o2dpg_sim_workflow.py -gen external -ini InjectedPolO2MC/generator/injectedPol.ini ...` then
   `o2_dpg_workflow_runner.py -f workflow.json -tt aod`.
2. **This directory** provides the external generator (Pythia8 + polarized-decay afterburner).
3. **O2Physics** (`lambdaJetPolarizationIons` → `lambdaJetPolarizationIonsDerived`) runs on the produced
   `AO2D.root` exactly as on local RAW-derived data.
4. **PhD_codes** post-processing compares injected vs reconstructed polarization.
