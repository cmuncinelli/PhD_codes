# InjectedPolO2MC -- injected-polarization full-detector MC (O2 / O2DPG)

> [!WARNING]
> **THIS IS A PERSONAL TEST / DESIGN SANDBOX. NOT AN OFFICIAL ALICE PRODUCTION SETUP.**
>
> Nothing in this directory has been reviewed by the ALICE Simulation/DPG experts or by the O2DPG
> maintainers. It changes the MC truth: Λ decays are done by a custom afterburner, not by Geant4.
> Its output must **never** be used as, or confused with, an official MC production.
> If this proves useful, the generator should be proposed upstream to O2DPG and reviewed there.

## References (what this work is built on)

Repositories inspected to design this (read-only):

| Repository | What we took from it |
|---|---|
| [AliceO2Group/AliceO2](https://github.com/AliceO2Group/AliceO2) | Generator event loop and units ([`Generators/src/Generator.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Generators/src/Generator.cxx), [`Generator.h`](https://github.com/AliceO2Group/AliceO2/blob/dev/Generators/include/Generators/Generator.h)); generator-level decay bookkeeping ([`Generators/src/EventDecayer.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Generators/src/EventDecayer.cxx), `EventDecayer::decayOne`); vertex handling and zero polarization of primaries ([`Generators/src/PrimaryGenerator.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Generators/src/PrimaryGenerator.cxx)); Pythia8 import and seeding ([`Generators/src/GeneratorPythia8.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Generators/src/GeneratorPythia8.cxx)); stack and process codes ([`Detectors/Base/src/Stack.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Detectors/Base/src/Stack.cxx), [`DataFormats/simulation/include/SimulationDataFormat/MCTrack.h`](https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/simulation/include/SimulationDataFormat/MCTrack.h)); physical-primary definition ([`DataFormats/simulation/src/MCUtils.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/DataFormats/simulation/src/MCUtils.cxx)); AO2D MC flags ([`Detectors/AOD/src/AODMcProducerHelpers.cxx`](https://github.com/AliceO2Group/AliceO2/blob/dev/Detectors/AOD/src/AODMcProducerHelpers.cxx), [`Framework/Core/include/Framework/AnalysisDataModel.h`](https://github.com/AliceO2Group/AliceO2/blob/dev/Framework/Core/include/Framework/AnalysisDataModel.h)); macro loading ([`Common/Utils/include/CommonUtils/ConfigurationMacroHelper.h`](https://github.com/AliceO2Group/AliceO2/blob/dev/Common/Utils/include/CommonUtils/ConfigurationMacroHelper.h)) |
| [AliceO2Group/O2DPG](https://github.com/AliceO2Group/O2DPG) (fork: [cmuncinelli/O2DPG](https://github.com/cmuncinelli/O2DPG)) | Workflow creation/running ([`MC/bin/o2dpg_sim_workflow.py`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/bin/o2dpg_sim_workflow.py), [`o2_dpg_workflow_runner.py`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/bin/o2_dpg_workflow_runner.py), [`o2dpg_sim_workflow_anchored.py`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/bin/o2dpg_sim_workflow_anchored.py), [`MC/run/ANCHOR/anchorMC.sh`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/run/ANCHOR/anchorMC.sh)); **the `importParticles()`-override pattern for forced decays** ([`MC/config/PWGEM/external/generator/Generator_pythia8_forcedDecays.C`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGEM/external/generator/Generator_pythia8_forcedDecays.C)); default Pythia card ([`MC/config/common/pythia8/utils/mkpy8cfg.py`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/common/pythia8/utils/mkpy8cfg.py)); Geant4 external decayer for Λ, the rejected alternative ([`MC/run/PWGLF/run_DoubleLambdaTriggered.sh`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/run/PWGLF/run_DoubleLambdaTriggered.sh), [`MC/config/PWGLF/pythia8/decayer/g4_ext_decayer_lambda.in`](https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGLF/pythia8/decayer/g4_ext_decayer_lambda.in)) |
| [AliceO2Group/O2Physics](https://github.com/AliceO2Group/O2Physics) / [cmuncinelli/O2Physics @ LambdaPolRing-v4](https://github.com/cmuncinelli/O2Physics/tree/LambdaPolRing-v4) | Analysis conventions to mirror: α values and ring axis (`PWGLF/Tasks/Strangeness/lambdaJetPolarizationIonsDerived.cxx`); jet/leading-particle definitions (`PWGLF/TableProducer/Strangeness/lambdaJetPolarizationIons.cxx`); V0 MC matching requirements (`PWGLF/TableProducer/Strangeness/strangenessbuilder.cxx`, daughters with `getProcess() == 4`); particle-level jet finders (`PWGJE/JetFinders/jetFinderMCPCharged.cxx`) |
| [gianniliveraro/ALICE_PhotonReconstruction](https://github.com/gianniliveraro/ALICE_PhotonReconstruction) | A working local O2DPG production setup (`GenProduction/`: external generator macro + ini + batch scripts). Its flags were designed for a different analysis (topological Σ⁰ reconstruction); see DESIGN.md §6 for which ones we do and do not reuse |
| [alisw/alidist](https://github.com/alisw/alidist) | `o2.sh` lists `fastjet` as a hard dependency, so FastJet is present in any O2/O2sim environment |
| [cmuncinelli/PhD_codes](https://github.com/cmuncinelli/PhD_codes) (this repo) | Polarized decay kernel from [`vorticity_hydro/lambda_pol_toy_model.cxx`](../vorticity_hydro/lambda_pol_toy_model.cxx) (`sample_P_angle_proton_from_cos_xi`, `Lambda_decay`) |

Upstream links point to the default branches. The exact commits inspected are listed at the top of [`DESIGN.md`](DESIGN.md).

## Goal

Produce full-detector pp MC (Pythia8 → our afterburner → Geant4 → digitization → reconstruction → AO2D)
in which the Λ/Λ̄ carry a **known, injected polarization**. We then analyse the AO2Ds with the same
O2Physics workflow used on real data (as if they came from Hyperloop) to run **full closure tests**.
These cover our ring observable and its correction framework, and also the correction methods used in
other polarization papers:

- **no polarization** → every estimator (R_z, STAR-like signal extraction, correction framework) must recover zero;
- **ring-like polarization** P* = P₀ t̂, with t̂ = (p̂_jet × p_Λ)/|p̂_jet × p_Λ| → compare injected vs reconstructed magnitude;
- **contaminants** (global, longitudinal, helicity-like, ...) → measure how they leak into the ring observable.

Reweighting a standard MC is explicitly *not* enough. Most of the fake polarization comes from
reconstruction, so the injected polarization has to go through the real decay topology, Geant4 and
reconstruction. See the paper roadmap, section 4, "Closure-test for the framework".

## Status

| Step | Status |
|---|---|
| Archaeology of O2 / O2DPG / O2Physics | done, see [`DESIGN.md`](DESIGN.md) |
| Afterburner design | done (v1), see [`DESIGN.md`](DESIGN.md) |
| Decay kernel + polarization models (`generator/InjectedPolKernel.h`) | written; **unit tests pass** (`tests/testKernel.cxx`, run with plain g++) |
| Pythia-side logic (`generator/InjectedPolPythia.h`) | written; compiled and parser-tested against a Pythia API stub only |
| O2 generator macro (`generator/generator_pythia8_injectedPol.C`) | written; **not yet run in O2** |
| Step 0: Pythia-only truth closure (`tests/standalonePythiaClosure.cxx`) | to run on the cluster |
| Step 1: o2-sim smoke test (`run/smokeTestSim.sh` + `checks/checkKine.C`) | to run on the cluster |
| Step 2: small full chain → AO2D (`run/runInjectedPol.sh` + `checks/checkAODFlags.C`) | to run on the cluster |
| Step 3: zero-polarization baseline vs standard Geant4 decays | not started |
| Step 4: injected-polarization productions + closure | not started |

## Layout

```
InjectedPolO2MC/
├── README.md                              <- this file
├── DESIGN.md                              <- archaeology findings, design decisions, simplifications, validation plan
├── generator/
│   ├── InjectedPolKernel.h                <- ROOT-free: polarized 2-body decay, decay vertex, polarization models, constants
│   ├── InjectedPolPythia.h                <- Pythia8-based: config parser, generator-level jet proxies, pre-selection, per-Λ decay
│   ├── generator_pythia8_injectedPol.C    <- O2 external generator (GeneratorPythia8 subclass), MC bookkeeping
│   ├── pythia8_pp_13600.cfg               <- Pythia card (same as O2DPG's default pp card; Λ left undecayed; no seed lines)
│   ├── injectedPol.ini                    <- example [GeneratorExternal]/[GeneratorPythia8] configuration
│   └── configs/                           <- model configs: none, ring (leadP / jet), global, longitudinal, helicity, 5-Λ test
├── run/
│   ├── envcheck.sh                        <- environment checks, path normalisation, log summaries, compiled checkKine (sourced)
│   ├── runTests.sh                        <- unit tests + step 0 (stand-alone Pythia closure), builds *.exe in ${OUT_BASE}/bin
│   ├── smokeTestSim.sh                    <- step 1: o2-sim only, 10 events, then checkKine.C
│   └── runInjectedPol.sh                  <- step 2+: full O2DPG chain to AO2D, then checkKine.C + checkAODFlags.C
├── checks/
│   ├── checkKine.C                        <- bookkeeping + truth-level closure on o2sim kinematics
│   └── checkAODFlags.C                    <- McParticles flags / process codes / eventPlaneAngle in AO2D.root
└── tests/
    ├── testKernel.cxx                     <- unit/closure tests of the kernel (g++ only)
    └── standalonePythiaClosure.cxx        <- step 0: afterburner on plain Pythia8 events, truth-level closure and cost
```

## Quick start (on the cluster)

**Environment.** The package `O2DPG` in aliBuild contains only scripts and has no runtime dependencies, so
`alienv enter O2DPG/...` alone does **not** give you O2 (no `o2-sim`, no `O2_ROOT`). What each step needs:

| Step | Needs | Example |
|---|---|---|
| kernel unit tests | g++ only | any shell |
| Pythia-only closure | Pythia8 (`pythia8-config`) | `alienv enter O2Physics/latest` |
| o2-sim smoke test | O2 | `alienv enter O2Physics/latest` |
| full chain to AO2D | O2 **and** O2DPG | `alienv enter O2Physics/latest O2DPG/latest-master-o2` |

`O2sim` (`aliBuild build O2sim --defaults o2`) is the official MC meta-package: O2Physics + O2DPG + QualityControl +
extra generators. Use it if mixing separately built O2Physics and O2DPG versions causes trouble.
The run scripts check this (`run/envcheck.sh`) and print the `alienv` command to use.

**Paths (no `export` needed).** Every script in `run/` has two hard-coded lines at the top:
```bash
REPO_DIR="/home/users/cicerodm/PhD_codes"            # this repository
OUT_BASE="/home/users/cicerodm/RingPol/InjectedPolMC" # ALL outputs (binaries, logs, simulations, AO2Ds)
```
Edit them there if the cluster paths change. A trailing `/` is stripped automatically. `INJPOL_ROOT` is derived from
`REPO_DIR` inside the scripts (the generator macro needs it); never set it by hand. **Nothing is written into the
repository**: binaries are `*.exe` in `${OUT_BASE}/bin`, ACLiC build products in `${OUT_BASE}/aclic`.

The scripts can be launched from any directory; each one changes to its own output folder:

| Step | Command | Output folder |
|---|---|---|
| unit + 0 | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/runTests.sh [modelCfg] [nEvents]` | `${OUT_BASE}/tests/` (logs), `${OUT_BASE}/bin/` (`*.exe`) |
| 1 | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/smokeTestSim.sh [modelCfg] [nEvents]` | `${OUT_BASE}/smoke/<cfg>/` |
| 2+ | `MODELCFG=<cfg> NSIGEVENTS=50 /home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/runInjectedPol.sh` | `${OUT_BASE}/full/<cfg>_seed<N>/` |

Model configs are in `generator/configs/` (default: `ring_leadP_0p30.cfg` for the tests, `none.cfg` for the full chain).
The very long o2-sim / workflow outputs stay in log files: o2-sim writes `o2sim.log` (steering), `o2sim_serverlog`
(**the generator, including all `[InjectedPol]` messages**), `o2sim_workerlog*` (Geant4) and `o2sim_mergerlog`. The
scripts print only the generator messages, the errors and the check results.

**Known harmless messages** (filtered out of the summaries):
- `Error in <TNetXNGFile::Open>: ... [3011] No servers are available`: a CCDB alignment file not found on one storage
  element. The next line shows it read from another one.
- `G4Exception : mat031 ... FDD_PMMA sum of fractional masses 1.001`: a standard Geant4 material-definition warning.
- `'ENABLE_UPGRADES' macro redefined`, `Particle with PDGcode=... already defined`: ROOT/O2 header noise.

## Where each repository fits

1. **O2DPG** (unmodified, latest release) builds and runs the workflow: `o2dpg_sim_workflow.py -gen external -ini ...`, then
   `o2_dpg_workflow_runner.py -f workflow.json -tt aod`.
2. **This directory** provides the external generator (Pythia8 + polarized-decay afterburner).
3. **O2Physics** (`lambdaJetPolarizationIons` → `lambdaJetPolarizationIonsDerived`) runs on the produced `AO2D.root`
   exactly as on data. A later stage adds dedicated MC tables and processes (DESIGN.md §8).
4. **PhD_codes** post-processing compares injected vs reconstructed polarization.
