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

## Notation

Same as in the analysis notes:
- **t̂ = trigger direction**: the leading jet or the leading particle (`model.trigger = jet | leadP`);
- **n̂ = (t̂ × p_Λ)/|t̂ × p_Λ|**: the ring direction, with the same cross-product order as the consumer task;
- **R = (3/α)⟨p̂*_B · n̂⟩**: the ring observable, where p̂*_B is the baryon direction in the Λ rest frame.

## Goal

Produce full-detector pp MC (Pythia8 → our afterburner → Geant4 → digitization → reconstruction → AO2D)
in which the Λ/Λ̄ carry a **known, injected polarization**. We then analyse the AO2Ds with the same
O2Physics workflow used on real data (as if they came from Hyperloop) to run **full closure tests**.
These cover our ring observable and its correction framework, and also the correction methods used in
other polarization papers:

- **no polarization** → every estimator (R_z, STAR-like signal extraction, correction framework) must recover zero;
- **ring-like polarization** P* = P₀ n̂, with n̂ = (t̂ × p_Λ)/|t̂ × p_Λ| → compare injected vs reconstructed magnitude;
- **contaminants** (global, longitudinal, helicity-like, ...) → measure how they leak into the ring observable.

Reweighting a standard MC is explicitly *not* enough. Most of the fake polarization comes from
reconstruction, so the injected polarization has to go through the real decay topology, Geant4 and
reconstruction. See the paper roadmap, section 4, "Closure-test for the framework".

## Status

| Step | Status |
|---|---|
| Archaeology of O2 / O2DPG / O2Physics | done, see [`DESIGN.md`](DESIGN.md) |
| Afterburner design | done (v1), see [`DESIGN.md`](DESIGN.md) |
| Decay kernel + polarization models (`generator/InjectedPolKernel.h`) | **unit tests pass** (here and on jarvis15) |
| Pythia-side logic (`generator/InjectedPolPythia.h`) | **runs with real Pythia** (step 0) |
| O2 generator macro (`generator/generator_pythia8_injectedPol.C`) | **runs inside o2-sim** (step 1) |
| Step 0: Pythia-only truth closure (`run/runTests.sh`) | **passed** (jarvis15, 2026-11-05): ring R recovered within 1σ, with and without a trigger requirement |
| Step 1: o2-sim smoke test (`run/smokeTestSim.sh` + `checks/checkKine.C`) | **passed** (jarvis15, 2026-11-05): B1–B5 all OK on 44 Λ |
| Step 2: small full chain → AO2D (`run/runInjectedPol.sh`) | first run found the O2DPG "separated" event-generation problem (DESIGN.md §4b); fixed with `--event-gen-mode integrated`; **re-run pending** |
| Step 3: zero-polarization baseline vs standard Geant4 decays (method systematic, DESIGN.md §5) | not started |
| Step 4: injected-polarization productions + closure | not started |
| Deployment on the cluster (parallel runs, then other machines) | planned after step 3: a step-by-step, machine-by-machine guide will be added (jarvis15 first) |

## Layout

```
InjectedPolO2MC/
├── README.md                              <- this file
├── DESIGN.md                              <- archaeology findings, design decisions, simplifications, validation plan
├── PHASE2_DESIGN.md                       <- DRAFT plan for what comes next (model suite v2, step 2b, publishing, statistics)
├── generator/
│   ├── InjectedPolKernel.h                <- ROOT-free: polarized 2-body decay, decay vertex, polarization models, constants
│   ├── InjectedPolPythia.h                <- Pythia8-based: config parser, generator-level jet proxies, pre-selection, per-Λ decay
│   ├── generator_pythia8_injectedPol.C    <- O2 external generator (GeneratorPythia8 subclass), MC bookkeeping
│   ├── pythia8_pp_13600.cfg               <- Pythia card (same as O2DPG's default pp card; Λ left undecayed; no seed lines)
│   ├── injectedPol.ini                    <- example [GeneratorExternal]/[GeneratorPythia8] configuration
│   └── configs/                           <- model configs: none, ring (leadP default / leadP requiring a trigger / jet), global, longitudinal, helicity, 5-Λ test
├── run/
│   ├── envcheck.sh                        <- environment checks, path normalisation, log summaries/moves, disk report, compiled checkKine (sourced)
│   ├── progress.py                        <- progress bar, ETA, per-task timing, disk usage, stall/runaway-loopers alarms
│   ├── killRun.sh                         <- stop ONE run (only processes working inside its output folder)
│   ├── runTests.sh                        <- unit tests + step 0 (stand-alone Pythia closure), builds *.exe in ${OUT_BASE}/bin
│   ├── recheckKine.sh                     <- re-run only checkKine.C on existing o2-sim output
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

**Shortcut:** if your `~/.bashrc` contains the block from
[`O2PhysicsTools/O2EnvHelper-bashrcBlockToAppend`](../O2PhysicsTools/O2EnvHelper-bashrcBlockToAppend), just type
**`o2envsim`** (same as `o2env sim`). It enters `O2Physics/latest,O2DPG/latest-master-o2` from the current folder,
which is the environment for every step below. `o2env` alone gives `O2Physics/latest`.

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
| unit + 0 | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/runTests.sh [modelCfg] [nEvents]` | `${OUT_BASE}/tests/logs/`, `${OUT_BASE}/bin/` (`*.exe`) |
| 1 | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/smokeTestSim.sh [modelCfg] [nEvents]` | `${OUT_BASE}/tests/smoke/<cfg>/` |
| 1 (re-check only) | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/recheckKine.sh [dir] [prefix]` | same folder as the checked output |
| 2+ | `MODELCFG=<cfg> NSIGEVENTS=50 /home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/runInjectedPol.sh` | `${OUT_BASE}/full/<cfg>_seed<N>/` |
| monitor | `/home/users/cicerodm/PhD_codes/InjectedPolO2MC/run/progress.py <full-chain folder> [--watch 30] [--final]` | prints; `--final` also writes `logs/run_report.txt` |

**Logs** always go to a `logs/` subfolder of the output folder (`tests/logs/`, `tests/smoke/<cfg>/logs/`,
`full/<...>/logs/`). For the full chain, our own logs are written there directly. The O2DPG task logs (`tfN/*.log`,
`*.log_done`, `sgn_serverlog`, ...) are moved there only **after a successful run**, because the O2DPG runner uses the
`*.log_done` files to skip finished tasks when a failed run is resumed in place.

**Progress, ETA, timing, disk.** While a full chain runs, open a second terminal (same environment not needed, plain
`python3`) and run `progress.py <folder> --watch 30`. It shows:
- a progress bar over the workflow tasks;
- event-level progress and seconds per event of the Geant4 transport (`sgnsim_N`), the dominant cost;
- an ETA once a first timeframe has finished (an upper bound: timeframes run in parallel when CPUs allow);
- at the end, wall time per task kind and disk usage by category (logs, AO2D, kinematics, hits, digits, other),
  total and per event.

Per-event numbers from small runs (< ~500 events) are dominated by fixed per-timeframe costs and **overestimate**
large productions; the report says so. `smokeTestSim.sh` and `runTests.sh` also print their wall times and disk usage.

Model configs are in `generator/configs/` (default: `ring_leadP_0p30.cfg` for the tests, `none.cfg` for the full chain).
**By default no generator trigger is required**: events without one are simulated too, with P* = 0, which gives the
data-like dilution (DESIGN.md §7b). `ring_leadP_0p30_requireTrigger.cfg` keeps only events with a generator trigger.
The very long o2-sim / workflow outputs stay in log files: o2-sim writes `o2sim.log` (steering), `o2sim_serverlog`
(**the generator, including all `[InjectedPol]` messages**), `o2sim_workerlog*` (Geant4) and `o2sim_mergerlog`. The
scripts print only the generator messages, the errors and the check results.

**Known harmless messages** (filtered out of the summaries):
- `Error in <TNetXNGFile::Open>: ... [3011] No servers are available`: a CCDB alignment file not found on one storage
  element. The next line shows it read from another one.
- `G4Exception : mat031 ... FDD_PMMA sum of fractional masses 1.001`: a standard Geant4 material-definition warning.
- `'ENABLE_UPGRADES' macro redefined`, `Particle with PDGcode=... already defined`: ROOT/O2 header noise.
- `Error in <TGeoPcon::ComputeBBox>: Duplicated section ...`: an O2 geometry-definition message, also in standard productions.
- `Failed to import loopers particles` / `ReadEvent failed in finalizeEvent` / `ReadEvent failed for generator ALICEo2`
  (full chain only): **harmless only with the loopers workaround or a recent O2.** In O2 versions from 2026-07-24 up to
  `daily-20260904-0700` each of these desynchronises the TPC-loopers generator, and the run later **hangs**, asking for
  ~1e9 loopers (DESIGN.md §4c). **Fixed in O2 (AliceO2) by [`c575d91`](https://github.com/AliceO2Group/AliceO2/commit/c575d910d57ae0f3fa0a73bf62be3dce683b0de1), first in `daily-20260905-0000`**; the related
  follow-ups [`c756dc9`](https://github.com/AliceO2Group/AliceO2/commit/c756dc946b622cea3d22a3f98b3c432cd8b1c4e6) and [`1f4567d`](https://github.com/AliceO2Group/AliceO2/commit/1f4567dbcfd04f266f6c648997bb5ea7cbffb6e2) are in `daily-20261031-0000` and later. With such an O2 nothing is changed; on older
  versions `runInjectedPol.sh` detects the O2 version and applies the workaround automatically.

**If a run hangs:** `run/progress.py <folder>` flags stalls (no log activity for 20 min) and runaway loopers.
`run/killRun.sh <folder>` lists, and with `--kill` stops, only the processes of that run (working directory inside the
folder); your other jobs and your terminals are never touched. Then delete the folder and start again.
- More events than requested (e.g. 97 for `NSIGEVENTS=50`): O2DPG adds "orbits-early" pile-up collisions before the
  timeframe start. All particles carrying the AO2D flag `FromBackgroundEvent` is also normal without embedding.
- Should you ever see `GeneratorFromO2Kine: Ran out of events`: the generator ran in O2DPG's *separated* mode, which
  breaks this generator. The scripts use `--event-gen-mode integrated` (DESIGN.md §4b).

## Naming used in this directory

| Name | Meaning |
|---|---|
| **afterburner** | Code that runs on a generator's output and modifies it before transport; here, the polarized Λ decays applied after Pythia (HEP jargon: a step "after" the main generator) |
| **closure test** | Run a measurement + correction on MC with a *known* input and check that the input is recovered ("the loop closes") |
| **smoke test** | The shortest run that exercises the whole chain, here 10 events of o2-sim plus the bookkeeping checks. From hardware testing: switch it on and see whether smoke comes out. It catches crashes and gross mistakes, not subtle physics |
| **trigger, t̂** | The direction the ring is defined around: leading jet or leading particle (`model.trigger`). Same meaning as in the analysis notes |
| **n̂** | Ring direction n̂ = (t̂ × p_Λ)/\|t̂ × p_Λ\| |
| **pre-selection** | Generator-level event filter applied before Geant4 (`presel.*`), to save transport time |
| **dilution** | Events/Λ that enter the measured average with zero true ring signal (e.g. a reconstructed trigger without a genuine jet) |
| **Kine / MCHeader** | o2-sim kinematics output (`<prefix>_Kine.root`, all MC particles) and per-event header (`<prefix>_MCHeader.root`) |
| **`o2sim` / `sgn` prefix** | o2-sim output prefix: `o2sim` in the smoke test; `sgn` ("signal") inside each `tfN/` of an O2DPG production |
| **TF** | Timeframe, the O2 unit of data (a fixed number of LHC orbits); O2DPG productions are organised per TF (`tf1/`, `tf2/`, ...) |
| **B1–B5, P1–P2** | Checks in `checks/checkKine.C`: B = bookkeeping, P = physics (truth-level) |
| **A1–A4** | Checks in `checks/checkAODFlags.C`, on the AO2D |
| **steps 0–4** | Validation plan in DESIGN.md §10 (unit tests and step 0: `runTests.sh`; step 1: `smokeTestSim.sh`; step 2+: `runInjectedPol.sh`) |

## Where each repository fits

1. **O2DPG** (unmodified, latest release) builds and runs the workflow: `o2dpg_sim_workflow.py -gen external -ini ...`, then
   `o2_dpg_workflow_runner.py -f workflow.json -tt aod`.
2. **This directory** provides the external generator (Pythia8 + polarized-decay afterburner).
3. **O2Physics** (`lambdaJetPolarizationIons` → `lambdaJetPolarizationIonsDerived`) runs on the produced `AO2D.root`
   exactly as on data. A later stage adds dedicated MC tables and processes (DESIGN.md §8).
4. **PhD_codes** post-processing compares injected vs reconstructed polarization.
