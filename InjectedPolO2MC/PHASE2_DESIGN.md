# InjectedPolO2MC -- phase 2 design (DRAFT for discussion)

*2026-11-05. Status: the production chain works end to end (validation steps 0–2 passed on jarvis15 with O2
`daily-20261104`). This document plans what comes next. **Implemented so far: §1 (both items).** Decisions taken
on 2026-11-05 are marked **[DECIDED]**. Decisions marked
**[DECIDE]** need your input. Physics references point to internal note "19 – PolarizationSourcesPaperReview"
(cited as *N19 §x*).*

---

## 0. What is left, in order

| # | Item | Type |
|---|---|---|
| A | Generator fixes found in step 2 (§1) | code, small |
| B | **Step 2b:** run your producer and consumer on the MC AO2D (§3) | test |
| C | Polarization model suite v2 (§2) | code + physics |
| D | Publish productions into the `RingPol` storage layout and registry (§4) | code, small |
| E | **Step 3:** method systematic, generator-level decays (`none.cfg`) vs a standard production (DESIGN.md §5) | production + analysis |
| F | Statistics and cost plan (§5), then **step 4:** injected productions and closure | production + analysis |
| G | Deployment guide: jarvis15 first, then machine by machine | docs + scripts |
| H | Later: O2Physics MC stage (MC tables in your datamodel, `processMC`, generator-level jets on McParticles) | O2Physics |

---

## 1. Generator fixes from the step-2 run

1. **Λ decaying outside the simulated world.** The full chain had decays up to 622 m from the vertex (very forward,
   high-momentum Λ). With `--skipModules ZDC` the world (`cave`) is a box of ±20 m × ±20 m × ±30 m
   (`AliceO2/Detectors/Passive/src/Cave.cxx:76-79`). In a standard production such a Λ simply leaves the world
   undecayed; with us its daughters start outside it and are kept in the MC record.
   **[DECIDED, implemented]** If the sampled decay vertex is outside a configurable box (`decay.onlyInsideWorld`,
   `decay.worldHalf{X,Y,Z}cm`, default = the cave above), the Λ is not decayed. It stays a non-transported
   physical primary without daughters, the same MC truth as an escaping Λ. The per-event count goes to the MC
   header (`injpol.nEscaped`) and to the generator summary. `checkKine.C` checks per event that the Λ without
   daughters are exactly those (B3e). `checkAODFlags.C` reports them separately from A2/A3.
2. **Known-harmless list:** `Error in <TGeoPcon::ComputeBBox>: Duplicated section ...` is an O2 geometry definition
   message (added to the README and the summary filter).

---

## 2. Polarization model suite v2

### 2.1 Direction fields (from N19 §3, "Symmetry constraints on the direction field")

The local triad of N19 §3.1, built around an axis â (trigger t̂, or parton ĉ):

- p̂_Λ (longitudinal, helicity);
- n̂ = (â × p_Λ)/|â × p_Λ| (ring, azimuthal);
- ê_r = n̂ × p̂_Λ (radial, in the plane of â and p_Λ, away from the axis).

Under parity conservation with unpolarized beams, **only n̂ is allowed** (N19 §3.1, proposition "Allowed direction
field"). The p̂_Λ, ê_r and â projections are parity-odd. They cannot feed the ring mean if the acceptance is
mirror-symmetric (N19 §3.5, proposition "Parity-odd sources feed only parity-odd projections"), which
makes them **controls**: injecting them tests the mirror symmetry of the full detector chain.

### 2.2 Proposed models (all combinable; each with its own amplitude)

| key | field | physics / purpose | default amplitude(s) |
|---|---|---|---|
| `ring` | P₀ n̂ (constant) | closure / response of the full chain (existing) | scan **0, 0.02, 0.05, 0.10, 0.30** (linearity and scale factor) |
| `ringPFF` | P_n(j_T) n̂, with P_n = A·(j_T/j₀)·exp(1 − j_T/j₀) | vacuum polarizing fragmentation function: vanishes on the axis, peaks at j_T = j₀, falls off (N19 §3.2, peak j₀ ≈ 0.3 GeV) | A = 0.02 (STAR: percent level, N19 §1.2, §4.2), j₀ = 0.3 GeV; sign as a parameter (STAR sees a sign change with jet p_T) |
| `ringMedium` | A·exp(−(\|Δφ\| − Δφ₀)²/2σ_φ²)·w(p_T)·n̂, with Δφ the azimuth of the Λ from the jet | hydro vortex ring around a quenched jet (Serenone 2021, Ribeiro 2024, as summarised in N19 §6.1) | **from N19:** Δφ₀ = 1 rad (Ribeiro: the jet-induced signal peaks about 1 rad in azimuth from the jet); w(p_T) = 1 for 0.5 < p_T < 1.5 GeV/c, else 0 (the window of the calculations); A = +0.003 realistic (Serenone: a ring observable of a few tenths of a percent, right-handed, no sign change) and A = +0.05 for closure. **[OPEN]**: σ_φ (proposal 0.4 rad), Δη extent, p_T shape inside the window, centrality/η/s scaling. These are not in N19; they need the papers themselves (N19 §9.3, L3) |
| `helicity` | P_H p̂_Λ | P-odd control (and the weak feed-down signature) | 0.10 |
| `radial` | P_r ê_r | P-odd control | 0.10 |
| `alongTrigger` | P_t â | P-odd control (jet-longitudinal) | 0.10 |
| `productionPlane` | A·f(y_Λ)·(ẑ × p̂_Λ)/\|ẑ × p̂_Λ\|, with f = tanh(y_Λ) | beam-related polarization, odd in rapidity in symmetric systems (N19 §3.4, "The beam"). It probes exactly the η-odd leakage you see in data (tanh(η_jet) dependence) | 0.10 |
| `global` | P_G (−ŷ′ of a random ψ_RP) | reaction-plane contaminant (existing). LHC value about 0.01% (N19 lit. table), so this is inflated for a leakage test | 0.05 |
| `longitudinal` | A sin(2(φ − ψ)) ẑ | reaction-plane contaminant (existing) | 0.05 |

Notes:
- The current v1 amplitudes (global 0.20, longitudinal 0.40) are only stress tests; v2 defaults are the table above.
- **Statistics guide** (single Λ): σ(R) ≈ √3/|α|/√N ≈ 2.3/√N. A 1% signal needs about 5×10⁴ reconstructed Λ
  for 1σ, and about 1.3×10⁶ for 5σ. Injections of 0.10–0.30 therefore calibrate the response with modest
  statistics, while the realistic 1–2% models need the large productions of §5.

### 2.3 Which axis the injection uses (`model.trigger`)

| value | axis | use |
|---|---|---|
| `leadP` | leading charged particle (current default) | matches your leadP proxy |
| `jet` | leading anti-kT charged jet (hadron level, as the producer) | matches your jet proxy |
| `parton` **(new)** | the hard-scattered outgoing parton ĉ from the Pythia event record (before pruning) | physically correct axis of the vacuum PFF (N19 §3.3). With a `jet`/`leadP` analysis this measures the axis dilution D_axis = ⟨n̂_c·n̂⟩ (STAR: 0.86–0.94) |

**[DECIDED]** `jet` for the closure models (`ring`), `parton` for the physics-like models (`ringPFF`, `ringMedium`).

**How ĉ is found (not with FastJet).** Pythia keeps the partons of every scattering in the event record. Before
our importParticles pruning, `computeAxes` runs on the full `mPythia.event`, where:
- status **23** marks the outgoing partons of the hardest subprocess;
- status **33** marks the outgoing partons of the subsequent (MPI) subprocesses. In minimum-bias pp the leading jet
  often comes from one of those.

ĉ is the direction of the outgoing parton (23 or 33) **closest in ΔR to the hadron-level trigger jet**, required
within ΔR < R. If no parton matches (soft or diffractive events), there is no ĉ, and the Λ of that event get P* = 0,
consistent with the dilution default. This uses the parton *before* the shower, i.e. the "initiating parton" of the
factorization picture. Option for later: the last copy before hadronization (Pythia `iBotCopy`), which includes the
recoil from final-state radiation.

### 2.4 Mother-aware injection (N19 §9.3 actionable S2; §5.6 "A correction recipe")

Pythia decays Σ⁰ and Σ(1385) before our afterburner sees the event, so their Λ daughters carry the mother in the
event record. Proposal:
- the models give the polarization of the **parent** (Σ⁰, Σ(1385), or the Λ itself if primary);
- the Λ from a Σ⁰ gets C = −1/3 of it, and from a Σ(1385) C = +1/3 (N19 §5.2, transfer table), applied to the mean spin
  vector as defined there;
- Ξ and Ω are decayed by Geant4, so their Λ stay unpolarized (as now). An option to decay Ξ at generator level too,
  with C_Ξ = (1 + 2γ_Ξ)/3 ≈ 0.9–0.95 and the induced Λ helicity, is listed as a later extension;
- the parent species of every injected Λ is recorded, so the coefficients C_i^eff of N19 eq. (recipe) can be measured
  per species from the MC (N19 §5.6).

**[DECIDED]** Default: **all Λ equal** (`model.motherAware = 0`). Mother-aware injection is a complementary workflow,
switched on per production, not part of the main machinery.

### 2.5 Truth bookkeeping

Unchanged principle (DESIGN.md §8): P* is a deterministic function of the event axes and the Λ (and parent)
kinematics, so it can be recomputed from McParticles plus the event header. The parton axis and the model parameters
go into the MC header (`injpol.*`). ψ_RP stays in `McCollisions.eventPlaneAngle`.

---

## 3. Step 2b: your producer and consumer on the MC AO2D

### 3.1 Is the current AO2D enough? Yes.

The O2DPG AO2D (`full/<prod>/AO2D.root`, merged by the `aodmerge` task, or `tfN/AO2D.root`) contains the same core
tables as a data AO2D: BCs, Collisions, TracksIU/Extra/Cov, FT0/FV0/FDD, and the **V0 and cascade index tables**
written by the secondary vertexer (`svfinder`). It also has the MC tables (McCollisions, McParticles, track and
collision labels).

**No extra production step is needed for V0s.** Your producer chain (`runProducerAsHY.sh`) already runs
`o2-analysis-propagationservice`, which performs the track propagation **and** the strangeness builder (V0s from
the svfinder index table), followed by event selection, mult/cent, TPC/TOF PID and `strangenesstofpid`, exactly as
on data.

### 3.2 What can differ on (unanchored) MC, checked in the O2Physics code (fork `LambdaPolRing-v4`, `34cca7e`)

| item | finding | expectation |
|---|---|---|
| event selection | `EventSelectionModule.h:217-218` explicitly supports unanchored Run 3 MC (`run < 500000`); MC mode is set from the AO2D metadata (`isMC = -1` = auto) | should work |
| RCT flags (`cfgRCTLabel = CBT`) | no RCT object for an MC run, so the service inserts dummy flags with only bit 31 ("RCT missing") set (`EventSelectionModule.h:293-297`). Your `RCTFlagsChecker` is initialised without `checkTableValidity` (default `false`), so bit 31 is **not** checked and MC events **pass** | passes; RCT is simply not exercised in MC |
| TPC PID (`pid-tpc-service`) | fetches `Analysis/PID/TPC/Response` at the BC timestamp. Unanchored MC is dated 2019-01-01 (`HBFUtils.startTime=1546300800000`) | **unknown, the first thing to check in the logs.** If no valid object exists, anchoring (DESIGN.md §9) or a fixed `ccdbTimestamp`/`recoPass` is needed |
| multiplicity/centrality | calibrations are looked up per run/timestamp | to check (pp: only used for QA in your producer) |
| your producer | `processData` runs on MC as on data; `processMC` (commented out) is the later O2Physics stage | works as "MC treated as data" |

### 3.3 Procedure

1. Create a local dataset folder in your RingPol storage, e.g. `/home/users/cicerodm/RingPol/InjPolMC_test/`, with
   `input_data_storage.txt` listing `file:/home/users/cicerodm/RingPol/InjectedPolMC/full/none_seed1/AO2D.root`.
2. Run your `LocalProducers` chain on it with the default producer config, then `run_local_consumers.sh`
   (it already writes a `local_registry_auto.conf` with a dummy train ID).
3. Check: no CCDB fatal in the logs (TPC PID, centrality); number of events after `sel8`; number of V0 candidates and
   Λ in the mass window; the consumer runs.

---

## 4. Publishing productions into your RingPol layout (item D)

### 4.1 Many AO2D files, not one merged file **[DECIDED]**

A single large merged file is fragile (one corruption loses everything) and serialises reading. O2DPG already writes
**one AO2D per timeframe** (`tfN/AO2D.root`). The top-level `AO2D.root` is only a convenience merge (task `aodmerge`,
which also runs `AODBcRewriter` on the merged BC table). Plan:
- **file size** is set by `NSIGEVENTS` per TF. At about 125 kB/event, 1000 events give about 125 MB per AO2D.
  Several TFs run per production job;
- **parallelism**: several production jobs with different seeds (`full/<model>_seed<k>/`), as many at a time as the
  machine allows. This is part of the deployment item G;
- **merge skipped** with `-tt 'aod_[0-9]+'` instead of `-tt aod`, saving disk. The runner regex-matches target names,
  so this keeps the per-TF `aod_N` tasks and drops `aodmerge`;
- `run/publishProduction.sh` collects `tfN/AO2D.root` of each job into
  `<STORAGE_BASE>/<DATASET>/rawAO2Ds/AO2D_<job>_<tf>.root` (hard links when on the same filesystem, else copies).
  Only files that pass `checkAODFlags.C` are published, and they are listed in `<DATASET>/input_data_storage.txt`
  (`file:` paths) for `LocalProducers`.

### 4.2 Registry **[DECIDED]**

A **separate registry file**, `RingPol_RAW_LocalHelpers/DerivedDataHY/mc_registry.conf`, with the same four columns
as `train_registry.conf`. `run_all_wagons.sh` already takes the registry as its first positional argument
(`run_all_wagons.sh [OPTIONS] [REGISTRY]`), so no change to your consumer framework is needed:

```bash
./run_all_wagons.sh mc_registry.conf
```

- **IDs:** reserved range **9000001, 9000002, …**, assigned by `publishProduction.sh` (next free number in
  `mc_registry.conf`). They can never collide with Hyperloop train IDs, and stay unique if a line is later copied into
  `train_registry.conf`.
- **Dataset names:** `InjPolMC_pp13p6_<model>_<tag>`, e.g. `InjPolMC_pp13p6_ring0p30_v1`.
- **Manifest:** each dataset folder gets a `production_manifest.txt` with the model config, seeds, O2/O2Physics/O2DPG
  versions, `confKey_used.txt`, event counts and the list of published files.

## 5. Statistics and cost (to refine after step 2b)

Measured on jarvis15 (`none_seed1`, 97 events, 8 CPUs, machine also busy): transport 1.13 s/event (7 CPUs booked),
4.5 s/event end-to-end. That end-to-end figure is dominated by fixed per-TF costs at this size. Disk: about 4.7 MB
per event in total, of which **AO2D is about 125 kB/event**. Hits (2 MB/event) and kinematics (1.1 MB/event) can be
deleted once the AO2D and the truth checks are done.

Needed from step 2b: **reconstructed Λ per event after your cuts** (N_Λ/ev). Then N_events = N_Λ / (N_Λ/ev), and with
the statistics guide of §2.2:

| goal | reconstructed Λ needed |
|---|---|
| response of a 0.30 injection to ±0.01 | about 5×10⁴ |
| 1% realistic signal at 5σ | about 1.3×10⁶ |

**Next measurement:** one production with `NSIGEVENTS=500`, `NTIMEFRAMES=4` on a quiet machine, to get the
asymptotic per-event cost before planning the large runs and the deployment (item G).

---

## 6. Step 2b results so far (2026-11-05)

- The producer chain runs on the MC AO2D, apart from `mult-cent-table`, which **segfaults in `Init()`**. That's not
  MC-specific; the same happened on your Pb–Pb data test. The JSON's `mult-cent-table.enabledTables` was written by
  an older O2Physics: 39 tables, including `FV0AOuterMults`. The current `MultModule.h` expects 42 (`FITExtraMults`,
  `CentFT0MAnchorCols`, `CentFT0MAnchorBCs`, `CentFT0MOuterAs`) and looks each one up by label; a missing label gives
  an invalid index, hence the crash in `Array2D::operator()`. The block was all `-1` (defaults), so the fix is to drop
  it: `O2PhysicsTools/dropDefaultEnabledTables.py` (run it after every O2Physics update). `propagation-service`'s
  list (37 tables) still matches the code.
- Still to read from the re-run: CCDB messages of `pid-tpc-service` and the centrality lookup on unanchored MC; the
  counts of events after `sel8`, V0 candidates and Λ in the mass window.

## 7. Item H preview: MC tables already within reach

Without any new datamodel, your current chain can produce the standard strangeness MC tables:
`propagation-service` with `processMonteCarlo: true` (instead of `processRealData`) fills `McV0Labels`, `V0MCCores`,
`V0CoreMCLabels`, `V0MCCollRefs` (all listed in your JSON's `enabledTables`). The AO2D also carries `McCollisions`,
`McParticles`, `McCollisionLabels` and `McTrackLabels`. These are the inputs for a later `processMC` in your producer:
true Λ ↔ V0 matching, particle-level jets on McParticles, and recomputing the injected P* from the generator header.
