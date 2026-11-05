# InjectedPolO2MC -- phase 2 design (DRAFT for discussion)

*2026-11-05. Status: the production chain works end to end (validation steps 0–2 passed on jarvis15 with O2
`daily-20261104`). This document plans what comes next. **Nothing here is implemented yet.** Decisions marked
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
   **Proposal:** if the sampled decay vertex is outside a configurable box (default = the cave above), do not decay
   the Λ: keep it as a non-transported physical primary without daughters, which is the same MC truth as an
   escaping Λ. This is counted in the generator summary, and `checkKine` B3 accepts 0 daughters only for those Λ.
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
| `ringMedium` | A·exp(−(ΔR − ΔR₀)²/2σ²) n̂ | hydro vortex ring around a quenched jet; angular extent still open (N19 §9.2, Q6) | A = 0.05, ΔR₀ = 0.4, σ = 0.2 **[DECIDE / literature, N19 §9.3, L3]** |
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

**[DECIDE]** the default trigger for each model. My proposal: `jet` for the closure models (`ring`), and `parton`
for the physics-like models (`ringPFF`, `ringMedium`).

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

**[DECIDE]** whether the default is "mother-aware" (physically motivated) or "all Λ equal" (cleanest closure). My
proposal: all Λ equal for the closure productions, and mother-aware for the feed-down study.

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

### 3.3 Procedure (proposal)

1. Create a local dataset folder in your RingPol storage, e.g. `/home/users/cicerodm/RingPol/InjPolMC_test/`, with
   `input_data_storage.txt` listing `file:/home/users/cicerodm/RingPol/InjectedPolMC/full/none_seed1/AO2D.root`.
2. Run your `LocalProducers` chain on it with the default producer config, then `run_local_consumers.sh`
   (it already writes a `local_registry_auto.conf` with a dummy train ID).
3. Check: no CCDB fatal in the logs (TPC PID, centrality); number of events after `sel8`; number of V0 candidates and
   Λ in the mass window; the consumer runs.

---

## 4. Publishing productions into your RingPol layout (item D)

The infrastructure exists already. `LocalProducers` takes `<STORAGE_BASE>/<DATASET_NAME>/input_data_storage.txt`,
and `run_local_consumers.sh` builds a local registry (dummy `TRAIN_RUN_ID 999999`) for any wagon folder with
`AO2Ds/`. Proposal:

- a script `run/publishProduction.sh <production folder(s)> <DATASET_NAME>` that writes
  `/home/users/cicerodm/RingPol/<DATASET_NAME>/input_data_storage.txt` (all `AO2D.root` of the productions), plus a
  `production_manifest.txt` (model config, seeds, O2/O2DPG versions, `confKey_used.txt`, event counts) for
  bookkeeping;
- dataset naming: `InjPolMC_pp13p6_<model>_<tag>`, e.g. `InjPolMC_pp13p6_ring0p30_v1`;
- in `train_registry.conf`: a section `# ---- InjPolMC (local MC, not Hyperloop) ----`. **[DECIDE]** either reuse the
  dummy `999999`, or a reserved range (e.g. `9000001`, `9000002`, ... one per production) so that several MC
  productions can coexist in the same registry without colliding. My proposal: the reserved range, written by
  `publishProduction.sh`.

---

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
