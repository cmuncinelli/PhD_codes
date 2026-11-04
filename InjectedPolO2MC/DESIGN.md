# InjectedPolO2MC -- archaeology findings and afterburner design

*v1, 2026-11-04. Source versions inspected: AliceO2 `a6674c4` (2026-11-03), O2DPG `e9d5f2d` (2026-11-03; your fork is
in sync with upstream), O2Physics fork `master` `8018361` and **`LambdaPolRing-v4` `ab79d01`**, alidist `master`.
File/line references are relative to those trees. Links to all repositories are in [`README.md`](README.md).*

---

## 1. TL;DR

- **Injection point:** an O2DPG-style **external generator macro**, i.e. a subclass of `o2::eventgen::GeneratorPythia8`
  that overrides `importParticles()`. After Pythia fills the event, every final-state Λ/Λ̄ is decayed by us with a
  polarized angular distribution, and the daughters are appended with a **displaced production vertex**. Geant4
  then transports the daughters (and everything else) normally. No change to O2, O2DPG or O2Physics is needed.
- **Recipe already exists in O2/O2DPG:** `AliceO2/Generators/src/EventDecayer.cxx:303-362` and
  `O2DPG/MC/config/PWGEM/external/generator/Generator_pythia8_forcedDecays.C` do generator-level isotropic decays.
  We follow their structure and replace the kinematics with the polarized kernel.
- **Bookkeeping fixed so that MC truth looks standard** (§4): the Λ stays a *physical primary*, and the daughters are
  *secondaries from a decay* (process 4). This is what the strangeness builder needs for MC matching.
- **Known deviation from standard productions:** the Λ does not interact with material before decaying.
  §5 gives a characterization plan.

---

## 2. How an O2DPG MC production works (the part relevant to us)

```
o2dpg_sim_workflow.py  ──► workflow.json ──► o2_dpg_workflow_runner.py -tt aod
                                               │
   per timeframe:  o2-sim (generator + Geant4) │ digitization │ reco (ITS/TPC/TOF/...) │ secondary vertexing
                   │                                                        │ AOD producer ──► tfN/AO2D.root
                   └─ generator = Pythia8 or external macro (.C, loaded with gROOT->LoadMacro, i.e. JIT by Cling)
                      By default O2DPG runs the generator in a SEPARATE o2-sim (task sgngen_N) and the transport re-reads
                      the events; for this generator we need --event-gen-mode integrated (§4b).
```

- An external generator is selected with `-gen external -ini my.ini` ([GeneratorExternal] `fileName`/`funcName`, plus
  [GeneratorPythia8] `config` for a Pythia card). `${ENV}` inside the ini is expanded (O2DPG uses it in both
  `fileName` and `funcName` arguments). Macros include shared headers with `R__ADD_INCLUDE_PATH($ENVVAR/...)`, which
  is what we do with `$INJPOL_ROOT`.
- **Event loop** (`Generators/src/Generator.cxx:308-384`): `generateEvent()` → `importParticles()` → `finalizeEvent()`
  → trigger (loop until accepted) → `decayParticles()` (optional `EventDecayer`, off by default) → `addTracks()`.
- **Units in `mParticles`:** those of the generator. Pythia uses mm and mm/c, converted with `mPositionUnit = 0.1`,
  `mTimeUnit = 3.3356410e-12` (`Generator.h:157-158`).
- **Vertices:** `PrimaryGenerator::AddTrack` adds the smeared interaction vertex to the particle's own production
  vertex (`PrimaryGenerator.cxx:148-150`), so displaced decay vertices set in the generator survive.
- **Transport flag:** `ParticleStatus::kToBeDone`. Non-transported generator particles are still stored (`Stack.cxx:228-247`).
- **Process code:** the generator's `TParticle::UniqueID` becomes the MC creation process (`Stack.cxx:240`,
  `MCTrack.h:399`). `MCTrack::isPrimary()` is `process == kPPrimary || no mothers` (`MCTrack.h:75`).
- **Seeding:** with `-gen pythia8`, O2DPG writes a per-TF card with the TF seed (`o2dpg_sim_workflow.py:913`).
  With `-gen external` it does not. `GeneratorPythia8` then seeds Pythia from `gRandom`, which follows
  `o2-sim --seed` (set per TF by O2DPG), provided the card has **no `Random:` lines** (`GeneratorPythia8.cxx:114-140`).
  The afterburner draws from Pythia's own `rndm`, so a production is reproducible from `-seed` alone.
- **What standard productions do with Λ:** Pythia leaves every particle with τ₀ > 10 mm/c undecayed
  (`ParticleDecays:limitTau0 = on`, `tau0Max = 10`; Λ has cτ = 78.9 mm). The same holds for K⁰_S, Σ±, Ξ, Ω, π±, K±.
  Geant4 transports them through the material, where hadronic interactions are possible, and decays them where they
  actually decay, **isotropically** (phase space, no polarization). The one exception in O2DPG is the PWGLF
  double-Λ production, which forces the charged decay channel through Geant4's external-decayer interface (§3.2).

---

## 3. Injection points considered

### 3.1 Generator-level decay (chosen for v1)

`importParticles()` override (as in about 30 O2DPG macros). The whole event is available (jet axis, multiplicity,
a reaction plane). It needs no core change, runs in a single process, and is reproducible.
Cost: the Λ does not interact with material before its decay (§5), and the decay length is sampled by us
(same cτ as Geant4).

### 3.2 Decay during Geant4 transport through a custom external decayer (possible v2)

There is a precedent: `O2DPG/MC/run/PWGLF/run_DoubleLambdaTriggered.sh` with
`MC/config/PWGLF/pythia8/decayer/g4_ext_decayer_lambda.in` (`/mcPhysics/setExtDecayerSelection lambda anti_lambda`).
Geant4 hands Λ decays to `o2::eventgen::DecayerPythia8`. A polarized `TVirtualMCDecayer` would keep the material
effects exactly standard. It is blocked because:
- `TVirtualMCDecayer::Decay(int pdg, TLorentzVector* p)` gets only the PDG code and momentum, with no event context;
- the polarization of primaries is **hard-coded to zero** (`PrimaryGenerator.cxx:165-167`), so it cannot travel with the Λ;
- `o2-sim` transports chunks of an event in **separate worker processes**, so a shared lookup table fails too.

This would need a small O2 core patch: pass the `TParticle` polarisation through `AddTrack`, and read it in the
decayer from the stack's current track. Worth doing only if §5 shows that the material effect matters.

### 3.3 Reweighting a standard production (rejected)

Weighting candidates by 1 + α P·p̂* (as in D. Scheirich's ATLAS thesis) does not change the topologies seen by the
reconstruction. Fake polarization comes largely from reconstruction, so this cannot replace a real injection.

---

## 4. MC bookkeeping (the critical part)

What the code downstream does:
- `MCTrackNavigator::isPhysicalPrimary` (`DataFormats/simulation/src/MCUtils.cxx:73-...`):
  - false for generator particles with HepMC status > 1;
  - for generator particles, true if "stable" (Λ, p, π are; `MCUtils.h:87-...`), whatever the mother;
  - for transport particles, false if the mother is a light-flavour hadron.
- AOD producer (`Detectors/AOD/src/AODMcProducerHelpers.cxx:229-235`): non-primaries get `ProducedByTransport`, and
  their `statusCode` is replaced by the process code. Decay daughters of a primary are kept (`isKeepPhysics` →
  `isFromPrimaryDecayChain`, process `kPDecay`).
- O2Physics strangeness builder (`strangenessbuilder.cxx:2020`): generated-V0 info only counts daughters with
  `getProcess() == 4`. It also uses the mother's `isPhysicalPrimary()` (`:1874`, `:1999`).
- `McParticles_001` (current) stores the flags as written. The vertex-based "protection" only exists in the old `_000` version.

| | Λ (mother) | p, π daughters |
|---|---|---|
| `EventDecayer` / `forcedDecays` recipe as-is | **not** PhysicalPrimary (HepMC 2) | **PhysicalPrimary**, process 0 |
| Standard production (Geant4 decay) | PhysicalPrimary | secondary, process 4, `ProducedByTransport` |
| **This generator** | PhysicalPrimary (HepMC 1, not transported) | secondary, process 4 (`UniqueID = kPDecay`), `ProducedByTransport` |

Implemented in `generator_pythia8_injectedPol.C` (`injectPolarizedDecays`, `pushDaughter`). It is verified by
`checks/checkKine.C` (B1-B5, on o2sim kinematics) and `checks/checkAODFlags.C` (A1-A4, on the AO2D).
Remaining question for the smoke test: does any downstream code assume "HepMC status 1 ⇒ transported"? None found in
O2 or the PWGLF code read so far.

---

## 4b. O2DPG event-generation mode: `--event-gen-mode integrated` is required

Found in the first full-chain test (jarvis15, 2026-11-05: B3–B5 and A3 failed for every Λ, smoke test fine).

- Current O2DPG defaults to `--event-gen-mode separated` (`o2dpg_sim_workflow.py:148`). Task `sgngen_N` runs the
  generator alone (`o2-sim --noGeant -g external ... -o genevents`), and the transport task `sgnsim_N` re-reads the
  events with `-g extkinO2` (`GeneratorFromO2Kine`).
- `GeneratorFromO2Kine::importParticles` keeps **only primary tracks** unless in "continue mode"
  (`AliceO2/Generators/src/GeneratorFromFile.cxx:282`). Our daughters are deliberately not primary (process `kPDecay`,
  with a mother; §4), so they are dropped. The non-transported Λ then never decays. The skipped particles also leave
  gaps that are not re-indexed, so the Λ's daughter indices point to unrelated particles in the AO2D (A2 passed, A3 failed).
- **Fix:** `--event-gen-mode integrated` (`run/runInjectedPol.sh`). The generator runs inside the transport `o2-sim`,
  exactly as in the smoke test, and the bookkeeping of §4 survives. `checkKine.C` and `checkAODFlags.C` print a hint
  if this ever happens again, and `checkAODFlags.C` A2 now also requires the daughters to point back to the Λ.
- Making the generator compatible with *separated* mode would require the daughters to be "primary" for
  `GeneratorFromO2Kine`. That is exactly the bookkeeping we avoided (it flips the physical-primary flags). Integrated
  mode is an official, supported O2DPG option, so we use it.

Other observations from the same run (standard O2/O2DPG behaviour, not related to our generator):
- *97 events for `-ns 50`:* the collision context holds up to `--maxCollsPerTF` (= `-ns`) collisions in the TF, plus
  "orbits-early" collisions before the TF start (pile-up from earlier bunch crossings), which are simulated too.
- *`Failed to import loopers particles` / `ReadEvent failed ...`:* **not harmless in the O2 version on jarvis15**; see §4c.
  In separated mode each retry also consumed one more pre-generated event, hence `GeneratorFromO2Kine: Ran out of events`.
- *All AO2D particles carry `FromBackgroundEvent` (0x2):* the AOD producer flags source 0 as background
  (`AODProducerWorkflowSpec.cxx:1215`), and `sgn` is source 0 without embedding.

## 4c. TPC loopers: an O2 bug that hangs the run, and the version-aware workaround

Found in the second full-chain test (jarvis15, 2026-11-05): the run hung for 4+ hours at event 94/97. The particle
server had logged `Flat gas loopers: 841227478` and `Time constraint for loopers: -2.5049e+14 ns`, then used ~70 GB of
memory. Diagnosis (code of AliceO2 `dev` as of 2026-08-22, i.e. jarvis15's `O2/daily-20260821`):

1. Since O2 2026-07, every generator derived from `o2::eventgen::Generator` adds **TPC loopers** (a machine-learned
   e±-looper background, `GenTPCLoopers`) in `Generator::finalizeEvent` whenever a collision context exists. That
   covers all O2DPG productions, ours included. O2DPG never switches it off; it only sets `GenTPCLoopers.colsys`
   (`o2dpg_sim_config.py:192`).
2. The geometrical protection added on 2026-07-24 (`51d8507`, `GenTPCLoopers.geomProtection = true`) drops candidate
   loopers outside the TPC gas. In versions before **2026-09-04** (`c575d91`, "Add protection against skipped loopers"),
   "all candidates dropped" is reported as a **failure** (`Failed to import loopers particles`).
3. A failure makes the primary server **retry** the event (`O2PrimaryServerDevice.cxx:219-227`). The loopers generator
   counts its calls (`mCurrentEvent++`) and uses that counter to index the collision list of the timeframe
   (`TPCLoopers.cxx:238`: `mInteractionTimeRecords[mCurrentEvent]`). Each retry pushes the counter one step ahead.
4. After k failures, the last k events of the timeframe read **past the end** of that list. A garbage time window
   gives a negative looper number, which becomes ~1e9 loopers once cast to `unsigned int`, so the run never finishes.
   Our run: 4 failures and 97 collisions, so event 94 used index 93 + 4 = 97, the first index past the end. The
   arithmetic matches the hang exactly.

**Not caused by our generator.** The looper count depends only on the collision-context timing, never on the event
content (Λ, nπ⁰ channel, photons from π⁰ → γγ: none of it enters). Which events fail depends only on where the ML
model puts looper vertices. Our pre-selection loop also lives inside `generateEvent()`, not in O2 triggers. That
matters: an O2 trigger would call `finalizeEvent` once per *rejected* event too, and desynchronise the counter in
every version. Standard productions with the same O2 version are exposed to the same bug.

**What we do** (`run/envcheck.sh: injpol_loopers_confkey`, used by `runInjectedPol.sh`; override with `INJPOL_LOOPERS`):

| installed O2 | setting | effect |
|---|---|---|
| before 2026-07-24 | default | no geometrical protection, so there is no failure path |
| 2026-07-24 … 2026-09-03 (jarvis15) | `GenTPCLoopers.geomProtection=false` | loopers **kept**; the protection, and with it the failure/retry chain, is off. Cost: some loopers start outside the TPC gas, i.e. the pre-July-2026 behaviour of the model |
| 2026-09-04 or later | default | fixed upstream (`c575d91`, plus the per-collision time windows of `c756dc9`/`1f4567d`, 2026-10-31) |
| `INJPOL_LOOPERS=veto` | `GenTPCLoopers.loopersVeto=true` | no looper background at all (last resort) |

The confKey actually used is written to `logs/confKey_used.txt`. **Best long-term fix: update O2/O2Physics** (and
O2DPG) to a version from 2026-10-31 or later. That is also what the official productions do: they run tagged O2
releases with these fixes rather than switching the feature off. `run/progress.py` raises an ALARM if it ever sees a
runaway looper request, and `run/killRun.sh` stops one run without touching other jobs.

## 5. Generator-level vs Geant4 decays: characterization plan

Our Λ skip whatever Geant4 would have done to them before the decay point: elastic/inelastic hadronic interactions
in the beam pipe, ITS and TPC inner field cage. Everything after the decay point is standard. Three QA steps,
from cheapest to most complete:

1. **Size of the effect, from a standard production.** For every generator-level Λ (|y| < 1) in a standard
   pp production (anchored LHC25 MC, or our own run with the same Pythia card and `-gen pythia8`), classify its
   fate from the Kine/McParticles daughters: (a) decayed by Geant4 (daughters with process 4); (b) interacted before
   decaying (daughters with hadronic process codes, or no decay daughters while inside the volume); (c) left the
   volume. The fraction (b), as a function of p_T and of the would-be decay radius, bounds the problem. Λ in (b)
   rarely give a reconstructable V0, so the bias on the V0 sample is smaller still.
2. **Method systematic (validation step 3).** Compare `configs/none.cfg` (our isotropic generator-level decays)
   with a standard production at equal statistics and settings: V0 yields vs p_T/η/φ, decay radius, cos PA,
   DCA between daughters, DCA to PV, invariant-mass peak, and the fake ring signal R itself. Any significant
   difference is a systematic of the injection method, and becomes an entry in the systematics table.
3. **If (1) or (2) is not negligible:** implement §3.2 (transport-level polarized decayer, with a small O2 patch)
   and repeat (2).

---

## 6. Simplifications and departures from defaults (complete list)

| # | Item | Choice | Standard practice / justification |
|---|---|---|---|
| 1 | Λ decay | generator level, polarized | Standard: Geant4, isotropic. Necessary for injection; characterized in §5 |
| 2 | Pre-selection | ≥ 1 Λ/Λ̄ with \|y\| < 1 (no trigger requirement by default, §7b) | Biases the sample (applied identically to every model); saves Geant4 time (as the O2DPG gap-triggered generators do) |
| 3 | `--skipModules ZDC` | yes | O2DPG convention (42 O2DPG configs/scripts). ZDC is ~112 m from the IP and does not affect the barrel |
| 4 | `--readoutDets` | **not restricted** | My earlier "speed" suggestion is withdrawn. Anchored productions set it to the detectors read out in the anchor run (`o2dpg_sim_workflow_anchored.py:750-751`); unanchored runs keep all |
| 5 | `Diamond.width[2]=6.` | yes | As in the O2DPG run scripts (unanchored); anchored runs take the vertex from CCDB |
| 6 | Interaction rate | 500 kHz (unanchored) | Anchored runs take it from the run |
| 7 | Anchoring | unanchored for the first tests | Goal: anchor to a 2025 pp run (§9) |
| 8 | Λ from Ξ/Ω | **unpolarized** (Geant4 decays them) | Realistic unpolarized feed-down; polarization transfer to be added later |
| 9 | Λ from Σ⁰ | polarized like primary Λ | Pythia decays Σ⁰ → Λγ, and the Λ counts as primary in ALICE. Σ⁰→Λ polarization transfer not modelled |
| 10 | BRs | PDG 2025: 64.1% pπ, 35.9% nπ⁰ | Geant4 uses 63.9/35.8%; negligible, but it enters the step-3 comparison |
| 11 | cτ | Pythia τ₀ = 78.9 mm | Same value as Geant4 |
| 12 | α parameters | 0.749 / −0.758 (pπ), 0.670 / −0.675 (nπ⁰) | Identical to `lambdaJetPolarizationIonsDerived.cxx:81-82` + PDG 2025 |
| 13 | Generator-level jets | SlowJet (fjcore) anti-kT | Same algorithm as FastJet in the producer (§7) |
| 14 | `--event-gen-mode` | **integrated** | O2DPG default is "separated", which drops our daughters (§4b) |
| 15 | Trigger requirement | **none** by default (`presel.requireLeadP = 0`) | Data-like dilution; about 2.5× more Geant4 events per triggered event (§7b) |
| 16 | Generator RNG | Pythia's own `rndm`, seeded by O2 from `-seed` per TF | Reproducible productions (could be randomised later if wanted) |
| 17 | TPC loopers | kept; `GenTPCLoopers.geomProtection=false` **only** on O2 versions 2026-07-24 … 2026-09-03 | Workaround for an O2 bug that hangs the run (§4c); best fix is a newer O2 |

**Your friend's setup** (`ALICE_PhotonReconstruction/GenProduction`), designed for topological Σ⁰ reconstruction, not
for this analysis:
- **Reused:** `--skipModules ZDC`, `Diamond.width[2]=6.`, 500 kHz, `-tt aod`, and the overall structure
  (external macro + ini + batches).
- **Not reused:**
  - `Random:setSeed = on` in the Pythia card: it overrides O2's seeding and makes runs non-reproducible;
  - `TRandom3(0)` inside the generator;
  - deleting `*_Hits*.root` inside the batch scripts. That is fine for disk space later, but we keep the Kine files for truth checks;
  - the particle-gun enrichment.
- He used no `--readoutDets` restriction.

---

## 7. Jets: generator side vs reconstructed side

- **Is FastJet available?** Yes. `alidist/o2.sh` lists `fastjet` as a hard dependency of O2. But no O2 or O2DPG code uses
  it: AliceO2 only has `find_package(FastJet)`, and GitHub code search finds no FastJet or SlowJet use in O2DPG.
  Loading the FastJet library from a Cling macro is possible but untested. We therefore use **Pythia8's `SlowJet`**,
  which delegates by default to **fjcore**, the FastJet core bundled inside Pythia. That gives the same anti-kT
  algorithm with no external dependency.
- **Definitions mirror the producer** (`lambdaJetPolarizationIons.cxx`, LambdaPolRing-v4):
  - candidates: final-state charged particles, p_T > 0.15 GeV/c, |η| < 0.9, pion-mass hypothesis;
  - leading particle: valid if p_T > 2 GeV/c;
  - jets: anti-kT R = 0.4, E-scheme, |η_jet| < 0.9 − R, no background subtraction.
  - **Thresholds are looser than in the analysis on purpose:** generator-level jet p_T > 8 GeV/c (charged-only, pp;
    the producer's reconstructed-jet default of 20 GeV/c would be far too tight at generator level), leading particle
    p_T > 2 GeV/c (4 GeV/c is the more jet-pure alternative, `axes.leadPtMin`). Generator-level cuts only pre-select
    and define the injection axis; the analysis thresholds are applied later, by FastJet on reconstructed tracks, so a
    loose generator cut avoids sculpting the reconstructed sample.
  At generator level Λ are neutral, so they are excluded, and their daughters do not exist yet. Charged Σ/Ξ/Ω
  are included, which is also what particle-level physical-primary selections do.
- **Generator-level vs reconstructed jets as a study:** the standard JE practice is to find **particle-level jets in
  O2Physics from `McParticles`** (`PWGJE/JetFinders/jetFinderMCPCharged.cxx`), with the same FastJet code as
  detector level. That is allowed and well trodden. For your producer: add a particle-level jet table to your
  datamodel, built with the same `jetsProcess` code on McParticles (physical-primary charged particles). The
  generator's own axis, when it computed one, is stored in the MC header (`injpol.jet.*`, `injpol.leadP.*` in `o2sim_MCHeader`) and should
  agree **exactly** with that O2Physics particle-level axis. That agreement is a good cross-check of both sides.
- **Jet finding in the generator is only done when needed** (`Config::needJets()`): for the ring model with the jet
  trigger, or `presel.requireJet = 1`. Leading-particle and contaminant-only configurations run no jet finder at all.
  All generator-vs-reconstructed jet studies are done in O2Physics on `McParticles` (above), not in the production.
- **Cost of the jet selection:** a fiducial 8 GeV/c charged jet still rejects a fraction of minimum-bias Pythia events
  in `configs/ring_jet_0p30.cfg`. That's cheap next to Geant4, but `standalonePythiaClosure` prints
  the acceptance and Pythia time per accepted event so we can decide. The alternative is a p̂_T-biased card
  (`HardQCD:all` + `PhaseSpace:pTHatMin`), with event weights in the analysis.

---

## 7b. Trigger t̂, pre-selection and signal dilution

**Notation** (as in the analysis notes): **t̂ = trigger direction** (leading jet or leading particle), and
**n̂ = (t̂ × p_Λ)/|t̂ × p_Λ|** = the ring direction. The injected ring polarization is P* = P₀ n̂.

**Default since 2026-11-05: no trigger requirement.** t̂ is the generator-level **leading charged particle**
(`model.trigger = leadP`, p_T > 2 GeV/c). The pre-selection keeps events with ≥ 1 Λ/Λ̄ (|y| < 1) **whether or not**
they have a trigger (`presel.requireLeadP = 0`). In events without a generator trigger the Λ get **P* = 0**.

**Why: dilution, as in data.** In data the reconstructed trigger is sometimes found in events with no genuine hard
scattering (fluctuations, combinatorics), and those events add R ≈ 0 to the average. A closure test that only
simulates events with a generator trigger would lack that dilution. Also, 2 GeV/c is the analysis threshold too, so
reconstructed triggers just above it can come from generator particles just below it. Requiring a generator trigger
would remove exactly those events.

**What it costs** (`run/runTests.sh`, step 0, jarvis15, 2000 accepted events each):

| configuration | Pythia events accepted | accepted events with a generator trigger | Pythia + afterburner per accepted event |
|---|---|---|---|
| `ring_leadP_0p30.cfg` (default, no trigger requirement) | 19.1% | 39.6% | 0.008 s |
| `ring_leadP_0p30_requireTrigger.cfg` | 7.8% | 100% | 0.020 s |

So the default transports about **2.5×** more events through Geant4 per injected-signal (triggered) event. Each
accepted event also carries fewer Λ (3.5 vs 4.3 per event). The requirement can be switched back on per
configuration (`presel.requireLeadP = 1`), e.g. for high-statistics closure checks of the signal itself.

**Why a trigger is still needed at all.** The leading particle (later the jet) defines n̂ for the injection, so the
model needs *a* trigger definition. What changed is only that events without one are no longer thrown away.

**Planned direction.** Use the **jet** (anti-kT, `model.trigger = jet`) as the default t̂, the same definition as the
FastJet jets of the producer, with the leading particle kept as a cross-check. Before that we must define carefully
**what a physical generator-level jet is**:
- which particles enter: physical-primary charged particles (as JE particle-level jets do), so weak-decay products
  and the Λ itself are excluded at generator level;
- the mass hypothesis: pion mass, mirroring reconstruction, or true masses;
- R, the p_T threshold (8 GeV/c for now, charged-only in pp) and the fiducial |η| < 0.9 − R;
- what to do with several jets: leading only, or all jets with each Λ attached to the closest;
- how generator and reconstructed jets are matched (ΔR) when quantifying t̂ smearing. That study belongs in O2Physics,
  on McParticles, not in the production.

Whatever is chosen, it must be identical in the generator (`InjectedPolPythia.h`) and in the O2Physics particle-level
jet table, and stated in the analysis note.

## 8. Afterburner physics

- **Decay kernel** (`InjectedPolKernel.h`). This is your `lambda_pol_toy_model.cxx` kernel, reviewed:
  - *Correct:* the inverse CDF c = (−1 + √((1−A)² + 4Au))/A of dN/dc ∝ 1 + A c (checked analytically and numerically);
    the random azimuth around P*, which makes the arbitrary transverse basis harmless; and the boost to the lab.
  - *Improved:* (i) rationalised form c = (A − 2 + 4u)/(1 + √(...)), exact at A = 0 and free of cancellation for
    small |A|; your version switches formula at |P| < 1e-9 but loses relative precision for small A·|P|. (ii) cos ξ*
    used directly, removing an `acos` plus `cos`/`sin` round trip per decay. (iii) both daughters returned,
    back-to-back in the rest frame. (iv) parent mass from its own 4-momentum, giving exact conservation. (v) no
    `TLorentzVector`/`TVector3` temporaries passed by value.
  - *Speed:* about 3 × 10⁶ decays plus estimators per second on one core (`tests/testKernel.cxx`), negligible next to Geant4.
    Your original was already efficient: inverse-CDF sampling, no rejection.
- **Channels:** Λ → pπ⁻ (α₋ = 0.749), Λ̄ → p̄π⁺ (α₊ = −0.758), Λ → nπ⁰ (α₀ = 0.670), Λ̄ → n̄π⁰ (ᾱ₀ = −0.675). The
  π⁰ and n are handed to Geant4 (π⁰ → γγ, as in standard productions).
- **Decay vertex:** proper length ~ exp(cτ), lab displacement βγ l p̂, lab time γ l.
- **Models** (`ModelConfig`; any combination; |P*| is capped at 1, and the number of capped cases is reported):

  | key | P* | notes |
  |---|---|---|
  | `model.ring` | P₀ n̂, n̂ = (t̂ × p_Λ)/\|t̂ × p_Λ\| | same cross-product order as the consumer task (`...Derived.cxx:2020`); t̂ = trigger: leadP (default) or jet (`model.trigger`) |
  | `model.global` | P_G (sin ψ, −cos ψ, 0) | random ψ_RP per event → AO2D `McCollisions.eventPlaneAngle` |
  | `model.longitudinal` | A sin(n(φ_Λ − ψ)) ẑ | n = 2 by default, n = 1 is the variant from your slides |
  | `model.helicity` | P_H p̂_Λ | fragmentation-like |
  | `model.antiSign` | multiplies P* of Λ̄ | +1: same sign as Λ (vorticity-like) |

  Axes come from lab vectors and are used as quantization axes in the rest frame. That is exactly what
  R = (3/α)⟨p̂*·n̂⟩ projects on, so a perfect detector returns R = P₀ (verified in `testKernel.cxx`).
  > **Reminder for the next model iteration:** please share your internal note **"19 - PolarizationSourcesPaperReview"**
  > so the model suite can cover all known theoretical sources.
- **Truth record:** `McParticles` has no polarization column. ψ_RP goes to `McCollisions.eventPlaneAngle` (the only
  header quantity propagated to the AO2D). The axes and model amplitudes go to `o2sim_MCHeader` (`injpol.*` keys).
  The models are deterministic given the Λ momentum and these event quantities, so a later O2Physics MC process can
  recompute the injected P* of every Λ.
- **Cheap 5-Λ test:** `configs/test_5lambdas_ring_leadP_0p30.cfg` (`presel.minLambdas = 5`). It selects
  high-multiplicity events, so it is a test of the machinery only. Injecting extra Λ "guns" into each event would be
  the other option, but those would not be correlated with the jets.

---

## 9. Anchoring to 2025 pp

- Anchored productions use `MC/run/ANCHOR/anchorMC.sh`. It needs the anchor run number, the data pass name,
  production tags (`ALIEN_JDL_LPM*` variables) and CCDB access. It takes the reconstruction settings of that async
  pass from O2DPG itself. **Our generator plugs in through `ALIEN_JDL_ANCHOR_SIM_OPTIONS="-gen external -ini <ini>"`**
  (default `-gen pythia8`), with nothing else to change. The template is `MC/run/ANCHOR/tests/test_anchor_2023_apass2_pp.sh`.
- `MC/run/ANCHOR/tests/test_anchor_cases.csv` has 2025 anchors for pO, OO and NeNe (LHC25ad/ae/af, apass1/apass2), but
  **no 2025 pp entry**. We need from the PWG-LF/DPG MC contacts (or the MonALISA/Hyperloop dataset pages): an LHC25 pp run
  number, its pass name (for your `LHC25_pass1_MinBias`-type data), and the matching O2PDPSuite async tag. We should
  also ask whether an anchored general-purpose 2025 pp MC already exists, since it would be the reference for step 3.

---

## 10. Validation plan (in order)

0. **Pythia only** (`run/runTests.sh`): truth-level R = injected P₀; pre-selection efficiency and cost.
   **Passed** on jarvis15 (2026-11-05), with and without the trigger requirement (§7b).
1. **o2-sim smoke test** (`run/smokeTestSim.sh` → `checks/checkKine.C`): B1-B5 bookkeeping, plus truth-level R in the Kine.
   **Passed** on jarvis15 (2026-11-05): B1-B5 OK for all 44 Λ (B5 after making its tolerances float-aware).
2. **Small full chain** (`run/runInjectedPol.sh` → `checkKine.C` per TF + `checks/checkAODFlags.C`): B1-B5 and A1-A4.
   Then run the strangeness builder with MC labels, and check that MC-matched Λ are found and flagged physical primary.
   First attempt (2026-11-05) found the separated event-generation problem (§4b); **re-run with
   `--event-gen-mode integrated` pending.**
3. **Method systematic** (§5.2): `none.cfg` vs a standard production.
4. **Injected productions:** ring scan in P₀, then ring + contaminants. Run them through `lambdaJetPolarizationIons` →
   `...Derived` as data. Then apply the correction methods of other polarization papers to the same samples.

---

## 11. Later

- O2Physics MC stage: MC tables in the `lambdaJetPolarizationIons` datamodel (particle-level jets, true Λ ↔ V0 links,
  injected P*), a producer `processMC`, and consumer processes that use them.
- Correction-factor families from several injected productions (General Framework note "16"), and a systematic from
  the choice of correction method.
- Ξ → Λπ (and Σ⁰ → Λγ) polarization transfer.
- Upstreaming to O2DPG (`MC/config/PWGLF/...`) through the PWG-LF MC contacts, if it proves useful.
