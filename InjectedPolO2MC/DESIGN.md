# InjectedPolO2MC -- archaeology findings and afterburner design

*First draft, 2026-11-04. Source versions inspected: AliceO2 `a6674c4` (2026-11-03), O2DPG (fork,
in sync with upstream) `e9d5f2d` (2026-11-03), O2Physics (fork) `8018361` (2026-10-13).
File/line references below are relative to those trees.*

---

## 1. TL;DR

- **Where to inject:** an O2DPG-style **external generator macro**, i.e. a subclass of
  `o2::eventgen::GeneratorPythia8` that overrides `importParticles()`. After Pythia fills the event,
  we decay every final-state Λ/Λ̄ ourselves with a polarized angular distribution, and we append the
  daughters to the particle list with a **displaced production vertex**. Geant4 then transports the
  daughters (and everything else) normally. No change to O2, O2DPG or O2Physics is needed.
- **Recipe already exists in O2:** `Generators/src/EventDecayer.cxx:303-362` (`EventDecayer::decayOne`)
  and `O2DPG/MC/config/PWGEM/external/generator/Generator_pythia8_forcedDecays.C` do this for
  isotropic Pythia decays. We copy their bookkeeping and replace the kinematics with the polarized
  kernel already in `vorticity_hydro/lambda_pol_toy_model.cxx` (`sample_P_angle_proton_from_cos_xi`
  and `Lambda_decay`).
- **Main pitfall:** if we copy that recipe verbatim, the **Λ loses `isPhysicalPrimary` and the p/π
  daughters gain it**. The daughters also get process code 0 instead of 4 (`kPDecay`), which
  breaks the strangeness builder's generated-V0 MC tables. Section 4 gives the fix: Λ keeps HepMC
  status 1 but is not transported, and the daughters get `UniqueID = kPDecay`.
- **Alternative rejected for now:** letting Geant4 hand Λ decays to a custom external decayer.
  O2DPG already routes Λ decays this way in one case, but polarization cannot reach the decayer
  without patching O2 core (section 3.2).

---

## 2. How an O2DPG MC production works (the part relevant to us)

```
o2dpg_sim_workflow.py  ──► workflow.json ──► o2_dpg_workflow_runner.py -tt aod
                                               │
   per timeframe:  o2-sim (generator + Geant4) │ digitization │ reco (ITS/TPC/TOF/...) │ secondary vertexing
                   │                                                        │ AOD producer ──► tfN/AO2D.root
                   └─ generator = Pythia8 or external macro (.C, JIT-compiled by ROOT/Cling)
```

- An external generator is selected with `-gen external -ini my.ini`, where the ini contains
  ```ini
  [GeneratorExternal]
  fileName=/abs/path/generator_pythia8_injectedPol.C
  funcName=generator_injectedPol(<args>)
  [GeneratorPythia8]
  config=/abs/path/pythia8.cfg
  ```
  Concrete local example: `gianniliveraro/ALICE_PhotonReconstruction/GenProduction/` (`configParticleGun.ini`,
  `generator_pythia8_gun.C`, `runSimulation.sh`). It uses exactly this mechanism plus batch scripts.
- **Event loop** (`AliceO2/Generators/src/Generator.cxx:308-384`): `generateEvent()` → `importParticles()`
  → `finalizeEvent()` → trigger (loop until accepted) → `decayParticles()` (optional `EventDecayer`)
  → `addTracks()`.
- **Units in `mParticles`:** whatever the generator uses. For `GeneratorPythia8` that means mm and mm/c.
  The base class converts with `mPositionUnit = 0.1` and `mTimeUnit = 3.3356410e-12`
  (`Generator.h:157-158`).
- **Vertices:** `PrimaryGenerator::AddTrack` adds the smeared interaction vertex on top of the
  particle's own production vertex (`PrimaryGenerator.cxx:148-150`). A displaced decay vertex set in
  the generator is therefore preserved.
- **Transport flag:** `ParticleStatus::kToBeDone` on each `TParticle` decides whether Geant4 tracks it.
  Particles with the flag off still go to the stack and the MC record (`Detectors/Base/src/Stack.cxx:228-247`).
- **Process code:** the generator's `TParticle::UniqueID` is stored as the MC process
  (`Stack.cxx:240`, `MCTrack.h:399`), and `MCTrack::isPrimary()` is
  `process == kPPrimary || no mothers` (`MCTrack.h:75`).
- **Pythia defaults already keep Λ undecayed** (`ParticleDecays:limitTau0 on`, `tau0Max = 10` mm/c;
  Λ has cτ = 78.9 mm), so in standard productions Geant4 decays Λ isotropically.

---

## 3. Injection points considered

### 3.1 Generator-level decay (chosen)

Override `importParticles()` (pattern used by about 30 macros in O2DPG, e.g. `PWGLF/pythia8/generator_pythia8_LF.C`,
`PWGEM/external/generator/Generator_pythia8_forcedDecays.C`):

```cpp
bool importParticles() override {
  if (!GeneratorPythia8::importParticles()) return false;
  computeEventAxes();      // jet proxy, (random) reaction plane, ...
  injectPolarizedDecays(); // decay every toBeDone Λ/Λ̄ with our kernel, append daughters
  return true;
}
```

Pros: the whole event (jet axis, multiplicity, a reaction plane) is available. It needs no core change,
runs in a single process, and is easy to make reproducible.
Cons: the Λ does not interact with material before its decay (absorption/scattering in the beam
pipe/ITS over a few cm is a small effect, but it should be listed as a systematic of the method),
and the decay length has to be sampled by us.

### 3.2 Decay during Geant4 transport through an external decayer (rejected for v1)

There is a precedent: `O2DPG/MC/run/PWGLF/run_DoubleLambdaTriggered.sh` with
`MC/config/PWGLF/pythia8/decayer/g4_ext_decayer_lambda.in` (`/mcPhysics/setExtDecayerSelection lambda anti_lambda`).
In that setup Geant4 hands Λ decays to `o2::eventgen::DecayerPythia8` (`G4.configMacroFile=...`).
A custom `TVirtualMCDecayer` could do polarized decays there, which would keep material effects and
all bookkeeping exactly as standard. However:

- `TVirtualMCDecayer::Decay(int pdg, TLorentzVector* p)` receives only the PDG code and momentum, with
  no event context;
- polarization is **hard-coded to zero** when primaries are pushed (`PrimaryGenerator.cxx:165-167`),
  so it cannot travel with the Λ;
- `o2-sim` transports chunks of the same event in **separate worker processes**, so a shared lookup
  table filled by the generator does not work either.

This would need an O2 core patch (propagate `TParticle` polarisation in `AddTrack`, read it in the
decayer from the current stack track). Keep it as a possible v2 if Λ-material effects turn out to matter.

### 3.3 Post-hoc reweighting of a standard production (not chosen)

We could reweight standard-MC Λ candidates by 1 + α P·p̂*. It's cheap, but it is not a real
"injected" sample. Reweighting does not change the topologies seen by the reconstruction. It is
still useful later as a fast cross-check of the full injection.

---

## 4. Bookkeeping: how to append the daughters (the critical part)

Relevant O2 code:

- `MCTrackNavigator::isPhysicalPrimary` (`DataFormats/simulation/src/MCUtils.cxx:73-...`):
  - returns false for any generator particle with HepMC status > 1;
  - for generator particles (`isPrimary()`), returns true if stable (Λ, p, π are all "stable" in
    `MCUtils.h:87-...`), with no check of the mother's identity apart from K0L;
  - for transport particles, returns false if the mother is a light-flavour hadron (Λ).
- AOD producer (`Detectors/AOD/src/AODMcProducerHelpers.cxx:229-235`): non-primaries get the
  `ProducedByTransport` flag, and their `statusCode` is replaced by the process code.
- O2Physics strangeness builder (`PWGLF/TableProducer/Strangeness/strangenessbuilder.cxx:2020`) only
  accepts V0 daughters with `getProcess() == 4` (`kPDecay`) when building generated-V0 MC info. It
  also uses the mother's `isPhysicalPrimary()` (`:1874`, `:1999`).

What `EventDecayer` / `forcedDecays` do, and what that would give us:

| | Λ (mother) | p, π daughters |
|---|---|---|
| Recipe as-is (HepMC 2 on mother, daughters `kPPrimary`) | **not** PhysicalPrimary | **PhysicalPrimary**, process 0, not `ProducedByTransport` |
| Standard production (Geant4 decay) | PhysicalPrimary | secondary, process 4, `ProducedByTransport` |
| **Proposed** | PhysicalPrimary | secondary, process 4, `ProducedByTransport` |

**Proposed encoding:**

```cpp
// mother Λ at index iL
encodeParticleStatusAndTracking(L, /*hepmc*/ 1, genStatus(L), /*toBeDone*/ false); // keep HepMC "final" -> stays PhysicalPrimary
L.SetFirstDaughter(iD1); L.SetLastDaughter(iD2);

// each daughter d
TParticle d(pdg, 0, iL, -1, -1, -1, px, py, pz, E,
            L.Vx() + dx, L.Vy() + dy, L.Vz() + dz, L.T() + dt);   // displaced decay vertex, generator units (mm, mm/c)
encodeParticleStatusAndTracking(d, 1, 91, /*toBeDone*/ true);
d.SetUniqueID(kPDecay);                // -> MCTrack::isPrimary() == false (has mother & process != kPPrimary)
d.SetWeight(L.GetWeight());
```

This must be **verified on a 10-event test**, as the first thing we do (see section 7), by dumping
`o2sim_Kine.root` and the AO2D `McParticles` flags. Things to watch:
- that the stack's index remapping keeps the daughter → Λ link (generator particles are kept
  1:1, `Stack.cxx:242`);
- that `PrimaryGenerator::AddTrack` does not complain (it only checks status encoding when
  `proc == kPPrimary`, `PrimaryGenerator.cxx:143`);
- whether any downstream code assumes "HepMC status 1 ⇒ transported" (none found so far).

---

## 5. Afterburner physics

**Decay kinematics:** reuse the kernel in `vorticity_hydro/lambda_pol_toy_model.cxx`:
- `sample_P_angle_proton_from_cos_xi` samples cos ξ* from dN/dcos ξ* ∝ 1 + α|P*| cos ξ* by
  inverse CDF;
- `Lambda_decay` builds p* in the rest frame and boosts it to the lab. The pion is p_π* = −p_p*.
  This needs adding (the toy only returns the proton).

**Constants:** these must match the analysis exactly (`lambdaJetPolarizationIonsDerived.cxx:81-82`):
α_Λ = **0.749**, α_Λ̄ = **−0.758**. The toy model uses 0.747, so use one shared header and no
hard-coded copies. A mismatch would be a guaranteed non-closure of about 0.3% plus a sign issue for Λ̄.

**Decay channels:** for the analysis only Λ → pπ⁻ (BR ≈ 63.9%) matters. Λ → nπ⁰ (≈ 35.8%) must
still be decayed, to keep yields right. Options are isotropic, or polarized with α₀ ≈ α. Take masses
and BRs from `TDatabasePDG` or Pythia's particle data, to stay consistent with the rest of the
simulation.

**Decay vertex:** sample the proper decay length from an exponential with cτ = 78.9 mm (taken from the
same particle table Geant4 uses). The lab displacement is βγ cτ p̂_Λ, with dt = γ τ. Daughters that start outside the
TPC are fine.

**Which Λ get polarized:** all generator-level, to-be-transported Λ/Λ̄. That includes Λ from Σ⁰
(decayed electromagnetically inside Pythia, and counted as primary in ALICE). Λ from Ξ/Ω are decayed
by Geant4 and stay **unpolarized**. That is a realistic, unpolarized feed-down contamination, and it should be
documented. Polarization transfer in Ξ → Λπ can be added later.

**Polarization models** (configurable, combinable, each with an amplitude):

| Model | P* direction (rest frame) | Notes |
|---|---|---|
| `none` | – | baseline, all estimators → 0 |
| `ring` | t̂ = (p̂_jet × p_Λ)/\|…\| | **same sign convention as the analysis** (`...Derived.cxx:2020-2021`) |
| `global` | −ŷ′ (normal to a random ψ_RP) | pp has no RP, so draw ψ_RP per event and record it |
| `longitudinal` | ẑ · A sin(2(φ − ψ₂)) | also the variant ẑ · A sin(φ − ψ_RP) from the talk |
| `helicity` | p̂_Λ | jet-fragmentation-like |

Frame convention: the axis is defined from lab vectors and used directly as the quantization axis in
the Λ rest frame. That is exactly what the measured R observable projects on. Add the option of a lab-frame
4-vector boost (`boost_polarization_to_rest_frame`) for vorticity-like models.

**Jet proxy at generator level:** start with the **leading charged final-state particle** (|η| < 0.9,
pT > threshold). This mirrors the leadP proxy of the analysis and is trivially Cling-compilable.
FastJet anti-kT can come later, if FastJet is available in the O2sim environment (check this; it is
not guaranteed there).

**Truth record:** `McParticles` has no polarization column. Write a small sidecar ROOT tree from the
generator, keyed by (run seed, event number, Λ index), holding the injected P*, the axes and the model
parameters. Because the models are deterministic functions of McParticles plus event-level random numbers,
the truth can also be recomputed in an O2Physics MC process.

**Event pre-selection (speed):** optionally loop `mPythia.next()` inside `generateEvent()` until the
event has a leading particle above threshold **and** at least one Λ/Λ̄ in acceptance. This saves
Geant4 time (the same idea as the gap-triggered generators in O2DPG). It biases the event sample, which
is fine for closure tests, but the bias must be stated and applied identically everywhere.

---

## 6. Production settings worth fixing early

- **Pin the software release.** The friend's setup uses `O2PDPSuite::MC-prod-2026-v1` / `daily-20260106-0000-1`.
  `EventDecayer` (2026) is fairly new, so check which features exist in the release we use.
- **Speed without changing physics:** restrict hit creation with `--readoutDets ITS,TPC,TOF,FT0,FV0,...`.
  This keeps the material budget (`o2dpg_sim_workflow.py:76`). **Do not** use `--skipModules` on detectors
  between the IP and the TOF, because it removes their geometry (`:114`). FT0 is needed for `sel8` (TVX).
  Also delete `*_Hits*.root` after digitization, as `ALICE_PhotonReconstruction/GenProduction/micro.sh` does.
- **Anchoring:** for a meaningful closure against LHC23 data, consider `o2dpg_sim_workflow_anchored.py` with
  the same anchor settings as LHC23k4i. That needs CCDB access from the cluster.
- **Seeds:** pass `-seed` per job, and seed the afterburner from the generator seed (not
  `TRandom3(0)`, as in the friend's macro) so that productions are reproducible.

---

## 7. Validation plan (in order)

0. **Pythia-only, no O2:** run the afterburner's kernel inside the existing `lambda_pol_PYTHIA.cxx`
   pipeline. Truth-level R must recover the injected P₀.
1. **o2-sim smoke test** (10 events, `-tt sgngen` or `o2-sim -g external`): in `o2sim_Kine.root`,
   check that the daughters exist with displaced vertices, the Λ is not transported, the mother/daughter
   links are right, and the process codes are right.
2. **Small full chain** (1 TF, about 50 events) → AO2D. Check the `McParticles` flags (Λ PhysicalPrimary;
   daughters `ProducedByTransport`, process 4), and check that the strangeness builder MC labels find the V0s.
3. **Method systematic:** zero-polarization generator-level decays vs. a standard production (Geant4
   decays) with the same settings and statistics. Compare V0 yields, decay radius, η/φ and DCA distributions,
   and R. Any difference is a systematic of the injection method itself.
4. **Injected productions:** `ring` scan in P₀, then `ring` + contaminants. Run them through `lambdaJetPolarizationIons`
   → `...Derived` exactly as data, and compare with the truth sidecar.

---

## 8. Open questions / later

- An MC process in `lambdaJetPolarizationIons` (currently commented out, `:2156`) for truth-level R and efficiency maps.
- Correction-factor families from several injected productions (General Framework note "16"). See the
  README and the paper roadmap; this is deliberately out of scope for v1.
- Upstreaming: if this works, propose the generator to O2DPG (`MC/config/PWGLF/...`) through the PWG-LF MC contacts.
