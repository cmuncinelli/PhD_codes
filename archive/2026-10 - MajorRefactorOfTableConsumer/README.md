## Consumer cleanup (lambdaJetPolarizationIonsDerived.cxx, October 2026)

The consumer went through a cleanup in a series of commits, before we added the new fixed-direction and event-plane analyses.

It consisted of removing old code and behaviors, and of fully superseding the kinematic-cut ring "families" that drove the flow of the code in previous iterations. The output folders were then reorganized around what each object is used for.

> Collaterally, this cleanup also reduced execution time by a factor of O(35%) (e.g., what ran in 100 s now runs in 65 s).

### Commits

In case you'd like to recover anything from the old behavior, these are the commits, in order:

1. [First name standardizations](https://github.com/cmuncinelli/O2Physics/commit/9b4f384bd183b25363619b0c11e5846c1c04b990)
2. [Duplicate removal and furthering the name conventions](https://github.com/cmuncinelli/O2Physics/commit/b1ba80eebe070f2ebabda76ed863aeeb23f1be93)
3. [COMPLETE removal of the ring families!](https://github.com/cmuncinelli/O2Physics/commit/931a249060156c75d333b2dc2282f05a6b8f59fc)
4. [Smaller changes: folder names and simpler code (some functions abandoned, smaller ones created)](https://github.com/cmuncinelli/O2Physics/commit/e8f2dd024b491c18f79cb663865e6bbd79df7b84)
5. [Small performance and organization changes](https://github.com/cmuncinelli/O2Physics/commit/69a28bed4b866593adac4edd4112369d2e35cec6)
6. [THE folder reorganization](https://github.com/cmuncinelli/O2Physics/commit/ab79d011683895539f57eb3fb36a58d3db84aa44)

The code as it was before any of these is the parent of the first commit (`9b4f384~1` in git).

Below you can find some more details on what changed in these commits.

### Species-agnostic objects

Objects that used to exist once per species (the AEE block in `FakePolDrivers/PhiAEE/`, `FakePolDrivers/hMassDecayGeom{Left,Right}`, `FakePolDrivers/PhiAEE/h2dMassVsPhiAEE`, `FakePolDrivers/CosHEE/ByMass/`) are now a single object without a species token. The population is chosen by the configuration, not by the object: BothHyperons, JustLambda and JustAntiLambda files, with matching variations. In the configs with random proxies the per-species files only agree statistically with the old objects.

### Kinematic-cut families removed

`Ring/` is the only family left (no extra cuts). Each old family is now a configuration of its own:

- Lambda pT and rapidity cuts: `analysisLevelCuts` (`v0MinPt`, `v0MaxPt`, `v0MaxRap`, with `doAnalysisLevelCuts`).
- Proxy acceptance: `analysisLevelCutsJets` (`maxLeadJetEta`, `maxSubJetEta`, `maxLeadPEta`, with `doAnalysisLevelCutsJets`, off by default), applied together with the minimum pT, before any V0 processing.

`IntegratedCuts/<Proxy>/` keeps its four bins (all Lambdas, Lambda cuts, proxy |eta| < 0.5, both), filled by hand with the old family definitions, so it can differ from the cuts of `analysisLevelCuts`.

### Configuration keys that changed

| Old key | New key |
|---|---|
| `familySwitches.*` | removed (O2 ignores the dead keys) |
| `minLeadParticlePt` | `analysisLevelCutsJets.minLeadPPt` |
| `minLeadJetPt` | `analysisLevelCutsJets.minLeadJetPt` |
| `minSubLeadJetPt` | `analysisLevelCutsJets.minSubJetPt` |

The minimum pT is always applied: `doAnalysisLevelCutsJets` gates only the |eta| acceptance.

### Output layout

Object names stayed (apart from the renames listed in the name map below), the folders changed. `Ring/` is proxy-first, so one proxy shows everything it has to offer. The QA folders are topic-first, so the proxies can be compared side by side.

| Folder | Holds |
|---|---|
| `SelectionFlow/` | the V0 selection cut flow and the Armenteros plots |
| `KinematicsQA/` | proxy kinematics (`<Proxy>/`) and V0 kinematics (`V0/`) |
| `IntegratedCuts/<Proxy>/` | the four-bin cut comparison, plus `h2dPairedNV0sVsCentrality` next to `pRingVsNV0s` (LeadJet) |
| `EventMixingQA/` | mixing QA (`AnglrCorrltns/` --> `AngularCorrelations/`) |
| `RzDiagnostics/` | longitudinal-ring diagnostics (only with `useRingZ` or `useRingZT`) |
| `RawPolarization/` | the polarization components (the old `Ring/QA/`) and `PolMaps/` |
| `Ring/<Proxy>/` | the ring observable: integrated objects at the root, then `AngularDependence/`, `PtDependence/`, `EtaDependence/`, `MassDependence/`, `EventLevelDependence/`, `RingMaps/`, `RingKernel/`, `KappaEff/`, `V0ProxyPairQA/` and `SigExtract3DInput/` (3D profiles with a mass axis) |
| `Ring/DeltaMethod/` | the leading-jet event tracker, read only by `extractDeltaErrors` |
| `FakePolDrivers/` | V0-related and proxy-related dependencies that greatly influence the full 3D ring signal: `CosHEE/`, `PhiAEE/`, `PhiStar/`, `CosHEEVsPhiStar/`, `ProxyEta/`, `ProxyZ/` (the old `HelicityEfficiencyQA/` and `EtaStudy/`) |

Rarely viewed folders come last: the registry writes its objects in booking order, so the booking order is the folder order in the file.

### Output changes

Every rename, move and removal, from the names before the cleanup to the final layout, is listed in [`nameMap_S0_to_C3b.txt`](nameMap_S0_to_C3b.txt), next to this README. One "oldPath newPath" pair per line, paths relative to the file root, with a newPath of "-" for an object that was removed on purpose. If an old macro or plot of yours reads an object that no longer exists, look its old path up there. The main changes:

- Proxy tokens: `LeadJet`, `LeadP`, `SubJet` everywhere, as folders inside `Ring/` and in the shared QA folders.
- Three renames in `V0ProxyPairQA/`: `hProxyPt`, `hMass`, `hLambdaPt` --> `hPairedProxyPt`, `hPairedV0Mass`, `hPairedV0Pt`.
- `ThreeBinMass/...Vs3BinMass` --> `FiveBinMass/...Vs5BinMass`: the CheapSigExtract mass regions now have five bins (left sideband, left gap, peak, right gap, right sideband), built from the NSigma configurables. The peak is bin 3. The   `CosHEE/ByMass/` profiles use the same five bins. `KappaEff/` keeps its two bins (sideband, peak).
- The 3D shadow counters (`h3d*`) and `hRingObservableNumMass` are gone: the TProfile3D bin entries hold the same counts.
- `h2dNV0sVsCentrality` --> `IntegratedCuts/LeadJet/h2dPairedNV0sVsCentrality`, filled once per collision with a valid leading jet.
- `h2dDeltaPhiVsProxyPhi`: the axes now follow the fill order (x = DeltaPhi, y = proxy phi).
- Two bug fixes: `smearSigma` is now computed in init (changes the forceJetDirectionSmudge configs with jetR != 0.4), and the LeadP t_z is taken after the proxy distortion, as for the leading jet (changes the LeadP `RingKernel/` objects in the configs with a distortion).
- The LeadP direction is built from the cached eta/phi, as for the jets (rounding-level change).

### How the changes were checked

Every step of the cleanup was checked against the outputs of the code before it with [`compareConsumerOutputs.cxx`](../../O2PhysicsTools/compareConsumerOutputs.cxx). It compares two consumer outputs object by object and bin by bin (profiles on their raw accumulators), after the name map, and sorts what it finds into content, structure and label differences, missing and extra objects:

```bash
root -l -b -q 'compareConsumerOutputs.cxx("before.root", "after.root", "nameMap_S0_to_C3b.txt", 1e-12)'
```

(or, an example of what I actually ran, with the proper paths in Jarvis15:)
```bash
root -l -b -q '/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs.cxx("/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/_archive/results_consumer-BeforeMajorCondensingAndSimplification/ConsumerResults_JustLambda.root", "/home/users/cicerodm/RingPol/LHC25ae_pass2/ITSandTPC/results_consumer/ConsumerResults_JustLambda.root", "/home/users/cicerodm/PhD_codes/O2PhysicsTools/compareConsumerOutputs-nameMap_S0_to_C3b-BeforeRefactorToAfterFolderChange.txt", 1e-12)' > compareConsumerLog.log
```

From the code before the cleanup to the final one, the only content differences were the expected ones: the LeadP objects (rounding of the LeadP direction, with entries migrating across bin edges), `h2dPairedNV0sVsCentrality`, the sums of the collapsed species objects in BothHyperons files and, in the configs they touch, the two bug fixes above. Everything else (LeadJet, SubJet and V0-only objects in particular) came out identical within 1e-12, up to titles. Configs with random proxies differ from run to run anyway, and were compared with a generous tolerance.