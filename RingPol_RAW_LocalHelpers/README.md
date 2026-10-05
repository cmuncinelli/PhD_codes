# RingPol RAW Local Helpers

Local analysis framework for the **ring polarization observable** in $\Lambda$--jet systems. It
operates on both raw AODs and Hyperloop-derived data from outside the Hyperloop environment,
complementing the O2Physics tasks that run inside it: `lambdaJetPolarizationIons` (producer) and
`lambdaJetPolarizationIonsDerived` (consumer).

Everything here is post-processing. The consumer produces `ConsumerResults_<SUFFIX>.root`; the macros
in this folder turn that into signal-extracted yields, uncertainty estimates, QA and comparison
canvases.

## Contents

1. [Layout](#layout)
2. [The post-processing chain](#the-post-processing-chain)
3. [Conventions](#conventions)
4. [The longitudinal ring, $R_z$ (`useRingZ`)](#the-longitudinal-ring-r_z-useringz)
5. [Macro reference](#macro-reference)
6. [Appendix -- AO2D bit forensics in detail](#appendix----ao2d-bit-forensics-in-detail)

---

# Layout

| Subfolder | Contents |
|---|---|
| [`Local_framework/`](Local_framework) | Running the O2 DPL workflow locally on raw AODs: download helpers, path generators, multi-threaded execution scripts. |
| [`DerivedDataHY/`](DerivedDataHY) | Downloading and analyzing Hyperloop-produced derived data. Holds the train registry and the coordinator script. See its [own README](DerivedDataHY/README.md). |
| [`JsonExamples/`](JsonExamples) | Example DPL workflow configs (`.json`): both-hyperon, $\Lambda$-only, $\bar\Lambda$-only, and full QA with permissive $p_{\rm T}$ cuts. |

The C++ macros live at the top level of this folder. All but one are driven by
[`DerivedDataHY/run_all_wagons.sh`](DerivedDataHY/run_all_wagons.sh); each has its own section under
[Macro reference](#macro-reference).

| File | One line |
|---|---|
| `signalExtractionRing.cxx` | Signal extraction of $\langle R\rangle_S$ from the $\Lambda$ invariant-mass spectrum, by sideband fit or by equal-width window counting |
| `extractDeltaErrors.cxx` | Delta-method uncertainties on $\langle R\rangle$, preserving the numerator--denominator covariance |
| `makeCumulativeDCAdauProfile.cxx` | Robustness of $\langle R\rangle$ against tightening DCA cuts (the AEE probe) |
| `auxiliaryPerConfigPlots.cxx` | Per-config derivative plots: polarization vector fields, ring-observable 2D maps |
| `auxiliarySummaryPlots.cxx` | Cross-config aggregation: systematics, MC and toy-model overlays |
| `zvtxBitForensics.cxx` | Bit-level and integrity QA of the raw derived AO2Ds, before any O2 abstraction |
| `signalExtractionRingTest.cxx` | Old development version of the signal extraction code |
| `mergeAODDerived.sh` | Merges AOD derived files prior to local analysis (this isn't really needed in any code. Just an early-version convenience) |

---

# The post-processing chain

`run_all_wagons.sh` runs seven steps for every (wagon, consumer config) pair enumerated from `train_registry.conf`. All macros are compiled ahead-of-time into native executables; none is run as an interpreted ROOT macro in production (`signalExtractionRing.cxx` absolutely cannot be ran in JIT-mode via root's CLI! It is just too heavy a code for that unoptimized execution)

| Step | Macro | Scope | Output folder |
|---|---|---|---|
| 1 | `runDerivedDataConsumer_HY.sh` | per config | `results_consumer/` |
| 2 | `extractDeltaErrors.cxx` | per config | `results_DeltaErr/` |
| 3 | `signalExtractionRing.cxx` | per config | `results_SigExtract/` |
| 4 | `makeCumulativeDCAdauProfile.cxx` | per config | `results_CumulativePlots/` |
| 5 | `auxiliaryPerConfigPlots.cxx` | per config | `results_AuxPerConfig/` |
| 6 | `auxiliarySummaryPlots.cxx` | per **wagon**, after all configs, choosing a specific subset of them | `results_consumer/` |
| 7 | `zvtxBitForensics.cxx` | per **wagon**, parallel, runs on AO2Ds themselves, but outside of O2Physics for (an honestly paranoid) QAing | `results_consumer/` |

Steps 2--5 are skipped for a pair whose step 1 failed. Step 7 reads the raw AO2Ds rather than any consumer output, so it also runs under `--post-process-only`.

---

# Conventions

> (this actually applies to most of the READMEs in this repo, but I realize I hadn't explicitly stated this anywhere, so here it goes!)

## Where documentation lives

**This README carries the physics.** Motivation, derivations, the reasoning behind a design choice, anything long enough to need equations (and to become nearly-ureadable in the scripts) -- all of it belongs here, where Markdown and MathJax (thankfully!) make it readable.

Each `.cxx` keeps a short header saying **what it is, how to invoke it, and what it needs to run**, plus a pointer back here. In-body comments are unaffected: they still explain most of the logic of the implementation, which is what a reader scrolling through the file actually wants, and are usually the most readable thing in it.

The rule exists because a sixty-line preamble is read once and skipped forever afterwards, while it pushes the real code below the fold on every subsequent visit. That is a really bad practice!

## Command-line signature

Every per-config macro takes the same two arguments, in essence:

```bash
./<macro>.exe <inputFilePath> <outputDirectory>
```

The output directory is explicit, not inferred from the input path, and is created recursively if absent. The output filename is derived from the input basename by substituting the prefix, e.g. `ConsumerResults_BothHyperons.root` becomes `ErrorPropagation_BothHyperons.root`.

The two per-wagon macros differ: `auxiliarySummaryPlots.exe` takes `<consumerDir> [mcRefDir] [ppRefDir] [toyModelPath] [cutFolder] [doIndividualComparisons] [sigExtractDir]`, and `zvtxBitForensics.exe` takes a manifest and an output file (see its [appendix](#appendix----ao2d-bit-forensics-in-detail)).

## Kinematic-cut families

The consumer books its ring-observable histograms under four sibling folders, one per selection
scenario:

| Folder | Selection |
|---|---|
| `Ring` | no kinematic cuts -- **always booked** |
| `RingKinematicCuts` | $\Lambda$ cuts: $p_{\rm T}^{\Lambda} \in [0.5,\,1.5]$, $\lvert y_{\Lambda}\rvert < 0.5$ |
| `JetKinematicCuts` | jet cuts: $\lvert \eta_{\rm jet}\rvert < 0.5$ |
| `JetAndLambdaKinematicCuts` | both |

Only `Ring` is unconditional. The other three sit behind consumer configurables and are simply
**absent from the output file** when switched off. That is a feature, not a fault: a config studying
one scenario should not pay for the other three.

Every macro that iterates these folders therefore does a **presence scan once, up front**, prints a
single summary line, and drives its sections from the resolved list:

```
 Families present : Ring, JetKinematicCuts
 Families absent  : RingKinematicCuts, JetAndLambdaKinematicCuts  (optional, not enabled in this config)
```

Two consequences worth knowing:

- A `WARNING` from inside a drawing or extraction section now always means something genuinely
  unexpected -- a family that *is* present but is missing a histogram it should contain. Absence of the family itself never reaches that code path.
- Nothing is written for an absent family. No empty `<Family>/` directory appears in the output, so a trimmed run and a complete one are distinguishable in a `TBrowser`.

`Ring` is marked **mandatory**: its absence exits `1`, because at that point the input is not the consumer output it claims to be.

### Inside a family

Every family has the same internal layout, organised by proxy:

| Folder | Holds |
|---|---|
| `LeadJet/`, `LeadP/`, `SubJet/` | everything that proxy's fill list fills, each with its own `QA/`, `EtaDependence/`, `ProxyPtDependence/`, `RingKernel/` (LeadJet, LeadP) and `KappaEff/` |
| `RingMaps/` | the ring planes, `p2dRingObservable[LeadP]Vs*` |
| `PolMaps/` | the polarization maps (proxy-agnostic, except `PrimeJet/`) |
| `DeltaMethod/` | the leading-jet event tracker, internal bookkeeping read only by `extractDeltaErrors` |

Object **names** carry their proxy as before (`LeadP/pRingObservableLeadPMass`); only the folders say it twice.

## Reporting: failures versus skips

`run_all_wagons.sh` prints two tables at the end, and they mean different things.

- **FAILURES** -- a step returned non-zero or crashed. Needs investigating; the log path is printed.
- **SKIPPED** -- a step had nothing to work on. The usual case is a wagon whose consumer has never
  been run, met under `--post-process-only`. Detected once per wagon rather than once per
  (config x step), which is what previously turned a single unprocessed wagon into a couple of dozen
  entries and buried the real failures underneath them.

Macro exit codes are binary on purpose: `0` for "ran", non-zero for "did not". A third code meaning
"ran, but some optional families were absent" was considered and rejected -- absent families are the
*common* case here, so encoding them as an exceptional status would make the exceptional path the
default one.

## Logging

Each step logs under its own output folder, alongside its results:

```
results_DeltaErr/logs/extractDeltaErr_<SUFFIX>.log
results_SigExtract/logs/sigExtract_<SUFFIX>.log
results_CumulativePlots/logs/cumulDCA_<SUFFIX>.log
results_AuxPerConfig/logs/auxPerConfig_<SUFFIX>.log
```

The exceptions are the consumer's own wrapper/batch logs and the two wagon-level logs
(`auxSummaryPlots.log`, `zvtxBitForensics.log`), which stay under `results_consumer/logs/`.

---

# The longitudinal ring, $R_z$ (`useRingZ`)

## What the switch does

The consumer configurable `useRingZ` (default off) redefines the ring observable, for every proxy, as its projection on the beam axis:

$$R_z \;=\; P_z\,n_z,\qquad \hat n=\frac{\hat t\times\vec p_\Lambda}{\lVert\hat t\times\vec p_\Lambda\rVert}.$$

Nothing else in the consumer changes. Every object filled with the ring -- the four cut families, `IntegratedCuts/`, `EtaStudy/`, the `HelicityEfficiencyQA/` AEE profiles and the Delta-method trackers -- then holds $\langle R_z\rangle$ under its usual name. The titles still read $\langle R\rangle$; the post-processing relabels them for $R_z$ files, and only where the meaning actually changes.

An $R_z$ output is identified **by its file name alone**: the consumer configs carry `_useRingZ` right after the family, as in `dpl-config-DerivedConsumer-BothHyperons_useRingZ_MixedEventProxies.json`, so every output follows as `ConsumerResults_BothHyperons_useRingZ_MixedEventProxies.root`. There is no marker object inside the file.

## Why it is cleaner

Write $\hat n=\cos\chi\,\hat\varphi+\sin\chi\,\hat\theta$ in the $\Lambda$ basis $(\hat p,\hat\theta,\hat\varphi)$, with $\chi$ the bearing of the proxy seen from the $\Lambda$ (the same $\chi$ as in the ring-geometry note). Because $\hat\varphi$ has no $z$ component, $n_z=-\sin\theta_\Lambda\sin\chi$, and the whole azimuthal group of the master equation -- the AEE fake, global polarization and $P_N$ -- drops out of $R_z$ **exactly**, candidate by candidate. Exactly, per candidate,

$$R_z=\sin^2\theta_\Lambda\,\bigl(\sin\chi\,P_\theta\bigr)\;-\;\sin\theta_\Lambda\cos\theta_\Lambda\,\sin\chi\,P_p .$$

What is left enters only with $\sin\chi$, which is odd between the two sides of the $\Lambda$ ("east" and "west"). The non-ring remainder therefore cancels up to the east--west imbalance $A_{EW}$, and $R_z$ has no large fake of its own. Two consequences follow:

- The denominator of the exact inversion is negligible, and the **linear subtraction** $R_{z,\rm true}=(R_{z,\rm meas}-R_{z,\rm fake})/\kappa_z$ is safe. Do not feed $R_z$ quantities into the full-ring Mobius formula: its denominator is written in terms of the full ring's fake.
- The response coefficient is $R_z$'s own, $\kappa_z$, not the full ring's $\kappa_{\rm eff}$ (see `KappaEff/` below).

The price is signal. $R_z$ keeps the ring with weight $n_z^2$, so near-side geometry retains about $0.45$ of the signal and $0.67$ of the significance.

## What is expected to vanish, and what is not

Several outputs of an $R_z$ run are **null tests**, and a non-zero value there is a finding, not a bug:

- **The ring kernel (`<Proxy>/RingKernel/`, and Sections 6--8 of `auxiliaryPerConfigPlots`).** The kernel exists to extract $B_\varphi$, which $R_z$ removes by construction. The profiles still fill and must be consistent with zero; the $B_\varphi$ and $M_0/M_1$ fits are meaningless and are read as nulls only.
- **$\langle R_z\rangle$ vs $\varphi_{\rm AEE}$.** It must be flat. Unlike the full ring, where binning in $\varphi_{\rm AEE}$ biases $\hat p^{*}\!\cdot\hat\varphi$ purely kinematically, $\hat p^{*}_z$ is independent of $\varphi^{*}_p$ for isotropic decays, so any structure is instrumental and E--W imbalanced.
- **`ringObservableOverJetZ`.** Built for the full ring's $\hat t_z$ invariance; not meaningful for $R_z$.
- **MixedEv** should give $R_z\approx0$. This is necessary but not sufficient: mixing cannot see fakes tied to the real jet's environment, and its $A_{EW}$ is not the data's.
- **pp** is zero only if vacuum fragmentation carries no ring. That is a physics question, not a closure test.

What survives in data besides the signal: the jet-$v_2$ times $v_2$-induced $P_z$ leak (entirely along $z$, so about twice as important relative to the signal as in the full ring), and any **field-odd** instrumental part. The solenoid field is the only thing in ALICE that breaks the mirror symmetry through the $(\hat z,\hat p_\Lambda)$ plane, which is what E--W cancellation relies on.

## East versus west, for free

Exactly, $n_z=\sin\Delta\varphi\,\sin\theta_t\sin\theta_\Lambda/\sin\Delta\theta$, so $\mathrm{sign}(n_z)=\mathrm{sign}(\sin\Delta\varphi)$, and zero is a bin edge of `axisDeltaPhi`. Hence:

- $A_{EW}$ comes from the halves of the existing `QA/hDeltaPhi*` counters;
- the east and west $\langle R_z\rangle$ come from integrating the halves of the existing `pRingObservable*DeltaPhi` profiles;
- the **odd part** of $\langle R_z\rangle(\Delta\varphi)$ is pure leak, and the even part is signal plus any field-odd residue.

No dedicated object is booked for any of this.

## `RzDiagnostics/` (booked only with `useRingZ`)

The full ring is computed alongside $R_z$ **only** in this mode, and only to fill this folder.

| Object | Content |
|---|---|
| `pRingVsChi<Proxy>`, `p2dRingVsChiVsMass<Proxy>` | $\langle R_z\rangle$ vs $\chi$, integrated and vs mass |
| `p2dRingVsDeltaPhiVsDeltaEtaLeadJet` | the same sectors in detector variables |
| `pRingPerpIntegrated`, `pRingPerp<Proxy>Vs{Mass,PhiAEE,DeltaPhi}` | the complement $R_\perp=R-R_z$, filled directly |
| `pCosSqThetaStarZLeadJetVsMass` | the standard, unweighted $\langle\cos^2\theta^{*}_z\rangle$ |

`<Proxy>` is `LeadJet`, `LeadP` or `SubJet`. The profile names read `Ring` rather than `RingZ` on purpose, since this whole folder only exists in $R_z$ mode.

**The $\chi$ sectors.** Under the reflections $(\chi\to-\chi,\ \chi\to\pi-\chi)$:

| Sector | Weight | What lives there |
|---|---|---|
| $(+,+)$ | $\sin^2\theta\sin^2\chi$ | the signal, plus field-odd instrumental parts only |
| $(-,+)$ | $\sin\chi$ | $\Lambda$-kinematic fakes ($B_\theta$, HEE) times $A_{EW}$ |
| $(-,-)$ | $\sin 2\chi$ | the in-plane "towards-the-jet" component $P_e$, $\hat e=\hat p\times\hat n$ |
| $(+,-)$ | -- | null (needs both the field and a $z$ asymmetry) |

$P_e$ and the helicity component $P_p$ are parity-odd, so strong production cannot generate them: both are purely instrumental. $P_e$ enters $R_z$ but **not** the full ring ($\hat e\perp\hat n$), and event mixing is blind to it.

**$R_\perp$ is filled, not subtracted.** $\langle R\rangle-\langle R_z\rangle$ has the right mean, but $R$ and $R_z$ share every candidate, so quadrature errors on that difference would be wrong.

## `KappaEff/` (always booked, in every cut family, one per proxy folder)

Both response coefficients share one form. With $u=R/c_\alpha$ the unitless ring ($c_\alpha$ the polarization prefactor) and $w$ its projection weight,

$$\kappa=3\,\frac{\langle u^2\rangle}{\langle w\rangle},\qquad w=\begin{cases}1 & \text{full ring: } u=\hat p^{*}\!\cdot\hat n\\[2pt] n_z^2 & R_z:\ u=n_z\,\hat p^{*}_z\end{cases}$$

so the same four profiles per proxy, all vs mass, serve both modes:

| Object | Filled with | Used for |
|---|---|---|
| `pKappaNum<Proxy>VsMass` | $u^2$ | $\kappa$ numerator |
| `pKappaDen<Proxy>VsMass` | $w$ | $\kappa$ denominator; the $R_z\to P_R$ conversion $R_z/\langle n_z^2\rangle$ |
| `pKappaNumTimesDen<Proxy>VsMass` | $u^2 w$ | $\mathrm{Cov}(u^2,w)$, for the error on $\kappa$ |
| `pRingTimesDen<Proxy>VsMass` | $R\,w$ | $\mathrm{Cov}(R,w)$, for the error on $R_z/\langle n_z^2\rangle$ |

In full-ring mode $w\equiv1$, so the last three are trivially redundant. They are filled anyway so the post-processing never has to know the mode.

The same moments are also filled **per CheapSigExtract mass region**, on a two-bin axis built from the consumer's own flags: bin 1 is the sideband (`v0InMassWindow`), bin 2 the peak (`v0InMassPeak`), and anything else goes to the underflow. These are `pRing<Proxy>VsMassRegion`, `pKappaNum<Proxy>VsMassRegion`, `pKappaDen<Proxy>VsMassRegion` and `pKappaNumTimesDen<Proxy>VsMassRegion`, all in `<Proxy>/KappaEff/`. They feed `Corrections/CheapSigExtract/` in the summary.

Three design points:

- **Unitless, on purpose.** $u$ carries no $\alpha$. $\langle R^2\rangle$ read off the existing profiles' spread would mix $\lvert\alpha_\Lambda\rvert\neq\lvert\alpha_{\bar\Lambda}\rvert$ in `BothHyperons`, and would not give the covariances.
- **Vs mass, per cut family.** $\kappa$ must be measured on exactly the candidates that make $R$, and the background's $\hat p^{*}$ distribution is not isotropic, so $\kappa_B\neq\kappa_S$ in general. The signal-region value comes out of the same sideband machinery as $\langle R\rangle_S$.
- **Ratios of means on one sample.** Each $\kappa$ and each $R_z/\langle n_z^2\rangle$ is a ratio of two means over the same candidates. First-order propagation needs only $\Sigma y$, $\Sigma y^2$ (which `TProfile` already stores) and the cross-product profile.

The full ring's $\kappa_{\rm eff}=3\langle(\hat p^{*}\!\cdot\hat n)^2\rangle$ is the data-driven second moment of the note, valid up to $O(\alpha P_{\rm true})$. $R_z$'s $\kappa_z=3\langle n_z^2(\hat p^{*}_z)^2\rangle/\langle n_z^2\rangle$ reduces to the standard $n_z^2$-weighted $3\langle\cos^2\theta^{*}_z\rangle$ for a $\phi^{*}$-independent, forward--backward-symmetric acceptance about $\hat z$; comparing it with the unweighted factor in `RzDiagnostics/` shows whether the weighting matters.

## Post-processing in $R_z$ mode

- **Pairing.** Data, MixedEv, pp and MC must be read in the same mode, which the suffix guarantees.
- **Toy Model.** Its overlay is a full-ring prediction, so it is dropped for $R_z$ files.
- **Titles.** Relabelled downstream from $\langle R\rangle$ to $\langle R_z\rangle$, only for $R_z$ files.
- **Every macro reads the mode from the file names**, so `run_all_wagons.sh` needs no change: the per-config macros see one file at a time, and the summary picks up any `_useRingZ` outputs in the wagon as extra families (`<Family>_RingZ/`). See [its section](#auxiliarysummaryplotscxx).

---

# Macro reference

## `signalExtractionRing.cxx`

**Step 3**, once per config. Output: `results_SigExtract/signalExtractionRing_<SUFFIX>.root`.
Takes the usual two positional arguments plus optional `--key=value` flags; `--help` lists them all
with their current defaults.

> Note to self: More documentation on this can be found on scribble "`10 - RingSignalExtractionSummary.pdf`" under my `Thesis/Summaries and LaTeX scribbles/` folder (local, not in this repo. Those are scribbles!)

### Why it exists

The polarization measured inside the $\Lambda$ mass peak is not exactly the signal polarization. The peak region contains combinatorial background, and that background is **not guaranteed to be unpolarized** (or to appear unpolarized) -- it can dilute the signal, be polarized differently, or carry a distorted angular structure of its own. What the peak gives is the mixture

$$\langle R\rangle_{\rm meas} = f_S\,\langle R\rangle_S + f_B\,\langle R\rangle_B ,
\qquad f_S + f_B = 1 .$$
(assuming a linear separation of the two, which seems good enough a separation when we see $\langle R\rangle_B\approx \text{cst}$, outside the mass peak window)

Signal extraction is therefore mandatory, not a refinement. Writing $N$ for yields in the peak window and solving for the signal,

$$\langle R\rangle_S = \frac{N_{\rm peak}\langle R\rangle_{\rm peak} - N_B\langle R\rangle_B}
{N_{\rm peak} - N_B} .$$

Both $N_B$ and $\langle R\rangle_B$ under the peak come from extrapolating a sideband fit inward.

### How

Per bin of the differential observable ($\Delta\phi$, $\Delta\theta$, the three jet-side $\eta$
proxies, the $\phi_\Lambda - \phi_p^*$ AEE probe, and 3D slices in $p_{\rm T}^{\Lambda}$ and
leading-jet $p_{\rm T}$):

1. Project the mass spectrum and fit the peak for $\mu$ and $\sigma$.
2. **Peak window**, either $[\mu - n_\sigma^{\rm peak}\sigma,\ \mu + n_\sigma^{\rm peak}\sigma]$ or,
   when `--massPeakHalfWidth` is positive, $[\mu - w,\ \mu + w]$ in absolute mass.
3. **Background**, by one of the two methods below.
4. Integrate over the peak window. In other words **bin-counting**, then background subtraction via
   an integrated background estimate -- much safer than fitting everything together, which is
   retained only as QA in `ResultsCombinedFit/`.

Every window is a command-line option; none is hardcoded, including in the plot labels. Run with
`--help` for the full list and the current defaults.

### Two background methods

`--bkgMethod` selects between them. They fill the **same booked histograms** and write to the same
paths, so the post-processing chain is unaffected; the method is recorded in bin 7 of
`hExtractionStatus` and tagged on the QA canvases, so a file always says which produced it.

**`sideband`** fits the sideband counts density with a `polN` and, separately, the sideband
$\langle R\rangle(m)$ with a low-order `polN`, then extrapolates both inward.

**`window`** (the default in `run_all_wagons.sh`) counts two windows of the **same width as the peak
window**, one on each side, separated by a gap. Nothing is extrapolated for $\langle R\rangle$.

The motivation is composition, not simplicity: equal-width windows sitting close to the peak share
its kinematic and topological makeup, so acceptance distortions -- the azimuthal and helicity
efficiency effects -- act on the background sample the same way they act on the signal sample. A
fitted extrapolation cannot promise that.

> **The band half-width is never configured.** It is taken from the *realized* peak half-width, after
> bin snapping, whatever defined the peak window. That makes the equal-width invariant hold by
> construction rather than by two configured numbers agreeing, and it is why window counting works
> with a sigma-defined peak as happily as with a mass-defined one. `massSidebandGap <= 0` likewise
> means "one peak half-width", so the whole geometry follows the peak.
>
> Because the band width depends on the realized window, the counting happens *inside*
> `ComputePeakWindowYields`, after the window has been resolved -- not at the call site.

#### What window counting does and does not replace

It gives the background's **polarization**. It does **not** give the background **yield** under the
peak, except in the linear case, and that distinction decides the design. Working in $t = m - \mu$
with a background $A + Bt + Ct^2$:

$$\int_{\rm peak} = 2Aw + \tfrac23 C w^3, \qquad
\int_{\rm bands} = 2Aw + \tfrac23 C w^3 + 2C g^2 w + 2C g w^2 .$$

The $A$ and $B$ terms cancel exactly -- symmetric equal-width bands are **exact for any linear
background**, which is the real content of the equal-width choice. The curvature does not cancel, and
the excess grows as $g^2$: already at $g = w$ the bands carry $4Cw^3$ against the peak's
$\tfrac23 Cw^3$, a factor of six too much curvature, and worse the further out the bands sit.

So the **yield always comes from the fitted counts density**, integrated over the peak window, in
both methods. Only $\langle R\rangle_B$ changes. Fit what a fit is good at -- a smooth yield -- and
count what counting is good at, a composition-matched polarization.

#### Pooling the two bands

From the sums, never by averaging two per-band means:

$$\langle R\rangle_B = \frac{\Sigma_L + \Sigma_R}{N_L + N_R},
\qquad \sigma^2_{\rm cond} = \sum_j (n_j\,{\rm SEM}_j)^2 .$$

This is exactly what a single wider bin would have given. Averaging the two means and adding their
errors in quadrature would both shift the central value and inflate the error whenever the bands have
unequal statistics -- on a realistic test with a 4:1 imbalance, by 3.5% and a factor 1.4 respectively.

The background numerator under the peak is then $M = \langle R\rangle_B \cdot B$ with conditional
error $B\,\sigma(\langle R\rangle_B)$, which is the same relation the fitted route arrives at and
exactly what the covariance structure below already assumes. The two methods therefore share every
derived quantity and every propagated difference: running both on one input is a genuine
cross-check, not a comparison of two implementations.

#### Geometry validation

Structural errors are **fatal**, checked once at startup before any extraction runs: a peak window
with no width, a negative gap, or an outermost band edge off the mass axis. Every one of those is a
typo rather than a property of the data.

Realized widths after bin snapping are **not** fatal. On a variable-width axis exact equality may be
unsatisfiable whatever is configured, and aborting there would leave the macro unrunnable for a
reason only a consumer rebin could fix. Those are reported and stored in `hAchievedWindow` instead.

### Fixed mass windows

Everything above is relative to the fitted `mu`, which is refitted in every angular bin. Six options
replace that with absolute boundaries in GeV/c^2:

```
--signalMassMin   --signalMassMax
--leftSidebandMin --leftSidebandMax
--rightSidebandMin --rightSidebandMax
```

All six or none; any one left unset falls back to the sigma-relative windows. When set they govern
the peak window, both sideband graphs, and the window-counting bands -- `mu` is not consulted for
any of them.

**Why this is not a convenience.** The peak is **not Gaussian in the tails**: the fitted Gaussian
falls faster than the real `dN/dm` does, so a window quoted in fitted sigmas understates how far the
peak actually reaches. Only around 6 sigma is one reliably clear of it. A sigma-relative window
therefore does not mean what it says, and it means something slightly different in every angular bin,
since `mu` and `sigma` move. Naming the intervals in mass removes both problems at once.

Two consequences follow:

- The windows are **reproducible** across bins, wagons and reruns.
- A bin is no longer **rejected** because its peak fit failed. With fixed windows the peak fit
  defines nothing, so it is demoted to QA and only the sideband background fit can still reject a
  bin. That removes most of the unused-angular-bin problem. The count of bins kept this way is
  reported, with the caveat that their `MassFits` canvases are not to be trusted even though their
  yields are.

A worked invocation, with `sigma = 0.0017127606` and `mu = 1.1156830`, giving a `mu +/- 1 sigma` peak
and `[6, 7] sigma` bands on each side:

```bash
SIGEXTRACT_OPTS=(
  --bkgMethod=window
  --signalMassMin=1.1139702   --signalMassMax=1.1173958
  --leftSidebandMin=1.1036937 --leftSidebandMax=1.1054064
  --rightSidebandMin=1.1259596 --rightSidebandMax=1.1276723
)
```

Startup prints both widths and notes any mismatch. Equal width is a property of the window-counting
method, not of the sideband fit, so a mismatch is **reported and never fatal** -- unequal bands may
well be wanted deliberately.

> **The band label follows the windows.** `SidebandBandLabel` once read the sigma fields
> unconditionally, so a run with fixed windows printed "8-10 sigma" in the log and in every saved
> canvas legend while integrating named mass intervals. A wrong-window run would have been
> indistinguishable from a right one in the output.

#### Narrow bands and the point count

A band 1 sigma wide holds about two bin centres on the current mass axis, so two bands give four
points. `minSidebandPoints` derives from the polynomial order, and a `polN` fitted with
`TF1::IntegralError`-grade covariance wants `N + 3`. Four against five rejects every bin, which is
what an over-tight configuration looks like: not a crash, but `NO RESULT` everywhere.

In **window mode** the `<R>` sideband is counted rather than fitted, so its point count constrains
nothing and only the counts fit sets the requirement: `minSidebandPoints` derives to `bkgPolOrder + 2`
there. The starvation warning prints the realized band, how many bin centres the axis offers inside
it, and the three remedies with their costs, because they are not interchangeable: lowering
`--bkgPolOrder`, widening the bands (costs the composition match with the peak), or refining the
consumer mass axis (costs counts per bin).

### The windows, and why the binning decides them

The mass axis is finite and variable-width, so a request for $\mu \pm n\sigma$ is snapped to bin
edges and comes back **asymmetric and wider than asked for**. The realized window is recorded per
bin in `<extraction>/Diagnostics/` and per extraction in `hAchievedWindow`, and every log line quotes
the achieved range next to the configured one. Plot labels carry the *configured* value; only the
integrated canvases annotate the achieved one, where a single number is well defined.

Two consequences worth internalising:

- **The Gaussian coverage of the window must follow it.** The signal yield is the fitted Gaussian
  integral times the fraction of it inside the realized window,
  $\tfrac12\left[\mathrm{erf}\!\left(\tfrac{x_{\rm hi}-\mu}{\sigma\sqrt2}\right) -
  \mathrm{erf}\!\left(\tfrac{x_{\rm lo}-\mu}{\sigma\sqrt2}\right)\right]$, which reduces to
  $\mathrm{erf}(n/\sqrt2)$ for a symmetric window. A hardcoded constant tied to one $n$ (the file
  once carried `0.9999366`, i.e. $\mathrm{erf}(4/\sqrt2)$) biases the yield by 4.7% at $n=2$.
- **A bounded sideband band can be unfittable.** On a coarse mass axis a $[3\sigma, 5\sigma]$ band may
  hold fewer points than the polynomial needs, in which case every bin is silently invalidated.
  `main` therefore resolves the band once at startup, against the reference spectrum, and refuses to
  run with a message naming the three ways out.

### Design constraints worth preserving

**The sideband fit uses `TGraphErrors`, not `TH1`.** The fit domain is *discontinuous* -- two
intervals with a gap between them. A `TH1` cannot express that: fitting over the full range would
include the peak, and a fit sub-range can only ever be one contiguous interval. A `TGraphErrors`
holding only the sideband points states the domain exactly.

**Counts are fitted as densities; $\langle R\rangle$ is not.** The mass axis is not uniformly binned,
so raw per-bin counts would hand the fit a spurious shape that follows the binning. Counts and
$\sum_i r_i$ are *extensive* -- they scale with bin width -- so dividing by the width is meaningful,
and integrating over the peak window recovers counts. $\langle R\rangle$ is *intensive*: it does not
scale with binning, and dividing it by a bin width would be meaningless.

**The background numerator comes from a fit to $\langle R\rangle(m)$, not to $\sum_i r_i$.** This is
a correction to an earlier design and the reasoning matters. In the sidebands
$\mathrm{d}(\sum_i r_i)/\mathrm{d}m = \langle R\rangle_B(m)\cdot \mathrm{d}N_B/\mathrm{d}m$. With the
counts modelled as a `pol2` and $\langle R\rangle_B$ roughly constant, that product is itself a
`pol2` -- so fitting $\sum_i r_i$ with a `pol1` asserted that a quadratic times a constant is linear,
and the fit absorbed the mismatch by forcing $\langle R\rangle_B(m) = {\rm pol1}(m)/{\rm pol2}(m)$, a
ratio with poles wherever the denominator crosses zero. Fitting $\langle R\rangle$ directly removes
the misspecification, makes `pol0` the honest default (it is exactly the assumption sideband
subtraction has always made), and puts that assumption on a QA canvas where it can be checked.

With the fit done on $\langle R\rangle$, the background numerator is written as the counts-weighted
average of $\langle R\rangle$ over the window,

$$R_B^{\rm eff} = \frac{\int R(m)\,b(m)\,\mathrm{d}m}{\int b(m)\,\mathrm{d}m} = \sum_k p_k w_k,
\qquad w_k = \frac{\int (m-m_0)^k\, b(m)\,\mathrm{d}m}{\int b(m)\,\mathrm{d}m},$$

so that $N_B\langle R\rangle_B = R_B^{\rm eff}\cdot N_B$ by construction. Because $R_B^{\rm eff}$ is a
*linear* functional of the fitted coefficients, its variance is exactly $w^{\rm T}\,\mathrm{Cov}(p)\,w$
-- no gradient integration, and therefore none of the `IntegralError` roundoff failures.

**Polynomials are fitted in $(m - m_0)$, not in $m$.** The $\{1, m, m^2\}$ basis is nearly degenerate
over a 70 MeV range centred on 1.115: coefficients come out around $10^{12}$ with alternating signs
and cancel to $\sim 10^{9}$ at the peak. Double precision absorbs the *value*, but the parameter
covariance does not survive it -- `TF1::IntegralError` was reporting *cannot reach tolerance because
of roundoff error* with integral errors of several thousand counts. Since $\mathrm{Var}(N_B)$ feeds
every downstream uncertainty, this was a correctness problem, not a cosmetic one.

**Background uncertainty comes from the full covariance, computed in closed form.** The polynomial
coefficients are strongly correlated, so adding their individual `GetParError` values in quadrature
would be wrong. `TF1::IntegralError` does use the full covariance, and it was used here -- but it
failed on this fit, twice over, with two GSL errors that are two symptoms of one cause:

```
Error 18  cannot reach tolerance because of roundoff error   (before recentring)
Error 11  number of iterations was insufficient              (after recentring, same window)
```

`IntegralError` builds its integrand from `TF1::GradientPar`, which differentiates **numerically**
with respect to each parameter. ROOT has no way to know the model is linear in its parameters, so it
takes finite differences whose step comes from the parameter magnitudes -- and with a density-scale
$p_0$ of order $10^9$ beside much smaller higher coefficients, that gradient is numerically noisy.
Adaptive quadrature then cannot converge on it, whatever basis is used. Recentring was necessary and
did help; it was never going to be sufficient.

None of that work is needed. A polynomial **is** linear in its parameters, so

$$\int_a^b {\rm pol}N = \sum_k p_k I_k, \qquad
I_k = \frac{(b-m_0)^{k+1} - (a-m_0)^{k+1}}{k+1},$$

is exact in closed form, and because the integral is a linear functional of the coefficients its
variance is exactly $I^{\rm T}\,{\rm Cov}(p)\,I$ -- full covariance, correlations included, which is
the only reason `IntegralError` was wanted. No quadrature, no gradient, no tolerance to fail to
reach. `ShiftedPolIntegral` does this, and there are now zero `IntegralError` calls in the file.

The same trick supplies the variance of $\langle R\rangle_B^{\rm eff}$ in the fitted method: it too
is a linear functional of the fitted coefficients, so $w^{\rm T}\,{\rm Cov}(p)\,w$ is exact.

**A `TProfile3D` must be converted before projecting.** `Project3D("yx")` on a `TProfile3D` returns a
`TProfile2D` whose every cell holds the *grand mean* over all $z$ bins in range -- not what a sliced
analysis wants, and with errors that are wrong for this purpose. `BuildNumFromProfile3D` therefore
converts to a `TH3D` first, in which

$$\text{bin content} = \bar R \cdot N = \sum_i r_i, \qquad
\text{bin error} = \frac{\sigma_R}{\sqrt N}\cdot N = \sigma_R\sqrt{N} ,$$

the absolute sum and its correct error, which *are* additive under projection. Do not reintroduce a
direct `Project3D` on the profile.

**A `TProfile2D` carries its own denominator.** `GetBinEntries` counts exactly the candidates that
went into the profile, so `CountsFromProfile2D` pairs a numerator and a denominator that cannot
describe different selections. That is what lets the AEE probe run with no separately booked counts
histogram, and it is the safer pairing even where one exists.

### Error propagation

All derived quantities come from four primitives measured in the peak window: the counts $T$, the
ring numerator $N = \sum_i r_i$, and their sideband-extrapolated backgrounds $B$ and $M$. The
unconditional covariance, in that order, is

$$\mathrm{Var}(N) = \sigma_{N|T}^2 + \langle R\rangle_{\rm peak}^2\,\mathrm{Var}(T),
\qquad \mathrm{Cov}(T,N) = \langle R\rangle_{\rm peak}\,\mathrm{Var}(T),$$

and likewise for $(B, M)$ with $\langle R\rangle_B$, the peak block independent of the sideband block
because they are measured in disjoint mass regions. The extra $R^2\mathrm{Var}$ terms convert the
*conditional* `TProfile` variances -- the spread of $\sum_i r_i$ at fixed candidate count -- into
unconditional ones. Dropping them gives an answer that looks plausible and is wrong.

`FinalizeDerivedQuantities` computes everything from those primitives in one place, which is also
what lets the **angle-combined** results reuse it: independent angular bins have independent
primitives, so summing the primitives and running the same function on the totals is algebraically
identical to a single extraction over their union.

**Differences between the three ring observables must not be added in quadrature.** They are three
functions of the same four primitives:

$$\mathrm{Var}(\langle R\rangle_{\rm peak} - \langle R\rangle_S) =
\left(\tfrac1T - \tfrac1S\right)^2\sigma_N^2 +
\frac{\sigma_M^2 + (\langle R\rangle_{\rm peak}-\langle R\rangle_S)^2\mathrm{Var}(T)
+ (\langle R\rangle_B-\langle R\rangle_S)^2\mathrm{Var}(B)}{S^2},$$

$$\mathrm{Var}(\langle R\rangle_{\rm peak} - \langle R\rangle_B) =
\frac{\sigma_N^2}{T^2} + \frac{\sigma_M^2}{B^2}, \qquad$$

$$\mathrm{Var}(\langle R\rangle_S - \langle R\rangle_B) =
\frac{\sigma_N^2 + (\langle R\rangle_{\rm peak}-\langle R\rangle_S)^2\mathrm{Var}(T)
+ (\langle R\rangle_S-\langle R\rangle_B)^2\mathrm{Var}(B)}{S^2}
+ \left(\tfrac1S + \tfrac1B\right)^2\sigma_M^2,$$

with $S = T - B$. On a representative primitive set, naive quadrature **overstates** the first by
about 60% and **understates** the third by about 20%. The second comes out as exact quadrature,
because those two are measured in disjoint mass regions -- a free consistency check on the whole
construction. All three are therefore computed here and stored, never re-derived downstream.

### Angle-integrated results: combine, do not project

Two ways to reduce a set of per-bin extractions to one number, and their **difference is itself a
measurement**:

- **Signal-weighted**: sum the primitives across bins and extract from the totals, giving
  $\sum_i S_i R_i / \sum_i S_i$. This is the *acceptance-weighted* answer -- the weights $S_i$ are
  exactly the non-uniform occupancies the azimuthal efficiency effect produces.
- **Flat-acceptance**: the unweighted bin average $\tfrac1n\sum_i R_i$. Bins hold disjoint candidates
  and are equal width, so this is a plain independent average and its error really is quadrature.

Their difference, stored as `hAEE_*`, is the azimuthal efficiency effect. Neither number alone shows
it. Writing $D = \sum_i (w_i - 1/n) R_i$ with $w_i = S_i/\sum_j S_j$, the $R_i$ are independent and
$\mathrm{Var}(D) = \sum_i (w_i - 1/n)^2\mathrm{Var}(R_i)$; the weights are treated as fixed, which is
safe while the bin yields are far better determined than the $R_i$ themselves.

For the $\phi_\Lambda - \phi_p^*$ extractions the projection route is **disabled outright**
(`IntegralMode::CombinePerBin`). Projecting every angular bin onto the mass axis before extracting
undoes the entire reason for splitting on a variable across which $\langle R\rangle$ changes sign: the
summed numerator sideband becomes a mixture of opposite-sign contributions, and a low-order
polynomial has no business describing it.

### The angle-combined output, and one label worth not misreading

`IntegratedCombined/` holds three labelled histograms rather than the fourteen single-bin ones it
started as. A `TBrowser` lists single-bin objects in whatever order it likes, and the grouping had to
be reconstructed from names.

| Object | Holds |
|---|---|
| `hCombinedSummary` | $\langle R\rangle_{\rm meas}$, $\langle R\rangle_S$, $\langle R\rangle_B$, $\langle R\rangle_S - \langle R\rangle_B$ |
| `hCombinedQA` | the other two differences, and the four "weighted-flat" rows |
| `hCombinedQuality` | purity, significance, signal and background counts, bins used against total |

They are split by **scale**, not by importance. Purity is of order $0.1$ and significance of order
$100$, while every $\langle R\rangle$ here is of order $10^{-2}$: on a shared axis the ring
observables collapse onto zero and nothing is readable.

> **"weighted-flat" is not a different $\langle R\rangle_S$.** It is the same quantity minus its
> flat-acceptance average -- what the number would be if the per-bin counts were ignored when
> integrating. The gap between the two is the acceptance weighting, which is what makes it an AEE
> probe. These rows were once labelled `AEE: <R>_S`, which read as "the $\langle R\rangle_S$ of the
> AEE observable". It is not that.

**What to quote.** $\langle R\rangle_S$ is the observable. $\langle R\rangle_S - \langle R\rangle_B$
is, in principle, a **diagnostic**: its error is propagated correctly, but a contrast is not a polarization. It is
worth watching because $\langle R\rangle_S \approx \langle R\rangle_B$ is the signature of a "signal"
that is really background leaking through. The remaining differences are cross-checks on the
subtraction, and the weighted-flat rows measure acceptance distortion, not polarization.

> A small warning though: if we assume that the background carries the same kind of fake polarization as the signal, but carries no true polarization because it is not composed of true Lambdas (and we assume that the proton daughters of the $\Lambda$'s are not polarized either), then $\langle R\rangle_S - \langle R\rangle_B$ becomes the result we actually want to report!!! After all, that would be equivalent to estimating the fake polarization from combinatorial V0s and then using that to correct the true V0s' polarization measurement! (Similar to Joseph Adams' method from his thesis' section "3.4 Polarization observable")

The flat-acceptance averages themselves live in `IntegratedCombined/NotForPhysics_FlatAcceptance/`,
with the warning in the histogram title so it travels with the object. The folder is deliberately not
called `QA/`: that name is used elsewhere for things that *are* trustworthy for what they claim,
whereas a flat average weights a thin bin exactly like a well-populated one. It answers "what would
this be with uniform acceptance", not "what is it", and it exists to be subtracted.

### Subtract first, then combine
> THIS SECTION IS PARTICULARLY IMPORTANT!

When the proxy-`eta` split is used, each half is extracted separately and the halves are recombined
afterwards. The **order** of subtraction and recombination is not a matter of taste.

Write the assumption explicitly. Within half `i`, the fake polarization `F_i` -- acceptance, AEE, HEE -- is common to signal and background:

$$R_S^i = T_i + F_i \quad (\text{true} + \text{fake}), \qquad R_B^i = F_i \quad (\text{fake only}),$$

so the per-half difference is $D_i = R_S^i - R_B^i = T_i$. **The fake part cancels exactly, inside
each half, with no weighting at all.** A combinatorial pair has no physics polarization; a real
Lambda carries the same detector distortions as one. That is what the subtraction is for.

**Method A, subtract then combine** gives the signal-yield-weighted true polarization and nothing
else:

$$D_A = \frac{S_P T_P + S_N T_N}{S_P + S_N}.$$

**Method B, combine each row then subtract** gives

$$D_B = D_A + \left( \langle F\rangle_S - \langle F\rangle_B \right), \qquad
\text{residual} = f\left[\frac{S_P - S_N}{S_P + S_N} - \frac{B_P - B_N}{B_P + B_N}\right].$$

That residual vanishes only if the halves are perfectly symmetric, **or** the purity is identical in both halves (If it were identical, then we could simply use $N^\Lambda_{\eta_{proxy}>0}$ and $N^\Lambda_{\eta_{proxy}<0}$ as weights instead of $S_P$ and $S_N$). Neither is guaranteed -- different acceptance on each side is precisely why the split exists. And it is multiplied by `f`, the **large** number. On a toy with a true polarization of 0.002 and a fake of 0.030:

| halves | A | B | B - A |
|---|---|---|---|
| symmetric | +0.002000 | +0.002000 | 0 |
| signal 5% asymmetric, background symmetric | +0.002000 | +0.003500 | **+75% of the true value** |
| signal 5% asym., background 8% asym. | +0.002000 | +0.001100 | **-45% of the true value** |

A 5% yield asymmetry between the halves corrupts the answer by most of its own size. That is the "summing two big numbers and getting zero" failure re-entering at the **recombination** step, after it had been successfully kept out of the fit. Method B undoes the benefit of splitting at the very last stage.

So: **subtract within a half, unweighted; combine across halves, signal-yield weighted, quadrature.**

There is an error-side argument too -- `D_B` would need $\mathrm{Cov}(R_S^{\rm comb}, R_B^{\rm comb})$, which is non-zero and is not stored, while `D_A` needs no covariance at any stage. Across halves the candidates are disjoint, so quadrature there is exact.

> This rests on `F` being common to signal and background *within* a half. If the fake polarization acted differently on a real Lambda's decay kinematics than on a combinatorial pair's, neither method would be clean. The `<R>_meas - <R>_B` comparison is one handle on that, and the unsplit-versus-recombined gap another.

> Notice that we make a hypothesis here, but if we assume that all the "fake" polarization is caused by the analysis cuts (DCAtoPV of the daughters, $p_T$ cuts, and all the cuts that were available in the toy model) and efficiency cuts, and assume also that these effects are the same for all V0 structures (be they $\Lambda$'s or just combinatorics), then it is a pretty reasonable hypothesis!

#### Per-row weights, and the closure test

Each row of `hCombinedSummary` is a mean over a **different** population, so each is recombined with
its own weight. Weighting them all by `S` is wrong for most of them:

| row | mean over | weight |
|---|---|---|
| `<R>_meas^FullMassRange` | every candidate, all mass, every angular bin | full-range counts |
| `<R>_meas^PeakWin,AccBins` | peak window, accepted bins only | `T = S + B` |
| `<R>_S` | signal candidates | `S` |
| `<R>_B` | background candidates | `B` |
| `<R>_S - <R>_B` | the true polarization carried by the signal | `S` |

> The first row is the **closure test**. It involves no fit, no mass window and no bin rejection, so recombining the two halves must reproduce `NoEtaSplit` exactly; if it does not, the weights are wrong. The second row need not close (with the first row), and **the gap between rows one and two is precisely what it cost to restrict the mass and to drop the bins that failed extraction**.

### Configuration: one set of numbers, four workflows

Four places perform the same measurement -- per angular bin, angle-integrated, denominator QA, and
the selection cut flow. They used to carry four independently tuned presets, which had drifted apart
on things that are not workflow-specific at all: one fitted raw counts while the others fitted
densities, one kept empty sideband bins, and the sigma limits disagreed. They now share **one**
`SidebandConfig` and differ only in `fallback`, the policy for a mass fit that does not converge,
which is genuinely semantic rather than numerical.

**Failure must be loud.** A zero-filled single-bin histogram is visually indistinguishable from a
genuine measurement of zero. Value histograms are therefore written *only* on success, while
`hExtractionStatus` is written unconditionally, so "the extraction failed" is always distinguishable
from "the extraction never ran". Every failure branch names itself in the log.

### The KappaEff moments

The response coefficient is $\kappa=3\langle u^2\rangle/\langle w\rangle$ (see [The longitudinal ring](#the-longitudinal-ring-r_z-useringz)), so its signal-region value needs only $\langle u^2\rangle_S$ and $\langle w\rangle_S$. Each is a per-candidate quantity profiled against mass, exactly like $R$, so each goes through **the same integrated extraction** -- same peak, same windows, same sideband model, same error propagation -- with nothing written for it but a folder name and a symbol for the labels:

```
<Family>/KappaEff/<Proxy>_Num/hIntegratedRSig_<Proxy>_Num   <-- <u^2>_S
<Family>/KappaEff/<Proxy>_Den/hIntegratedRSig_<Proxy>_Den   <-- <w>_S, R_z files only
```

The object names are the usual ones, so `hIntegratedRSig_<name>` is $\langle y\rangle_S$ for whatever $y$ was profiled; the titles carry the right symbol. The ratio is taken by the summary, in `Corrections/SignalExtracted/`.

$\langle w\rangle_S$ is extracted for `_useRingZ` files only, the mode being read from the file name as everywhere downstream of the consumer. For the full ring $w\equiv1$, so $\langle w\rangle_S=1$ identically and there is nothing to extract; a zero-spread profile would also leave the sideband fit with no errors to weigh. A consumer output that predates `KappaEff/` is reported in one line and skipped.

### Assumptions this rests on

1. Background polarization varies smoothly in invariant mass.
2. No strong mass--observable correlation exists (**must be validated**, not assumed).
3. Sidebands are representative of the background under the peak.

If these fail, the escape route is a simultaneous mass--polarization fit,

$$\text{Numerator}(m) = S(m)\,R_S + B(m)\,R_B, \qquad \text{Denominator}(m) = S(m) + B(m),$$

which avoids sideband extrapolation entirely. The `ResultsCombinedFit/` output directory is
groundwork for exactly that. It is **QA only** and must not be read as a result: the joint fit assumes
a *constant* $\langle R\rangle_B$ across the mass window while the sideband method fits its mass
dependence, so a disagreement between the two may be nothing but that modelling difference.

### V0 selection cut flow

Run once per config, from `AnalysisResults_merged.root` (one directory *above* the consumer output),
into `V0SelectionCutFlow/`. The producer's flow counter is laid out as (and this may need updating if we ever change the counter in the TableProducer, but here goes the current state FYI):

| ROOT bins | Block |
|---|---|
| 1--31 | generic V0 cuts, shared by both mass hypotheses |
| 32--42 | $\Lambda$-specific cuts |
| 43--53 | $\bar\Lambda$-specific cuts |

> **The two hypothesis blocks are not nested in one another.** A candidate goes down one branch or
> the other; only the generic block is a common ancestor. Chaining across the branch boundary would
> break the nesting assumption the retention uncertainties depend on, and reading one hypothesis
> block out of the other species' mass histogram means looking at candidates selected under the wrong
> hypothesis. Each species is processed along its own chain: 1--31 then 32--42, or 1--31 then 43--53.

`kNGenericV0Cuts` and `kNHypothesisCuts` at the top of the file **must** be kept in sync with
`V0SelectionFlowCounter` in the producer task.

---

## `extractDeltaErrors.cxx`

**Step 2**, once per config. Output: `results_DeltaErr/ErrorPropagation_<SUFFIX>.root`.

### Why it exists

The ring observable is a **ratio of sums** over candidates,

$$R = \frac{\sum_i r_i}{\sum_i n_i} \equiv \frac{S_r}{S_n} ,$$

and numerator and denominator are built from the *same* candidates, so they can be strongly correlated.
A standard error of the mean treats them as if they were not, dropping precisely the term that
matters. The consumer stores the five raw accumulators $S_r$, $S_n$, $\sum r_i^2$, $\sum n_i^2$ and
$\sum r_i n_i$ so that the correlation can be restored here, offline.

### The formula, and why the uncentred moments are exact

The delta method for a ratio gives

$$\operatorname{Var}(R) \simeq \frac{1}{S_n^{2}}
\Big[\operatorname{Var}(S_r) + R^{2}\operatorname{Var}(S_n) - 2R\operatorname{Cov}(S_r, S_n)\Big] .$$

`ComputeDeltaError` implements

$$\operatorname{Var}(R) = \frac{1}{S_n^{2}}
\Big[\textstyle\sum_i r_i^{2} + R^{2}\sum_i n_i^{2} - 2R\sum_i r_i n_i\Big],$$

which at first sight looks like it forgot to centre the second moments. It did not. For $N$ iid
candidates the centred quantities are

$$\operatorname{Var}(S_r) = \sum_i r_i^2 - \frac{S_r^2}{N}, \qquad
\operatorname{Var}(S_n) = \sum_i n_i^2 - \frac{S_n^2}{N}, \qquad
\operatorname{Cov}(S_r,S_n) = \sum_i r_i n_i - \frac{S_r S_n}{N},$$

so the terms the implementation omits sum to

$$-\frac{1}{N S_n^{2}}\Big[S_r^{2} + R^{2}S_n^{2} - 2R\,S_r S_n\Big]
= -\frac{\big(S_r - R\,S_n\big)^{2}}{N S_n^{2}} = 0 ,$$

because $S_r - R S_n = 0$ identically, by the definition of $R$. **The uncentred form is not an
approximation of the centred one; the two are algebraically the same number.** Checked numerically on
correlated toy data: the two expressions agree to machine precision, and both land within about 1.5%
of a bootstrap of the ratio.

The only surviving approximation is the delta method itself -- first order in the fluctuations of
$S_r$ and $S_n$ -- which is more than adequate at the statistics involved here.

### What it actually bought

Honestly: **not much.** The delta-method errors come out close to the SEM ones on real data. The
output keeps both side by side (`DeltaMethod_err/` and `SEM_method_err/` per family, plus
`pRingCuts_Delta` and `pRingCuts_SEM` at the top level) precisely so that this can be re-checked
rather than taken on trust. The value is in knowing that the agreement is a result and not an
accident.

Bootstrapping and jackknifing are the natural next additions to this file, but I wouldn't expect too much from them right now.

### Note on the summary histogram (and on the new behavior of error-reporting for the missing consumer families)

`pRingCuts_Delta` has four fixed, labelled bins, one per kinematic-cut family, and bin $i+1$ always means the same family whether or not its neighbours were produced. An absent family leaves its bin empty rather than shifting the others.

---

## `makeCumulativeDCAdauProfile.cxx`

**Step 4**, once per config. Output: `results_CumulativePlots/CumulativeProfiles_<SUFFIX>.root`.

### Why it exists

This is a robustness probe for the **Azimuthal Emission Efficiency (AEE)**: the daughter DCA cuts, combined with the magnetic field, make the reconstruction efficiency depend on $\phi^{*}$, which manufactures a fake polarization signal. If AEE drives what we measure, tightening the minimum DCA cut should move $\langle R\rangle$ systematically. If $\langle R\rangle$ is flat against that cut, the observable is robust against the dominant AEE source.

The consumer books $\langle R\rangle$ differentially in $(\phi^{*}, \mathrm{DCA})$. This macro turns those into **cumulative** profiles, in which bin $j$ holds $\langle R\rangle$ over every candidate with $\mathrm{DCA} >$ the lower edge of bin $j$ -- that is, the value the analysis *would* have measured under that minimum-DCA cut. Reading along the axis is then reading a cut scan.

Three DCA definitions are treated separately, which is what localises the effect: DCA between the V0 daughters, DCA of the proton-like daughter to the PV, and DCA of the pion-like daughter to the PV. Each is additionally split by the sign of $\eta_{\rm jet}$ and of $\eta_{\Lambda}$, since the ring fake signal is known to appear linearly in jet pseudorapidity (maybe a $\tanh\eta_{\text{Jet}}$, but this is very preliminary).

### Design constraints worth preserving

**Accumulate raw statistics, not means.** The running sums are $\big(N,\ \sum Y,\ \sum Y^2\big)$, written straight into the destination profile's own per-bin storage via `SetBinEntries`, `SetBinContent` and `GetSumw2()->SetAt`. ROOT then derives $\langle R\rangle$ and its error itself. Storing a precomputed mean and error instead would look identical on screen but would break `TBrowser` rebinning, which recombines the underlying statistics.

**Output stays a `TProfile`/`TProfile2D`, never a `TH1F`/`TH2F`.** Same reason.

**Axes are transposed on output**: input is $(\phi^{*}, \mathrm{DCA})$, output is $(\mathrm{DCA}, \phi^{*})$. DCA is the variable being scanned and belongs on $x$; $\phi^{*}$ is supporting information and reads better on $y$.

**`Project3DProfile("yx")` is the correct convention here** -- $y$ vertical, $x$ horizontal -- and has been confirmed against the source histograms. Do not "fix" it to `"xy"`.

> **Adjacent points are not independent.** Cumulative bins are nested by construction: the sample at
> a given minimum-DCA cut contains every sample to its right. The per-point error bars are correct as
> per-point uncertainties, but they must **not** be read as though they supported a point-to-point
> comparison, and a $\chi^2$ against a flat line computed from them would be meaningless. Judge
> flatness from the overall trend, not from the scatter.

### Output

`Cumulative_Counts/`, `Cumulative_2D/`, `Cumulative_1D/`, and `Comparison_Canvases/` -- the last
holding 11 canvases: one $\phi^{*}$-integrated comparison of all three DCA types, three jet-$\eta$
splits, three $\Lambda$-$\eta$ splits, and four fixed-hemisphere canvases overlaying all three DCA
types.

---

## `auxiliaryPerConfigPlots.cxx`

**Step 5**, once per config. Output: `results_AuxPerConfig/AuxiliaryPerConfigPlots_<SUFFIX>.root`.

The general home for **per-config derivative plots** -- anything cheap to build from a finished `ConsumerResults_*.root` that does not belong inside the O2Physics consumer and does not need cross-config aggregation. It is deliberately not scoped to one quantity; new sections are expected (see the `ADD MORE POST-PROCESSING SECTIONS HERE` marker in `main()`).

### Coordinate systems used by the polarization maps

The consumer books the polarization maps in **three** rotated frames, under `<folder>/PolMaps/<frame>/`. They are not competing definitions of the same picture: each answers a question the others cannot, and two of them are *not* vector fields at all. This subsection is the reference for why.

Throughout, $\vec p^{\,*}$ is the unit vector along the proton-like daughter's momentum boosted into the $\Lambda$ rest frame, still expressed in **lab Cartesian axes** (a pure boost does not rotate the axes). $\vec P^{*} = (3/\alpha)\,\vec p^{\,*}$ is what the profiles actually hold. $\hat p_\Lambda$ is the unit vector along the $\Lambda$ **lab momentum**, $\hat t$ the jet direction, and

$$\Phi_{AEE} \equiv \phi_{\Lambda} - \phi^{*}_{p}$$

is the azimuthal-efficiency angle (`deltaPhiLambdaProtonStar` in the consumer).

#### The one rule that decides whether to rotate

Every rotated frame here is built from vectors reconstructed per candidate, so the binned plane is invariant under some group of lab rotations. **Any polarization component drawn on that plane must be invariant under the same group, or it averages to zero bin by bin.**

The lab azimuth of an event is random, so the relevant group always contains rotations about $\hat z$. That single observation settles every design question below:

- Components expressed in **lab** axes are *not* $z$-rotation invariant. Overlaying them on a rotated plane gives an arrow field consistent with zero, whose residual measures only the detector's azimuthal non-uniformity.
- Components expressed in the **rotated** axes *are* invariant -- but if the frame was built out of the polarization itself, they are invariant *trivially*, because the rotation froze them.

So: **rotate when the new frame is built from vectors independent of $\vec P^{*}$; do not when the frame is built from $\vec P^{*}$.**

| Frame | $\hat x$ along | Built from | What it is for |
|---|---|---|---|
| `Lab` | detector $x$ | -- | the detector-frame picture, unchanged |
| `Aee` | $\Phi_{AEE} = 0$, i.e. $\vec p^{\,*}_{T}$ | the proton | acceptance maps (**scalars only**) |
| `PrimeV0` | $\vec p^{\,\Lambda}_{T}$ | the $\Lambda$ | the $\Lambda$ production plane (vector field) |
| `PrimeJet` | beam, orthogonalised against the jet | jet + beam | the ring measurement (vector field) |

#### `Aee` -- and why it carries no arrows

The rotation is by $-\phi^{*}_{p}$ about $\hat z$, so that $\hat x_{AEE}$ is the transverse direction of the rest-frame proton. The plane is then polar in disguise: **radius is $p_T^{\Lambda}$ and azimuth is $\Phi_{AEE}$.** The counts map on it is the azimuthal efficiency modulation, resolved in $p_T$ -- the AEE, drawn directly.

But $\hat x_{AEE}$ is aligned with the transverse polarization *by construction*. With $\hat y_{AEE} = \hat z \times \hat x_{AEE} = (-\sin\phi^{*}_{p},\ \cos\phi^{*}_{p},\ 0)$,

$$\vec p^{\,*}\cdot\hat y_{AEE} = p^{*}_{T}\big(-\cos\phi^{*}_{p}\sin\phi^{*}_{p} + \sin\phi^{*}_{p}\cos\phi^{*}_{p}\big) = 0 \quad \text{identically,}$$

so in this frame the polarization is exactly

$$\big(P^{*}_{x'},\ P^{*}_{y'},\ P^{*}_{z'}\big) = \big(|\vec p^{\,*}_{T}|,\ 0,\ P^{*}_{z}\big) \;\propto\; \big(\sin\theta^{*},\ 0,\ \cos\theta^{*}\big),$$

with $\theta^{*}$ measured from the **beam** -- which is the axis $\vec B$ lies along, hence the AEE-relevant one. Every arrow would point along $+\hat x_{AEE}$ with the same frozen direction; only the length would vary.

Leaving the components *unrotated* does not rescue it either: the bin coordinates are $z$-rotation invariant while $P^{*}_x$ and $P^{*}_y$ are not, so both average to zero.

The general statement, which no choice of $\hat x$ escapes: **the transverse plane has one angular degree of freedom, and $\Phi_{AEE}$ already spends it.** Any transverse vector built from $\{\hat p_{\Lambda,T},\ \vec p^{\,*}_T,\ \hat z\}$ has an azimuth that is a fixed function of $\Phi_{AEE}$, in any frame.

What survives is not nothing. $\langle\sin\theta^{*}\rangle$ and $\langle\cos\theta^{*}\rangle$ per $(p_T, \Phi_{AEE})$ bin are the first two moments of the rest-frame polar-angle acceptance in that bin. For a perfect, isotropic detector $\langle\sin\theta^{*}\rangle = \pi/4 \approx 0.785$ and $\langle\cos\theta^{*}\rangle = 0$, both flat; every deviation is acceptance sculpting. Together with the counts map that is a complete low-order characterisation, which is exactly what this family books.

#### `PrimeV0` -- the $\Lambda$ production plane

The rotation is by $-\phi_{\Lambda}$, so $\hat x_{V0}$ is the $\Lambda$'s own transverse direction and $\hat y_{V0} = \hat z \times \hat x_{V0}$. The components are

$$P^{*}_{x'V0} = p^{*}_{T}\cos\Phi_{AEE}, \qquad P^{*}_{y'V0} = -\,p^{*}_{T}\sin\Phi_{AEE}, \qquad P^{*}_{z} \ \text{unchanged}.$$

Two things make this frame worth having.

**First**, $\hat y_{V0} = \widehat{\hat z \times \hat p_{\Lambda}}$ is the ring observable's normal with the **beam substituted for the jet**. So

$$\langle P^{*}_{y'V0}\rangle \;=\; \langle R\rangle \ \text{computed with}\ \hat z\ \text{as the jet proxy},$$

identically, not approximately. That is a control observable of the same kind as the leading-particle and mixed-event proxies: it carries no jet-correlated ring signal, but the full AEE bias. It costs two multiplications.

**Second**, $(p_z,\ p_T)$ *is* the production plane, and $\hat z$ and $\hat x_{V0}$ are its own in-plane unit vectors. Axes and arrows therefore share one basis, with no mixed-frame comparison hiding anywhere:

| | |
|---|---|
| horizontal axis | $p_z^{\Lambda}$, in-plane, along $\hat z$ |
| vertical axis | $p_T^{\Lambda} \ge 0$, in-plane, along $\hat x_{V0}$ |
| arrows | $(\langle P^{*}_z\rangle,\ \langle P^{*}_{x'V0}\rangle)$ |
| COLZ | $\langle P^{*}_{y'V0}\rangle$, out of plane |

#### `PrimeJet` -- the jet frame

$\hat z_{Jet} = \hat t$, and $\hat x_{Jet}$ is the beam direction **orthogonalised against the jet**. With $t_z = \hat t\cdot\hat z$ and $t_T \equiv \sqrt{1-t_z^{2}} = 1/\cosh\eta_{jet}$:

$$\hat z_{Jet} = \hat t, \qquad \hat x_{Jet} = \frac{\hat z - t_z\,\hat t}{t_T}, \qquad \hat y_{Jet} = \hat z_{Jet}\times\hat x_{Jet} = \frac{\hat t \times \hat z}{t_T},$$

giving, for any vector $\vec v$,

$$v_{z'} = \vec v\cdot\hat t, \qquad v_{x'} = \frac{v_z - t_z\,v_{z'}}{t_T}, \qquad v_{y'} = \frac{v_x t_y - v_y t_x}{t_T}.$$

$t_z$ and $1/t_T$ depend only on the jet, so the consumer resolves them **once per collision**, right after `applyProxyDistortion()` (which rewrites `leadingJetUnitVec` in place -- computing the basis before it would use the undistorted jet). The same rotation is applied to $\vec p^{\,*}$ and to $\vec p^{\,\Lambda}$.

**Why the orthogonalisation is not optional.** The naive triad $\{\hat z,\ \hat t\times\hat z,\ \hat t\}$ fails twice: $\hat z\cdot\hat t = \tanh\eta_{jet} \neq 0$ for any jet off midrapidity, and $|\hat t\times\hat z| = 1/\cosh\eta_{jet} \neq 1$. The second is only a normalisation; the first is structural. With $\hat x = \hat z$ the "plane transverse to the jet" is tilted by $\tanh\eta_{jet}$ *along* the jet, so the projection is skewed by an amount correlated with the jet kinematics, which would smear the ring pattern in an $\eta_{jet}$-dependent way. Note that $\hat y_{Jet}$ comes out in the same direction either way -- only $\hat x_{Jet}$ actually changes -- and the whole construction reduces to the naive one at $\eta_{jet}=0$.

**No sign ambiguity.** A frame axis fixed by an arbitrary algorithmic choice would average $\langle P^{*}_{x'}\rangle$ to zero. That does not happen here:

$$\hat x_{Jet}\cdot\hat z = \frac{1-t_z^{2}}{t_T} = t_T > 0 \quad\text{always},$$

and $\hat x_{Jet}$ is a smooth, branch-free function of $\hat t$ alone. It is the unique unit vector perpendicular to the jet, lying in the beam-jet plane, on the $+\hat z$ side. The only degeneracy is $t_T \to 0$, a jet exactly along the beam, where the beam-jet plane does not exist; unreachable for the $|\eta_{jet}|$ the jet finder accepts, and guarded anyway.

### The ring observable *is* the azimuthal component of the jet-frame arrows

This is the reason the $(x'_{Jet},\ y'_{Jet})$ panel is the headline plot.

In the jet frame $\hat t = \hat z'$, so $\widehat{\hat t\times\hat p_\Lambda}$, the ring normal, is the azimuthal unit vector $(-\sin\phi',\ \cos\phi',\ 0)$ with $\phi'$ the azimuth of $\hat p_\Lambda$ in the jet frame. Therefore

$$R \;=\; -P^{*}_{x'}\sin\phi' + P^{*}_{y'}\cos\phi' \;=\; P^{*}_{\phi'}.$$

Consequences worth internalising:

- A genuine ring signal appears as arrows **circulating** about the origin. A radial or divergent pattern is something else entirely.
- The radius of that plane is $|\vec p_{\Lambda}|\sin\Delta\theta_{jet}$ -- literally the radius of the ring.
- It is a **bin-by-bin closure test**: the tangential projection of the drawn arrows must reproduce `p2dRingObservableVsPxPyPrimeJet`, up to binning effects.

### Parity relations, and which ones are supposed to break

Momenta are polar vectors, $\vec P^{*}$ is axial, and $\vec B = B\hat z$ is axial too. Mirror symmetries therefore relate the maps to themselves for free -- and, crucially, the magnetic field breaks exactly one of the two, which is what makes the pair diagnostic.

#### Mirror through the beam-jet plane ($M_1$, the $x'z'$ plane)

The jet frame is invariant under it. For a parity-conserving production mechanism and an $M_1$-symmetric acceptance:

$$\langle P_{x'}\rangle,\ \langle P_{z'}\rangle \ \text{odd in } p_{y'}; \qquad \langle P_{y'}\rangle,\ \langle R\rangle \ \text{even in } p_{y'}.$$

$\vec B$ lies **in** this mirror plane and is axial, so $M_1$ flips it. **This relation is designed to break**: a violation is a magnetic-field-induced bias, which is the thing the AEE study is hunting.

#### Mirror through $z=0$ ($M_z$)

Here the frame does not simply follow the mirror. Working it through: $\hat z'_{img} = M_z\hat z'$ and $\hat y'_{img} = M_z\hat y'$, but $\hat x'_{img} = -M_z\hat x'$, because the definition always points beam-ward while the mirror image of $+\hat z$ is $-\hat z$. The result is

$$\langle P_{x'}\rangle \ \text{even in } p_{x'}; \qquad \langle P_{y'}\rangle,\ \langle P_{z'}\rangle \ \text{odd in } p_{x'}; \qquad \langle R\rangle \ \text{even in } p_{x'}.$$

$\vec B$ is **normal** to this mirror plane and axial, so $M_z$ preserves it. This relation is broken only by genuine detector $z$-asymmetry or collision-system asymmetry, neither of which applies to O--O in the central barrel. **This one should simply hold**, and a violation means something is wrong with the machinery, not with the physics.

So $\langle R\rangle$ is even under both: the ring map should be symmetric about **both** axes of the $(x',y')$ panel. Free closure test on the headline plot, nothing extra to book.

#### In the production plane (`PrimeV0`)

Take $M_{\Lambda}$, the mirror through the $\Lambda$ production plane. It leaves $\vec p_{\Lambda}$ alone and sends $\Phi_{AEE}\to-\Phi_{AEE}$. At fixed $(p_z, p_T)$ we integrate over $\Phi_{AEE}$, so

$$\langle P^{*}_{x'V0}\rangle = \langle P^{*}_{z}\rangle = 0, \qquad \langle P^{*}_{y'V0}\rangle \ \text{unconstrained}.$$

This is the Basel-convention statement: a parity-conserving single-spin polarization must be normal to the production plane. **The arrows in that panel must vanish identically**, so any non-zero arrow field there is B-induced bias -- the cleanest pure null test in the whole set. The COLZ is *not* a null: $\langle P^{*}_{y'V0}\rangle$ is allowed by parity and is the classic transverse $\Lambda$ polarization normal to the production plane.

#### On the AEE planes

Differential in $\Phi_{AEE}$, so parity forces a symmetry rather than a vanishing:

$$\text{counts},\ \langle\sin\theta^{*}\rangle \ \text{even in } \Phi_{AEE}; \qquad \langle\cos\theta^{*}\rangle \ \text{odd in } \Phi_{AEE}.$$

Note this corrects a natural but wrong expectation: the null for the $\langle P^{*}_z\rangle$ map is **odd in $\Phi_{AEE}$**, not flat -- it integrates to zero rather than being zero everywhere. And it gives a quantitative prescription:

$$\text{AEE} \;\equiv\; \text{the odd-in-}\Phi_{AEE}\text{ component of the counts map.}$$

Folding the map about $\Phi_{AEE}=0$ leaves an antisymmetric residue that *is* the effect, isolated from acceptance that is merely non-uniform for geometric reasons. See Section 5.

### Section 1 -- polarization vector fields

$\langle \vec{P}^{*}\rangle$ drawn over momentum planes, in the same visual language as the Helicity Efficiency Toy Model's plotter ["PhD_codes/ToyModels/plotHelicityEfficiency.cxx"](../ToyModels/plotHelicityEfficiency.cxx), so that real data and toy model can be compared by eye.

One canvas per coordinate system, per kinematic-cut folder.

**`cVectorFieldLab`** -- read from `PolMaps/Lab/`:

| Panel | COLZ background | Arrows |
|---|---|---|
| X--Y plane | $\langle P^{*}_z\rangle$ | $(\langle P^{*}_x\rangle,\ \langle P^{*}_y\rangle)$ |
| Z--X plane | $\langle P^{*}_y\rangle$ | $(\langle P^{*}_z\rangle,\ \langle P^{*}_x\rangle)$ |
| Y--Z plane | $\langle P^{*}_x\rangle$ | $(\langle P^{*}_y\rangle,\ \langle P^{*}_z\rangle)$ |

**`cVectorFieldPrimeJet`** -- read from `PolMaps/PrimeJet/`, the same three planes in primed components. The X'--Y' panel is the measurement itself (see above).

**`cVectorFieldPrimeV0`** -- read from `PolMaps/PrimeV0/`: the production-plane panel described above, plus its counts map.

Arrows are **block-averaged** over `arrowBlockSize` x `arrowBlockSize` tiles of the underlying `TProfile2D`, because one arrow per histogram bin is unreadable (and too unstable in error!). Bins below `minEntries` are excluded from the tile average rather than diluting it.

Arrow *length* is normalised to a **percentile** of the tile-magnitude distribution (95th by default) rather than to the maximum: a single noisy tile would otherwise set the scale and shrink every real arrow to nothing. Lengths are then capped at that reference, so an outlier shows up as a full-length arrow instead of one running off the pad, and the reference magnitude is printed on the panel so the scale never has to be guessed.

The COLZ range is symmetric about zero, so the sign of the out-of-plane component reads directly from the colour.

#### The axis-aspect correction, and the bug it fixes

The arrow direction $(b_x, b_y)$ lives in **polarization-component space**, where both components carry the same units, so the drawn angle is meaningful. The arrow tip, however, is placed in **data coordinates**. Displacing by the same number of data units along both axes only renders the correct angle when both axes map the same number of data units per pixel -- and they generally do not.

`DrawVectorFieldPanel()` therefore now converts the vertical displacement through the pad's pixel geometry,

$$\text{aspect} = \frac{y_{\text{range}}/H_{\text{px}}}{x_{\text{range}}/W_{\text{px}}}, \qquad y_2 = y_c + \frac{b_y}{|\vec b|}\,\ell\,\cdot\,\text{aspect},$$

so the arrow has a fixed pixel *length* and the true polarization *angle* on screen. $W_{\text{px}}$ and $H_{\text{px}}$ are the frame extents, i.e. the pad size times the margin-free fractions; `gPad->Update()` must have run first.

This was a genuine, silent bug. With $p_z$ spanning $[-4,4]$ against $p_x$ spanning $[-3,3]$ on a non-square pad, a 45-degree polarization drew at roughly 34 degrees. Every non-square panel ever produced was systematically misread, and since the whole point of these panels is reading directions off the page, that mattered. `aspect = 1` recovers the old behaviour exactly, which is what a square panel with equal ranges would give anyway, so the X--Y panels are unaffected.

Relatedly, all rotated-frame momentum axes now share `axisLambdaPRot`, and `axisLambdaPtRot` is binned at 0.2 GeV/c to match `axisLambdaPz`, so the panels are isotropic in the first place rather than relying on the correction to rescue them.

#### On the out-of-plane sign convention

Every vector-field canvas uses the three **cyclic** planes, read as an ordered pair (horizontal, vertical):

$$(x,y) \to +\hat z, \qquad (y,z) \to +\hat x, \qquad (z,x) \to +\hat y.$$

All three are right-handed by construction, so one rule holds on every panel without exception: **a positive colour points out of the page, toward the reader.** No sign is flipped anywhere in the plotting code -- each COLZ is the raw booked profile -- and the price is only cosmetic: $p_z$ is horizontal on one panel and vertical on another.

The convention lives in the consumer's booking (`_vsPyPz`, `_vsPyAeePz`, `_vsPyPzPrimeJet`), not here. Anyone adding a plane should keep it: a non-cyclic pair such as $(z,y)$ has out-of-plane direction $-\hat x$, and would silently invert the reading of that one panel.

**What to check when adding a vector panel:** the arrow's horizontal component must be the profile of the *horizontal* axis. On a $(y,z)$ panel that is $\langle P^{*}_y\rangle$, not $\langle P^{*}_z\rangle$; carrying over the $(z,x)$ ordering by habit transposes the whole field about the diagonal, and a transposed field still looks plausible. The three component letters of `(hA, hB, hC)` must form a cyclic triple.

### Section 2 -- ring-observable 2D maps

$\langle R\rangle$ over the momentum planes. No arrow overlay: the ring observable is already one scalar per candidate, so the COLZ background *is* the whole result. No arrows here!

$R$ is a scalar and invariant under every rotation above, so **only the binning plane changes** between these canvases -- nothing is recomputed:

| Canvas | Planes |
|---|---|
| `cRingObservable2D` | lab: X--Y, Z--X, Y--Z |
| `cRingObservableAee2D` | Aee: $(x_{AEE}, y_{AEE})$, $(z, x_{AEE})$, $(y_{AEE}, z)$ |
| `cRingObservablePrimeJet2D` | jet frame: X'--Y', Z'--X', Y'--Z' |
| `cRingObservableLeadP2D` | lab, leading-particle proxy |

### Section 3 -- AEE acceptance maps

Three canvases, one per Aee plane, each holding candidate counts, $\langle P^{*}_{T}\rangle$ and $\langle P^{*}_{z}\rangle$. The counts panel of `cAeeMapsXY` is the headline: radius $p_T^{\Lambda}$, azimuth $\Phi_{AEE}$, so any azimuthal structure there is the azimuthal efficiency modulation resolved in $p_T$.

Counts are drawn one-sided (floor pinned at zero) rather than with the diverging symmetric range used for the polarization maps: they are a one-sided quantity, and the symmetric treatment would waste half the palette and put the empty regions in mid-scale.

### Section 4 -- candidate counts

`cCountsLab` and `cCountsPrimeJet`, three planes each. Plain occupancy, useful as the denominator behind every map above and as the first place a detector hole shows up.

### Section 5 -- AEE fold

The parity relations above say what each Aee map must look like under $\Phi_{AEE}\to-\Phi_{AEE}$. On the $(x_{AEE}, y_{AEE})$ plane that reflection is simply $y_{AEE}\to-y_{AEE}$, so the fold is a pairing of bin $i_y$ with bin $N+1-i_y$. `cAeeFold` puts the effect and its three parity nulls side by side:

| Panel | Content | Expectation |
|---|---|---|
| 1 | $\mathcal{A} = \dfrac{N(x,y)-N(x,-y)}{N(x,y)+N(x,-y)}$ | **the AEE itself**; non-zero if present |
| 2 | odd part of $\langle P^{*}_{T}\rangle$ | null: $\langle\sin\theta^{*}\rangle$ is even |
| 3 | even part of $\langle P^{*}_{z}\rangle$ | null: $\langle\cos\theta^{*}\rangle$ is odd |
| 4 | odd part of $\langle R\rangle$ | null: $R$ is even |

The folded histograms are written next to the canvas, so the asymmetry can be projected or integrated downstream rather than only looked at.

**Design constraints worth preserving:**

- **The counts are folded into a normalised asymmetry, not a raw difference.** $\mathcal{A}$ is dimensionless, bounded in $[-1,1]$, reads directly as a fractional efficiency modulation, and for independent Poisson counts has the exact variance $\sigma^{2}_{\mathcal{A}} = 4ab/(a+b)^{3}$.
- **The profile folds need *both* partners populated.** A pair enters only if both bins clear `minEntries`. A half-populated pair would otherwise masquerade as a parity violation -- precisely the failure mode this canvas exists to detect.
- **The output fills the full plane.** It is antisymmetric (or symmetric) by construction, so drawing both halves makes that manifest by eye, at the cost of some redundancy.
- **The axis is checked, not assumed.** Folding needs the vertical axis symmetric about zero with an even bin count. `axisLambdaPRot` qualifies, but it is a `ConfigurableAxis`, so `HasFoldableYAxis()` verifies it at runtime and skips the canvas with a message otherwise.

**What to check:** panels 2--4 are nulls only for a mirror-symmetric setup, and the magnetic field that produces panel 1 is exactly what breaks that mirror. Read the four together; a residue in a null is not automatically a bug.

### Direction cosines, and why the kernel uses them

Sections 6 and 7 bin in direction cosines rather than in pseudorapidity or angle. For pseudorapidity the identity

$$\tanh\eta = \cos\theta$$

is exact, so the $z$ component of any unit vector already *is* $\tanh\eta$. The consumer reads them straight off the vectors it already builds -- `jetZ` $= \hat t_z$, `leadPZ` $= \hat t^{\,\rm LeadP}_z$, `lambdaZ` $= (\hat p_\Lambda)_z$ -- with no transcendental evaluated anywhere. Throughout, $t_z$ is the reference-axis direction cosine and $\lambda_z = \cos\theta_\Lambda$, with $\lambda_T = \sin\theta_\Lambda = \sqrt{1-\lambda_z^{2}}$.

Three further reasons make $\cos\Delta\theta_{jet}$, rather than $\Delta\theta_{jet}$, the right binning variable: the kernel's numerator is linear in it; it is the flat variable of the phase space, which keeps the occupancy even across bins; and it turns the $\Delta\theta\leftrightarrow\pi-\Delta\theta$ pairing into a plain fold about zero.

### Section 6 -- ring projection kernel

> In an $R_z$ file, Sections 6 to 8 are null tests: see [The longitudinal ring](#the-longitudinal-ring-r_z-useringz).

Sections 6 to 8 rest on the geometry derived in full in the ring-geometry note (`18 - RingSymmetriesAndTanh.pdf`). This section summarises what the code relies on, and nothing more; every statement below is proved there.

#### The kernel

Notation: $t_z$ is the reference-axis direction cosine (`jetZ` or `leadPZ`, accordingly), $\lambda_z=\cos\theta_\Lambda$ (`lambdaZ`), $\lambda_T=\sin\theta_\Lambda$, and $c=\cos\Delta\theta$. Since $\tanh\eta=\cos\theta$ exactly, $t_z=\tanh\eta_t$ and $\lambda_z=\tanh\eta_\Lambda$.

After an exact average over the decay, the mean ring at fixed kinematics is the projection $\vec P_{\rm meas}\cdot\hat n$ of one effective polarization vector, for *any* acceptance. Decomposing $\vec P_{\rm meas}$ in the $\Lambda$'s local triad $(\hat p_\Lambda,\ \hat\theta,\ \hat\varphi)$, with $\hat\varphi=\widehat{\hat z\times\hat p_\Lambda}$ (the `PrimeV0` frame's $\hat y_{V0}$) and $\hat\theta=\hat\varphi\times\hat p_\Lambda$, gives, exactly and configuration by configuration,

$$\langle R\rangle \;=\; P_\varphi\cos\chi \;+\; P_\theta\sin\chi, \qquad \cos\chi \;=\; \frac{t_z-\lambda_z c}{\lambda_T\sqrt{1-c^{2}}},$$

where $\chi$ is the bearing of the reference axis as seen from the $\Lambda$, with the beam as north. $\cos\chi$ is the spherical law of cosines, and it is the origin of every $\tanh\eta_{\rm jet}$ in the ring. The sign of $\sin\chi$ records on which side of the $\Lambda$'s meridian the reference axis lies (east or west).

Three structural facts follow, all properties of the projection rather than of any dataset:

- **The helicity component never enters.** $\hat p_\Lambda\perp\hat n$ identically.
- **$p_T^\Lambda$ never enters the geometry.** $\hat n$ is built from unit vectors only.
- **$P_\theta$ enters only through the branch sign.** In a cell of fixed $(t_z,\lambda_z,c)$, $|\sin\chi|$ is fixed and only its sign varies, so $P_\theta$ contributes $|\sin\chi|\langle\sigma P_\theta\rangle$. It drops out only if the polarization does not know where the reference axis is *and* the two branches are equally populated. A genuine ring puts part of its signal into $P_\theta$, so this term is not negligible in general.

**At fixed $\lambda_z$**, for an efficiency that does not know where the reference axis is and balanced branches, the kernel is exactly affine in $t_z$:

$$\langle R\rangle \;=\; \bar P_\varphi(\lambda_z)\,\frac{t_z-\lambda_z c}{\lambda_T\sqrt{1-c^{2}}}.$$

This is what Section 7 uses.

**Integrated over $\lambda_z$** -- which is what the 2D profile does -- the $P_\varphi$ term still reduces exactly to two moments,

$$\langle R\rangle(c,t_z) \;=\; \frac{M_0\,t_z-M_1\,c}{\sqrt{1-c^{2}}} \;+\; (\text{branch term}) \;+\; S, \qquad M_0=\big\langle P_\varphi\cosh\eta_\Lambda\big\rangle_{c,t_z},\quad M_1=\big\langle P_\varphi\sinh\eta_\Lambda\big\rangle_{c,t_z},$$

but the moments are **cell-local**: averages over the $\Lambda$s in that $(c,t_z)$ bin. They vary across the surface even for a constant $P_\varphi$ and a perfectly symmetric detector, because at fixed $(c,t_z)$ a $\Lambda$ enters with phase-space weight $1/\sqrt{G}$, where $G$ is the Gram determinant of beam, reference axis and $\Lambda$ -- and $G$ depends on $t_z$. Two consequences:

- the 2D surface is not exactly affine in $t_z$;
- $M_1$ has a purely geometric part, odd in $t_z$ and in $c$, which vanishes at $t_z=0$ and at $c=0$. A non-zero $M_1$ is therefore not by itself a detector asymmetry. Its value at $t_z=0$, where the geometric part vanishes, can come from an $\eta_\Lambda$-asymmetric yield, an $\eta_\Lambda$-asymmetric acceptance, or a genuine polarization normal to the production plane ($P_N$, odd in $\eta_\Lambda$) -- and only the first two are detector effects.

#### Input

Booked by the consumer under `<folder>/LeadJet/RingKernel/` and `<folder>/LeadP/RingKernel/`, jet-gated via `RING_OBSERVABLE_FILL_LIST` and `RING_OBSERVABLE_LEADP_FILL_LIST`:

| Histogram | Axes |
|---|---|
| `p2dRingObservableCosDeltaThetaVsJetZ` | `axisCosTheta` x `axisJetZ` |
| `h2dCountsCosDeltaThetaVsJetZ` | same |
| `p3dRingObservableCosDeltaThetaVsJetZVsLambdaZ` | `axisCosThetaCoarse` x `axisJetZ` x `axisLambdaZ` |
| `p2dRingObservableLeadPCosDeltaThetaVsLeadPZ` | `axisCosTheta` x `axisLeadPZ` |
| `h2dCountsLeadPCosDeltaThetaVsLeadPZ` | same |
| `p3dRingObservableLeadPCosDeltaThetaVsLeadPZVsLambdaZ` | `axisCosThetaCoarse` x `axisLeadPZ` x `axisLambdaZ` |

The two proxy axes differ on purpose. FastJet at $R = 0.4$ inside $|\eta| < 0.9$ confines jets to $|\eta| < 0.5$, i.e. $|t_z| \le \tanh 0.5$; leading particles reach $|\eta| < 0.9$. `axisJetZ` and `axisLeadPZ` follow each proxy's fiducial range rather than wasting bins on structurally empty space. `axisProxyZ` keeps the older full-range binning for `pRingVsJetZcomponent`. Every axis used by Sections 7 and 8 must be symmetric about zero with an even bin count; both sections verify this at runtime.

The consumer only books and fills these. **All fitting happens here**, per the ALICE convention that no post-processing lives in an O2Physics task.

#### What Section 6 does

`FitRingKernelSlices()` groups `cosRebin` bins of $\cos\Delta\theta$, collapses each group into a $t_z$ profile with `ProfileY`, and fits a straight line. Because the moments depend on $t_z$, the slope is not $M_0/\sin\Delta\theta$ but

$$\text{slope}\times\sin\Delta\theta \;=\; M_0 \;+\; t_z\,\frac{\partial M_0}{\partial t_z} \;-\; c\,\frac{\partial M_1}{\partial t_z} \;+\; (\text{branch term}),$$

and since the geometric $M_1$ is odd in $c$, the correction $-c\,\partial M_1/\partial t_z$ is **even** in $c$. The fit range in $t_z$ is symmetric about zero, so the intercept is the value at $t_z=0$, where the geometric $M_1$ vanishes.

`cRingKernel` shows three panels:

| Panel | Content | Prediction |
|---|---|---|
| 1 | $\langle R\rangle$ over $(\cos\Delta\theta,\ t_z)$, leading jet | a tilted surface, steepening toward $\vert\cos\Delta\theta\vert\to1$ |
| 2 | slope $\times\sin\Delta\theta$, jet and leading particle overlaid | **a bowl**, even in $\cos\Delta\theta$, with its minimum at $\cos\Delta\theta=0$ close to $M_0$ |
| 3 | intercept $\times\sin\Delta\theta$ | $-M_1(0,c)\,c+S(c)\sin\Delta\theta$: odd in $\cos\Delta\theta$ whenever $S=0$ |

The dashed lines on panels 2 and 3 mark each proxy's **fully allowed window**, $|\cos\Delta\theta|\le-\cos(\vartheta_t+\vartheta_\Lambda)$, with $\vartheta$ the smallest polar angle in each acceptance. Inside it no $\Lambda$ in the acceptance is kinematically excluded; outside it, kinematic truncation steepens the bowl further, in a way that depends on how the bins sample the edge of the domain, and no simple shape should be expected there. The acceptance edges are read off the 3D profile's own occupancy -- conservatively, at the outer edge of the outermost populated bin -- so the window follows whatever the consumer actually accepted instead of repeating its cuts here.

Per-slice fit canvases go to `RingKernel_Canvases/Fits_<folder>/`, and the four `TGraphErrors` are written alongside for downstream use.

**Design constraints worth preserving:**

- **A straight-line fit, with no Minos.** At fixed $\cos\Delta\theta$ the model is a line in $t_z$ (exactly so at fixed $\lambda_z$), so there is no model choice to make. The $\tanh$ fits in `auxiliarySummaryPlots.cxx` need Minos because their two parameters collapse onto a curved valley at small argument; a line has no such valley, and with a $t_z$ axis symmetric about zero the slope and intercept are nearly uncorrelated. Chi-squared is the right method, since the bins are means carrying SEM errors, not Poisson counts.
- **Do not fit a $\tanh$ here.** The $\tanh$ shape is what the $t_z$-linear kernel turns into *after* integrating over $\Delta\theta$ and the $\eta_\Lambda$ distribution. Fitting it to this surface would fit a consequence of the model instead of the model.
- **Panel 2 is not a flatness test.** The bowl is geometry, present for a constant $P_\varphi$ and a perfect detector. What *is* informative is its symmetry about $\cos\Delta\theta=0$.
- **$M_1$ is not separated from $S$ here.** Doing so by parity in $\cos\Delta\theta$ would need $S$ to be even about $\cos\Delta\theta=0$, which nothing guarantees -- the recoil jet sits near $\pi-\Delta\theta$. Section 8 performs the separation properly, at fixed $\lambda_z$, where it is exact.

Defaults (`cosRebin = 5`, `maxAbsCos = 0.95`, `minEntries = 200`) are parameters of `MakeRingKernelCanvas()`. The $|\cos\Delta\theta|$ ceiling exists because the kernel diverges as $1/\sin\Delta\theta$ there, and the acceptance has already emptied those bins.

**What to check:**

- Panel 2 symmetric about $\cos\Delta\theta=0$. An odd component points at a near-side correlation between reference axis and $\Lambda$, or at an efficiency that knows where the reference axis is.
- The minimum of panel 2 against Section 7's $M_0$. They should be close but need not be equal: the minimum is the cell-local $M_0$ at $(c,t_z)\approx(0,0)$, weighted by the phase-space Jacobian there, while Section 7's is a sample-wide average. The exact comparison predicts the whole bowl from Section 7's $K(\lambda_z)$ and the 3D occupancy (not yet implemented).
- Jet against leading particle. For azimuthally uncorrelated proxies their minima should agree, but their bowls need not: the fitted slope depends on the $t_z$ range, which differs, and on each proxy's correlation with the $\Lambda$. Proxy independence is exact only at fixed $\lambda_z$ (Section 7).
- Panel 3 odd in $\cos\Delta\theta$. An even part is signal, or a branch term.
- The per-slice $\chi^{2}/\text{ndf}$ in `Fits_<folder>/`, and the counts histogram, to see where the plane is kinematically forbidden before reading anything into its edges.

### Section 7 -- kernel moments

Section 6's bowl brackets $M_0$ without pinning it down. Section 7 goes behind it: it recovers the axis-independent azimuthal component **differentially in $\cos\theta_\Lambda$** from the 3D profile, where the kernel is exactly affine in $t_z$, and builds both moments from it.

The code calls this component `Bphi`. It is the whole axis-independent azimuthal component of $\vec P_{\rm meas}$: the detector-induced part *and* any genuine $P_N$. The two are not separated here.

#### The extraction

At fixed $(\lambda_z,\ c)$ the kernel's slope in $t_z$ is $m = K/\sqrt{1-c^{2}}$, with $K = \bar P_\varphi/\lambda_T$. Together with $\cosh\eta_\Lambda = 1/\lambda_T$ and $\sinh\eta_\Lambda = \lambda_z/\lambda_T$, this gives

$$m\sqrt{1-c^{2}} = K \quad \text{(constant in } c\text{)}, \qquad M_0 = \langle K\rangle_{w}, \qquad M_1 = \langle \lambda_z K\rangle_{w},$$

with $w(\lambda_z)$ the $\Lambda$ occupancy of each slab over the whole surface. **The 3D profile alone carries everything**: the moments through its $\lambda_z$ axis, the weights through its own per-bin entry counts, and the model test through the redundancy across $\cos\Delta\theta$ columns. Nothing else needs booking.

These moments are **sample-wide**, unlike Section 6's cell-local ones: the weights are whole-slab occupancies, so the geometric, $t_z$-odd part of the cell-local $M_1$ cancels. What remains in $M_1$ is the $\eta_\Lambda$-odd content of the sample: an A/C yield asymmetry, an A/C acceptance asymmetry, or a genuine $P_N$.

Each $\cos\Delta\theta$ column gives an independent estimate of the same $K$. They are combined with inverse-variance weights, and their scatter about the combination is kept as a $\chi^{2}$. That scatter *is* the model test: at fixed $\lambda_z$ the kernel says $m\sqrt{1-c^{2}}$ does not depend on $c$ at all. Like Section 6, the extraction is immune to a genuine ring, because $K$ comes from a $t_z$ slope and a ring of constant magnitude is flat in $t_z$.

**Symmetrising.** Replacing $w(\lambda_z)$ by $\min\big(w(\lambda_z),\ w(-\lambda_z)\big)$ gives the largest mirror-symmetric weight set that never up-weights a slab beyond what was recorded. It removes the yield asymmetry, and nothing else. The residual $M_1^{\rm sym}$ is non-zero if $K$ has a part odd in $\eta_\Lambda$, which can be an A/C *acceptance* asymmetry or a genuine $P_N$. Field reversal separates those two: a detector-induced $P_\varphi$ is odd in the field, a genuine $P_N$ is not.

#### What Section 7 does

`cKernelMoments` shows three panels:

| Panel | Content | Prediction |
|---|---|---|
| 1 | $\bar P_\varphi(\lambda_z)$, leading jet, with its own mirror image overlaid | coincide if detector-induced with an A/C-symmetric acceptance; an odd part is an A/C acceptance difference or a genuine $P_N$ |
| 2 | $\bar P_\varphi(\lambda_z)$, jet and leading particle | coincide -- at fixed $\lambda_z$ the component does not know the proxy |
| 3 | $K = \bar P_\varphi\cosh\eta_\Lambda$, with $M_0$, $M_1$, $M_1^{\rm sym}$, the yield asymmetry and $M_0^{\rm LeadP}$ | read together with panel 1 |

The same numbers are printed to stdout per folder, since they are what the cross-config stage consumes.

**Design constraints worth preserving:**

- **The per-column fits are closed-form, not Minuit.** `WeightedLineFit()` solves the two-parameter weighted least squares exactly. There are hundreds of these per folder and none is meant to be inspected, so routing each through `TF1` would buy only overhead and transient ROOT objects. Section 6's fits go through ROOT precisely because their canvases *are* meant to be looked at.
- **The $\cos\theta_\Lambda$ axis must mirror-pair.** Symmetrising pairs slab $i_z$ with slab $N+1-i_z$, which `axisLambdaZ` supports. It is a `ConfigurableAxis`, so `ComputeKernelMoments()` verifies the pairing and reports $M_1^{\rm sym}$ as unavailable rather than computing it on a mis-paired axis.
- **The validation used a known truth.** The extraction was checked on a synthetic profile carrying a prescribed azimuthal component, a deliberate yield asymmetry, and a non-trivial ring. Re-run such a check after any change to the extraction: the ring must not leak into $K$, and symmetrising must remove the yield asymmetry without moving $M_0$ appreciably.

Defaults (`maxAbsCos = 0.95`, `minEntriesCell = 30`) are parameters of `MakeKernelMomentsCanvas()`. The per-cell threshold is lower than Section 6's per-slice one because a 3D cell holds far fewer entries.

**What to check:**

- Panel 1: filled and open markers overlie. A systematic gap is an A/C acceptance difference or a genuine $P_N$; compare field polarities to tell which.
- Panel 2: jet and leading particle agree. Here, unlike Section 6, this is exact.
- $M_0$ and $M_0^{\rm LeadP}$ agree, for the same reason.
- $M_1^{\rm sym}$ small compared with $M_1$ when the $M_1$ comes from yield alone.
- The per-slab scatter $\chi^{2}$ does not grow systematically with $|\lambda_z|$.

**An independent cross-check, not yet implemented.** $\bar P_\varphi$ is also measured directly by `p2dPyStarPrimeV0_vsPzPt`, from which $\cosh\eta_\Lambda = |p|/p_T$ and $\sinh\eta_\Lambda = p_z/p_T$ follow bin by bin. That route needs no fitting at all. Two caveats: `PolMaps/PrimeV0` is filled from the ungated `POLARIZATION_PROFILE_FILL_LIST`, so it includes events without a leading jet; and it is an *axis-blind* map, into which a genuine ring leaks unless the reference axes accompanying each $\Lambda$ are antipodally balanced.

### Section 8 -- symmetry decomposition of the kernel

Two reflections act on the $(c,\ t_z,\ \lambda_z)$ cells of the 3D kernel profile:

| Operation | Action on the cell | Physical meaning |
|---|---|---|
| $M_z$ | $(c,\ t_z,\ \lambda_z)\to(c,\ -t_z,\ -\lambda_z)$ | reflection of the event through $z=0$ |
| $\mathcal A$ | $(c,\ t_z,\ \lambda_z)\to(-c,\ -t_z,\ \lambda_z)$ | reference axis replaced by its antipode, $\hat t\to-\hat t$ |
| $M_z\mathcal A$ | $(c,\ t_z,\ \lambda_z)\to(-c,\ t_z,\ -\lambda_z)$ | both |

With the identity they form a four-element group, so any cell function $F$ splits exactly into four sectors, labelled by the sign under $M_z$ and under $\mathcal A$:

$$F_{ab} \;=\; \tfrac14\big(F + a\,F\circ M_z + b\,F\circ\mathcal A + ab\,F\circ M_z\mathcal A\big), \qquad a,b=\pm1.$$

What lands in each sector follows from two facts. A detector-induced polarization is a *polar* vector and a genuine one an *axial* vector, so under any reflection that is a symmetry of the setup the first gives an odd ring and the second an even one. And $\hat n$ reverses under $\mathcal A$, so anything that does not know where the reference axis is gives a ring that is odd under $\mathcal A$ -- with no condition on the detector at all.

| Sector | $M_z$ | $\mathcal A$ | Contains |
|---|---|---|---|
| `(+,+)` | even | even | a genuine ring, averaged between $\Delta\theta$ and $\pi-\Delta\theta$ |
| `(+,-)` | even | odd | a genuine $P_N$; the near/away-odd part of a genuine ring; detector effects leaking through an A/C-asymmetric acceptance |
| `(-,-)` | odd | odd | detector-induced contributions that do not know where the reference axis is |
| `(-,+)` | odd | even | **nothing -- a null test** |

The assignments rest on two assumptions, one per reflection. $M_z$ requires a barrel whose two sides have equal acceptance, for the detector part. $\mathcal A$ requires an efficiency that does not know where the reference axis is. Each tolerates exactly what the other forbids -- $\mathcal A$ is indifferent to A/C asymmetry, $M_z$ to a jet-aware efficiency -- which is what makes the null sector informative. It was validated on a synthetic profile carrying a detector-induced component, a genuine $P_N$ and a genuine ring with different near- and away-side magnitudes: each landed in its sector, and the null stayed at the finite-bin level.

#### What Section 8 does

For each proxy, `cKernelSymmetry{Jet,LeadP}` shows the four sectors over $(\cos\Delta\theta,\ t_z)$, each averaged over $\lambda_z$ with the cells' entry counts as weights. The four sector maps are written alongside, and the $\chi^{2}$ of each sector against zero -- counted once per orbit of the group, since the four cells of an orbit carry the same information -- is printed to stdout.

**Design constraints worth preserving:**

- **The decomposition is done at fixed $\lambda_z$, and averaged only afterwards.** At fixed $\lambda_z$ the statements are exact per cell. Integrating over $\lambda_z$ first would not be: the phase-space weights of a cell and of its antipodal image differ as soon as the reference axis is correlated in azimuth with the $\Lambda$. Averaging the *sectors* afterwards is harmless, because a weighted average of sector values stays in its sector.
- **A cell is used only if all four members of its orbit clear `minEntriesCell`.** A partially populated orbit would put a missing mirror into a sector as a spurious signal.
- **All three axes are verified to mirror-pair.** The maps pair bin $b$ with $N+1-b$ on each axis; `IsMirrorPairable()` checks every edge, so a variable-width `ConfigurableAxis` is caught as well.

**What to check:**

- `(-,+)` consistent with zero. A non-zero null is a detector effect that knows where the reference axis is, or a finite-bin residual, which appears wherever the population varies across a cell differently from its images -- a near-side correlation does this, and it shrinks with the bin width.
- `(-,-)` against Sections 6 and 7: the same detector-induced component, isolated here without fitting.
- `(+,+)` as the signal estimator: free of every axis-independent detector effect, under the $\mathcal A$ assumption alone.
- `(+,-)` separated further by field reversal: a genuine $P_N$ is even in the field, a detector leak through an A/C asymmetry odd.
- The same four sectors for jet and leading particle. The `(-,-)` sector should agree between them; the genuine sectors need not.

### Section 9 -- KappaEff

Written for both ring definitions, one canvas per family, one panel per proxy: the **raw** $\kappa$ per mass bin, $3\langle u^2\rangle/\langle w\rangle$, signal and background mixed. It is a QA of $\kappa$'s mass dependence, not a result -- $\kappa_S$ belongs to the signal extraction, like $\langle R\rangle_S$.

Both means come from the same candidates, so the first-order error keeps their covariance,

$$\mathrm{Var}(\kappa)=\frac{9}{n}\,\frac{\sigma_U^2-2\,(U/W)\,C+(U/W)^2\sigma_W^2}{W^2},\qquad C=\langle u^2w\rangle-UW,$$

with $U=\langle u^2\rangle$, $W=\langle w\rangle$ and the spreads $\sigma$ read from the profiles themselves. For the full ring $w\equiv1$, so $\sigma_W=C=0$ and this is the plain standard error of $3\langle u^2\rangle$. The histograms are written as `hKappaVsMass<Proxy>_<Family>`.

### Section 10 -- $R_z$ diagnostics

Only for `_useRingZ` files, into `RzDiagnostics_Canvases/`. The physics is in [The longitudinal ring](#the-longitudinal-ring-r_z-useringz); what the code does:

- **East/west fold** (per family and proxy). The existing $\langle R_z\rangle(\Delta\varphi)$ is folded into its even (signal) and odd (leak) parts vs $\lvert\Delta\varphi\rvert$. The two halves are merged with `TProfile::Rebin`, which is exact, and summarised as: all, $\Delta\varphi>0$, $\Delta\varphi<0$, symmetrized $\tfrac12(R_++R_-)$, half-difference $D=\tfrac12(R_+-R_-)$, and the leak $A_{EW}D$. The identity $\langle R_z\rangle-\text{symmetrized}=A_{EW}D$ is exact, and $A_{EW}=(N_+-N_-)/(N_++N_-)$ comes from the profile's own entries, so it describes exactly the candidates that make $\langle R_z\rangle$.
- **$\chi$ sectors** (per proxy, task level). $f_{s_a s_b}(\chi)=\tfrac14\bigl[f(\chi)+s_af(-\chi)+s_bf(\pi-\chi)+s_as_bf(\chi-\pi)\bigr]$ on $\chi\in[0,\pi/2]$, fitted with $\sin^2\chi$, $\sin\chi$, $\sin2\chi$ and a constant. The shapes are those of the weights with $\theta_\Lambda$ averaged, so the amplitudes are summaries, not measurements. The axis must be uniform on $[-\pi,\pi]$ with a multiple of 4 bins, so every orbit partner is a whole bin.
- **$\varphi_{\rm AEE}$ flatness** (LeadJet and LeadP). $\langle R_z\rangle$ and $\langle R_\perp\rangle$ overlaid, with a constant fitted to $\langle R_z\rangle$; its $\chi^2$ is the test. Both profiles fill on the same candidates, but the $R_z$ one lives under `HelicityEfficiencyQA/`, so the canvas is skipped quietly when that QA is off.

The symmetrized estimator removes the leak **without external input** under one assumption: at fixed $\lvert\Delta\varphi\rvert$, the fake does not know which side of the $\Lambda$ the proxy sits on. That holds exactly for mixed events. In data only a mirror-breaking effect -- the field, again -- can violate it, and the $(+,+)$ and $(+,-)$ sectors are where it would show.

**Mode check.** `_useRingZ` in the file name and the presence of `RzDiagnostics/` must agree; a mismatch exits `1`. A mislabelled file would otherwise be drawn, and compared downstream, as the wrong observable.

### On the duplicated `DrawVectorFieldPanel()`

It is copied almost verbatim from `plotHelicityEfficiency.cxx`, deliberately. The toy-model plotter and this post-processor are independent workflows, separately compiled and separately run, and the function is small, self-contained and stable. Sharing it through a common header would couple two otherwise unrelated build targets for very little. If it ever needs to diverge meaningfully, or grows non-trivially, that is the trigger to revisit though.

The axis-aspect correction is **not** such a divergence -- it is a bug fix, and both copies carry it. Keep it that way: a fix applied to only one copy would make the toy model and the data stop being comparable by eye on every non-square plane, silently.

### Panel registry

Canvases are described as data, not code: a `PanelSpec` table (kind, histogram names, axis titles, banner) consumed by a single `MakePanelCanvas()` that handles fetching, missing-object reporting, layout and labelling. A new coordinate system therefore costs one table, not one function. Axis-title fragments are named once and shared, so the Lab, Aee, PrimeV0 and PrimeJet versions of the same plane cannot drift apart; the ring panels are generated from a prefix plus plane suffixes rather than written out four times.

Fetching is **all-or-nothing per canvas**. The folder is already known to exist (`ScanPresentFolders()` resolved it), so a missing histogram is a genuine anomaly -- most likely a consumer/macro name drift -- and half a canvas would hide it rather than show it. The warning names the first offender so the drift is one `grep` away.

---

## `auxiliarySummaryPlots.cxx`

**Step 6**, once per wagon, after every config has been processed. Output:
`results_consumer/auxiliarySummaryPlots.root`.

The cross-configuration aggregator. Where each per-config macro looks at one
`ConsumerResults_*.root`, this one reads *all* of a wagon's outputs plus external references and
draws them on shared canvases, so that systematic variations can be compared rather than merely
inspected one at a time.

It combines:

- the **hyperon selections** ($\Lambda$, $\bar\Lambda$, both), each its own set of consumer outputs;
- the **systematic variations** within a selection (data-like jet, randomised jet, the   minimum-$p_{\rm T}$ gate study on artificial proxies, ...);
- a **Monte Carlo reference** and a **pp baseline**, when their directories are configured;
- the **Helicity Toy Model** curve, when its path is configured.

Observables are declared in one table (`kObservableGroups` in the main function). Each entry carries its own in-file directory, so profiles living in `EtaDependence/`, `ProxyPtDependence/`, the cut folder root, or task-level folders such as `EtaStudy/` all flow through the same machinery. Adding an observable means adding a row, not a code path.

Beyond the overlays it produces subtracted curves against the data reference, folded versions where
symmetry makes that meaningful, $\tanh$ fits of the $\eta$ dependence, mass signal-versus-background
splits, the azimuthal- and helicity-efficiency cross-checks, the signal-extraction results read back
from `signalExtractionRing`, and a per-family integrated summary with a cross-family canvas above it.

`sigExtractDir` is the seventh argument and defaults to `<consumerDir>/../results_SigExtract`, which
matches the pipeline layout, so `run_all_wagons.sh` needs no extra argument. Pass `"none"` to skip
every signal-extraction plot.

### $R_z$ families

The summary reads the ring definition of every input from its **name**, so one call covers both. The consumer configs carry `_useRingZ` right after the family, so an $R_z$ output is simply another family whose suffix ends in the tag: `BothHyperons_useRingZ` plus `_MixedEventProxies` is exactly the file an ordinary variation lookup builds. Each family whose $R_z$ data file exists therefore gets an **$R_z$ sibling**, flagged `isRingZ` and written to `<Family>_RingZ/`, and everything -- variations, the MC and pp references, the signal-extraction outputs -- resolves for it with no separate code path. Nothing changes for a wagon without $R_z$ outputs.

The flag changes four things, and nothing else:

- **Absent variations.** Only a subset is run in $R_z$. Their absence is expected, so it is recorded once per $R_z$ family, quietly, with one log line, instead of warned about per observable.
- **Toy Model.** Off for $R_z$ families: it predicts the full ring. The references list keeps its toy entry with an empty path, so every per-reference row stays aligned.
- **Cross-family blocks** (`CrossFamily_Integrated/`) stay full-ring. $R_z$ is compared against the full ring in `RingVsRingZ/` instead.
- **Titles.** They come from many literals (`"<R>"`, `"Integrated <R>_{S}"`, a bare `"R"` y title and `"#Delta"` plus that title). Rather than touching each, an $R_z$ family's folder is **relabelled in one pass** once it is complete: every canvas, histogram and graph under it is read back, its titles, axis titles, bin labels, legend entries and text lines are swapped to $R_z$, and it is rewritten in place. Only whole-token forms change, so a jet-radius "R = 0.4" is left alone, and no other folder is touched.

Two folders need both definitions at once. They are drawn at the top level from the full-ring families, pairing each file with its `_useRingZ` sibling, and each uses only the files it finds.

#### What both folders are built from

Everything is computed from candidate moments -- $n$, $\sum y$, $\sum y^2$ -- rebuilt from each cell's mean and spread of mass-binned profiles, every bin and flow included. All of a proxy's profiles are filled on exactly its candidates, so they share one candidate set: the ring comes from each proxy's own 1D mass profile (`<Proxy>/pRingObservable*Mass`), and the $\kappa$ moments from `<Proxy>/KappaEff/`. **Nothing here is signal-extracted:** these are all candidates, signal and background mixed. The signal-extracted versions come from `signalExtractionRing`.

A ratio of two means over the same candidates keeps its covariance,

$$\mathrm{Var}(r)=\frac{\sigma_a^2-2rC+r^2\sigma_b^2}{n\,\langle b\rangle^2},\qquad r=\frac{\langle a\rangle}{\langle b\rangle},\qquad C=\langle ab\rangle-\langle a\rangle\langle b\rangle ,$$

which is used for both $\langle R_z\rangle/\langle n_z^2\rangle$ and $\kappa=3\langle u^2\rangle/\langle w\rangle$.

#### `RingVsRingZ/` -- the headline check

Per family and proxy, one column per config present in **both** definitions (Data, every variation, MC and pp), three numbers:

| Series | Meaning |
|---|---|
| $\langle R\rangle$ | the full ring, fake included |
| $\langle R_z\rangle$ | the longitudinal ring |
| $\langle R_z\rangle/\langle n_z^2\rangle$ | the ring density per candidate, comparable with $\langle R\rangle$ once the latter's fake is removed (under ring alignment) |

If $R_z$ does what it should, $\langle R_z\rangle$ is about zero for MixedEv and pp and survives only in data, while $\langle R\rangle$ carries its large fake everywhere. A second canvas repeats $\langle R_z\rangle/\langle n_z^2\rangle$ alone, with its Data $-$ Var differences.

#### `Corrections/` -- $\kappa$ and the corrected ring, raw and signal-extracted

The fake is subtracted and the response divided out, reducing both definitions to the same quantity, the ring density $P_R$:

$$P^{\rm full}=\frac{R_d-R_f}{\kappa_{\rm eff}-\tfrac{\alpha^2}{3}R_dR_f}\quad\text{(exact inversion)},\qquad P^{z}=\frac{R_{z,d}-R_{z,f}}{3\langle u^2\rangle}=\frac{R_{z,\rm true}}{\langle n_z^2\rangle}\quad\text{(linear)},$$

with $d$ the data sample, $f$ its mixed-event partner, and $\kappa$ and $\langle u^2\rangle$ measured on the data sample. The pairs are Data/MixedEv, the AN-cuts pair, and each MC or pp reference against its own `_MixedEventProxies` run when that exists. $R_z$ uses the linear form, and **not** because its own fake is small. The exact inversion is $R_{z,\rm true}=\bigl[R_{z,d}(1+\langle\delta\rangle_0)-R_{z,f}\bigr]/\kappa_z$, where $\langle\delta\rangle_0\approx\tfrac{\alpha^2}{3}R_{\rm true}R_{\rm fake}$ is the normalization of the whole sample, made of **full-ring** quantities whatever $R_z$'s own fake is. It multiplies $R_{z,d}$, so the linear form drops $\langle\delta\rangle_0R_{z,d}\sim10^{-4}R_{z,d}$: negligible even when $R_z$'s fake is as large as its signal. The full ring cannot drop the same term, because there it multiplies $R_d\approx0.8$. Plugging $R_z$ quantities into the full-ring formula would be wrong in form, since its denominator is written in terms of the full ring's fake.

Since $\kappa_{\rm eff}=3\langle u^2\rangle$ for the full ring and $\kappa_z\langle n_z^2\rangle=3\langle u^2\rangle$ for $R_z$, both are **one expression**, $P=(R_d-R_f)/(3\langle u^2\rangle-k\,R_dR_f)$, with $k=\alpha^2/3$ for the full ring and $k=0$ for $R_z$. One routine computes it; only its inputs change:

| Folder | $\langle R\rangle$ | $\langle u^2\rangle$, $\langle w\rangle$ |
|---|---|---|
| `Corrections/Raw/` | all candidates, from the consumer's mass-binned profiles | all candidates, from `KappaEff/` |
| `Corrections/SignalExtracted/` | $\langle R\rangle_S$, from `IntegratedSummary/` | $\langle u^2\rangle_S$, $\langle w\rangle_S$, from `KappaEff/` of `signalExtractionRing` |
| `Corrections/CheapSigExtract/` | peak $-$ sideband, from `<Proxy>/KappaEff/*VsMassRegion` | the same, from the same profiles |

A reference's signal extraction is looked for where the pipeline puts it, `<base>/../results_SigExtract`; when it is absent, that reference simply has no signal-extracted column.

**The cheap estimator.** For any per-candidate quantity $y$ ($R$, $u^2$ or $w$),
$$\langle y\rangle_S=\frac{\sum_Py-\sum_By}{N_P-N_B},$$
with $P$ the peak window and $B$ the sideband. The consumer enforces equal total widths, so for a background whose density and $\langle y\rangle_B$ are both linear in mass, the sideband holds exactly the background under the peak. What this misses is the curvature of their product, a linear density times a linear $\langle y\rangle_B$; separate left and right sidebands would recover it, and are the natural extension. The two regions are disjoint and each count is treated as Poisson, so with $A=\sum_P-\sum_B$, $D=N_P-N_B$ and $r=A/D$,
$$\mathrm{Var}(\langle y\rangle_S)=\frac{\mathrm{Var}(A)-2r\,\mathrm{Cov}(A,D)+r^2\,\mathrm{Var}(D)}{D^2},\quad \mathrm{Var}(A)=\textstyle\sum_Py^2+\sum_By^2,\ \mathrm{Var}(D)=N_P+N_B,\ \mathrm{Cov}(A,D)=\sum_Py+\sum_By .$$
With no background this is the standard error of $\langle y\rangle_P$, and for $y\equiv1$ it gives exactly 1 with no error, so the full ring's $\langle w\rangle_S$ is handled with no special case. A Monte Carlo check (Poisson signal and background, sideband of equal population) confirmed both the mean and the error. Note the size of that error: about $\sigma_y\sqrt{S+2B}/S$, so the background costs precision as $\sqrt{1+2B/S}$, whatever estimator is used.

**$\kappa$ per region.** `Canvas_KappaByRegion_<Proxy>`, also in `Corrections/CheapSigExtract/<Family>/`, shows the raw $\kappa_{\rm eff}$ and $\kappa_z$ of the peak and of the sideband for every data-like sample, each with the ratio-of-means error. Where the two disagree, the background does not share the signal's response, which is why $P$ divides by the signal-region value and not by either of these.

**$\alpha$ is the family's own**, from the same decay constants as the consumer ($\alpha_\Lambda=0.749$, $\alpha_{\bar\Lambda}=-0.758$). `BothHyperons` mixes the two while the exact term is per species; it uses the mean of the two $\alpha^2$, which is within 1.2% of either, on a term that is itself a $\sim10\%$ correction, so the choice moves $P$ by $\lesssim0.2\%$.

Per family and proxy each folder holds two canvases: $\kappa_{\rm eff}$ and $\kappa_z$ for every data-like sample, and $P^{\rm full}$ and $P^z$ for every pair. The full ring is always drawn; $R_z$ joins it wherever its files exist.

**Errors.** They are first order, with the full derivative of the exact inversion.
- The data and fake rings come from separate consumer runs and are added in quadrature. The mixed events reuse the data's $\Lambda$s, so they over-cover slightly, as everywhere in this macro.
- The error of $3\langle u^2\rangle$ is propagated, but its covariance with $R_d$ is not, because there is no $\langle Ru^2\rangle$ profile. This is safe: every term in it carries $\partial P/\partial D=-P/D$, so the dropped term is at most $\sim2P\sigma_D/\sigma_R\sim10^{-3}$ of the leading one.
- $\kappa$ itself is only displayed, never used in $P$. The raw $\kappa$ keeps the covariance of $\langle u^2\rangle$ and $\langle w\rangle$ (the ratio-of-means formula above). The signal-extracted $\kappa_z$ takes their errors in quadrature: the two are positively correlated, both scaling as $n_z^2$, so this slightly overstates it.

### Output layout

Per family ($\Lambda$, $\bar\Lambda$, both), in this order:

| Folder | What | Sweeps variations? |
|---|---|---|
| `EtaProxy/`, `EtaSplitStudies/`, `EtaV0/` | the main $\eta$ dependences | yes |
| `AEE/` | $\langle R\rangle$ vs $\phi_\Lambda - \phi_p^*$, plus `Mass_Selection/` per observable | yes / **no** |
| `AEE_SignalExtracted/` | the same axis after signal extraction, plus `CrossSystem/` and `QA/` | yes |
| `HEE/` | $\langle R\rangle$ vs $\cos\theta_{\rm HEE}$, and the same split by the sign of the AEE angle | **no** |
| `CheapSigExtract/<proxy>/` | integrated $\langle R\rangle$ in and out of the mass peak, and their difference | yes |
| `IntegratedSummary/<proxy>/` | the integrated value per proxy, all variations on one axis | yes |
| `SignalExtraction/<proxy>/` | $\langle R\rangle_{\rm meas}$, $\langle R\rangle_S$, $\langle R\rangle_B$ side by side and across systematics | yes |
| `ProxyPt/`, `PVz/`, `QA_Mult/`, `OtherAngDepndncs/`, `BruteForce/` | the rest | yes |

> **TODO.** Two blocks still read only `ConsumerResults_<dataSuffix>.root` and draw Data alone:
> `HEE/` and the `Mass_Selection/` canvases under `AEE/`. That wastes the one thing the variation
> machinery is for. Both are differential rather than categorical, so they want the variations
> overlaid through `DrawComparisonCanvas` rather than the integrated drawer. (Actually not a big "TODO": these are postponed because their plots would be really messy and wouldn't really contribute to the bigger picture right now)

Blocks that need a specific position in the file **must run from inside the observable loop**:
directory order in a ROOT file is *creation* order, so a block written after the loop lands after
every group regardless of where its table entry sits. `runEtaSplitStudies()` and `runAeeHeeBlocks()`
are hooked that way.

### Naming

Consumer object names encode their ROOT type and a long "what is this" prefix, which produced folders
like `Lambda/EtaProxy/pRingObservableEtaLeadP` -- most of it repeated from the group folder above.
`DeriveObservableName` strips mechanically:

1. drop the type marker, `p2d` or `p`;
2. drop a leading `RingObservable`, **only** if something remains and it does not start with `Vs`.

That guard is what keeps `pRingVsCentrality` readable as `RingVsCentrality` instead of the meaningless
`VsCentrality`. Where the derived name is still unwieldy -- the AEE observables especially -- a
`displayName` field on `ProfileConfig` or `CategoricalObservable` overrides it.

> **Fix a bad name with an override, not by extending the rule.** The rule is meant to stay small
> enough to hold in your head; a pile of special cases in it is worse than nine explicit strings in
> the table, where the judgement is visible and reviewable.

Any new field on those structs goes at the **end**, since the tables use positional aggregate
initialisation.

### Which variations appear where

Three flags on `VariationConfig` decide membership, so adding a variation to a canvas is a one-line
table edit: `inRedux`, `inPtGateStudy`, `inEtaGateStudy`.

The full systematics canvas had become unreadable -- nine curves plus data, several of them
gated/ungated pairs whose interest is a comparison against a specific *partner* rather than against
the data. Those move to focused canvases and drop off the reduced one. Data-like is also off it: it
is the crudest of the mixing proxies, cheaper even than prevJet, and does not carry the correlations
between the angular distributions in full.

The two gate studies have **opposite polarity**, which is the thing most likely to be misread later.
The $p_{\rm T}$ study asks whether *adding* a gate changes an artificial proxy. The $\eta$ study asks
what is lost by *removing* one, since $\eta$ gating is now the default -- an invented direction that
ignores the experimental acceptance produces an observable the data could never have produced. So the
pair is (gated parent, ungated variation), and `Data - Var` reads as the cost of not gating. Only
RandJet and PerpToJet exist in both forms; the data-like proxy gates on $\eta$ by construction.

> **Every canvas drawn from the full variation set must emit its redux twin from the same call site.**
> Emitting them apart is how two loops end up disagreeing about which variations they show -- exactly
> the drift that once let an all-V0s canvas into the $\Lambda$ folder. `ReduxSubset()` is the single
> definition of the reduced set.

### Integrated values: read them, do not recompute them

The consumer already publishes the integrated $\langle R\rangle$ in **bin 1** of the `TProfile1D`s
`IntegratedCuts/pRingCuts`, `pRingCutsLeadingP` and `pRingCutsSubLeadingJet`. The Toy Model publishes
its own in `WithEtaGate/BothCuts/All/pRingProxyJet`, a single-bin `TProfile`. Those are authoritative
and are read directly.

> Note `pRingProxyJet`, **not** `pRingProxy`: the latter uses $\hat z$ as the axis rather than a jet
> direction, and is a different observable entirely.

This macro used to re-integrate each differential profile with `GetIntegratedProfile`, once per
observable folder, which recomputed a number that already existed by a different route -- a recipe
for two values quietly disagreeing after an innocent rebin. Those per-observable folders have been
removed in favour of one `IntegratedSummary/` per family.

The two routes can legitimately differ in one way: re-integrating covers only the histogram's bins,
so anything in the **underflow or overflow is silently excluded**, while the consumer's value counts
every candidate. Angular axes are safe by construction -- they cover the full interval -- but the
proxy-$p_{\rm T}$ and `EtaV0` axes do overflow, confirmed, with a small shift in some fake-ring
estimators (RandJet without the $\eta$ gate among them). A small unexplained difference on the
centrality axis was also observed. In all of these the consumer's value is the one to report.

`GetIntegratedProfile` survives only for the brute-force and cross-family canvases, which put the
*observable* on the axis rather than the variation.

> **TODO (needs the consumer's code in hand):** confirm every reported observable's integral is
> available from one of the `IntegratedCuts` X bins, so that no observable needs a per-observable
> integral at all.

### Reading the signal-extraction results

`sigExtractDir` is the seventh argument, defaulting to `<consumerDir>/../results_SigExtract`; pass
`"none"` to skip those plots entirely. Two things every fetch must respect:

- **Check `hExtractionStatus` bin 1 first.** A failed extraction writes no value histogram, but a
  zero-filled one would be indistinguishable from a real measurement of zero.
- **Never re-derive a stored difference.** `hRSigMinusRBkg`, `hIntegratedDiff*` and `hAEE_*` are
  computed where the primitives and their covariance are in scope. Differences between *systematic
  variations* are a different matter -- those come from independent consumer runs, so quadrature there
  is correct, and that is what the systematics canvases do.

$\langle R\rangle_{\rm measured}$ is **not** the same number as the "Data" column: the summary
framework averages over the consumer profile's whole mass range, while $\langle R\rangle_{\rm meas}$
is restricted to the achieved peak window. That is why it deliberately gets no systematics canvas of
its own -- two views of the same candidates through different code paths should not be read as
independent measurements.

### Reading the AEE signal extraction

`AEE_SignalExtracted/<proxy>_<species>/` carries the per-bin extraction against
$\phi_\Lambda - \phi_p^*$: `Canvas_SigVsBkg` (and its pull), `Canvas_SigMinusBkg`, the two
angle-combined digests, and `CrossSystem/` with each quantity across every consumer variation.

> **These objects live at the TOP LEVEL of the extraction file, not under a cut folder**, because the
> consumer books the source profiles with a bare path rather than `(folder + "/...")`. A path built
> the way the per-proxy `IntegratedSummary` ones are built finds nothing -- and fails *silently*,
> since every fetch is null-guarded. That is why nothing here uses `cutFolder`.

`CrossSystem/` reads bins 2, 3 and 4 of `hCombinedSummary`, so it is **coupled to the row order** of
that histogram in `signalExtractionRing.cxx`. Changing the order there without changing it here
produces plausible-looking plots of the wrong quantity.

`Canvas_SigMinusBkg` is read from `hRSigMinusRBkg`, never recomputed: $\langle R\rangle_S$ and
$\langle R\rangle_B$ share the sideband primitives, and quadrature understates their difference by
roughly a fifth.

### The cheap cross-checks

`CheapSigExtract/<proxy>/` and the `Mass_Selection/` canvases answer "how far does the background
pull the number?" without fitting anything. The former reads
`IntegratedCuts/p2dRingCuts*V0MassPeak` at X bin 1 ("All Lambda") crossed with the mass flag, where
**Y bin 1 is strictly out of the mass peak and Y bin 2 strictly in it** -- bin numbers, not fill
values, since the consumer fills 0 and 1 on a `{2, 0, 2}` axis. It produces the two states as
separate redux-style axes, plus `Canvas_InVsOutOfMassPeak` and `Canvas_InMinusOutOfMassPeak`.

> **That difference IS plain quadrature**, unlike the differences inside the signal extraction. In-
> and out-of-peak candidates occupy disjoint mass regions, so the two measurements share no
> candidate and are independent. The two kinds of difference now sit close together in the output,
> so the distinction is worth carrying: quadrature is correct between independent *samples* and
> wrong between correlated *functions of the same primitives*. They are **deliberately independent** of `signalExtractionRing`,
so that a disagreement between the two is informative rather than circular.

Inside an AEE observable's `Mass_Selection/`, two canvases with different jobs:

- `Canvas_MassSelection` -- in-peak against out-of-peak, built from the `massVariations` files
  (`_excludeOutOfPeak` / `_excludeInPeak`), and carrying a third curve, the unselected data.
- `Canvas_SidebandAgreement` -- left sideband against right, which the `massVariations` route cannot
  produce since it only knows in and out. **Read this one first**: a sideband subtraction is only
  meaningful if the two sides agree, so this is a precondition for trusting anything above it.

The two also read *different inputs* -- separate mass-cut consumer files versus a slice of the
three-bin mass axis in one nominal file. They should agree, and a disagreement is worth knowing about.
All mass-region canvases draw their colours from one set of constants so the legend need not be
re-learned between them.

### Presentation

Curves are drawn as zero-x-width point graphs with a small per-curve horizontal offset, not as
superposed histograms. ROOT's default `ErrorX` gives every `"PE"` histogram a horizontal bar half a
bin wide, which carries no information the axis does not already show, and eight curves put eight
sets of them inside every bin with all markers at identical $x$. The offset index is **shared between
the upper and lower pads**, so a curve sits at the same $x$ in both and can be followed down.

Categorical points are placed by index within their axis, so a filtered subset must be **rebuilt**
rather than have entries dropped -- dropping would leave gaps where the removed columns used to be.

Species-specific observables appear only on the family they describe (`SpeciesMatchesFamily`). A
per-species profile read from a single-species file is either empty or a duplicate of that species'
own curve. `BothHyperons` keeps all three, since there the mixed "Lambda-like" set is the QA for
species competition: with balanced yields its dependence should cancel, so a residual measures the
imbalance rather than a physics asymmetry.

### Reference paths

`MC_REF_DIR`, `PP_REF_DIR` and `TOY_MODEL_PATH` are **hardcoded at the top of
`run_all_wagons.sh`**, not discovered. Update them there if the local storage layout changes. Set any of them to `""` to drop that overlay.

---

## `zvtxBitForensics.cxx`

**Step 7**, once per wagon, run in parallel across `FORENSICS_JOBS` workers. Output: `results_consumer/zvtxBitForensics.root`.

Bit-level and integrity QA on the **raw derived AO2Ds**, before any O2 abstraction. It validates the inputs the event-mixing machinery depends on -- storage truncation, repeated stored values, duplicated collision rows, index-column integrity and mixing-bin occupancy -- reading raw ROOT with
no O2 headers, so its verdict is independent of the code it is validating.

It needs the O2Physics environment (`alienv enter O2Physics/latest`), because it uses `hadd` to merge its parallel batch outputs.

It is a genuinely different kind of tool from the rest -- a diagnostic rather than a measurement, with its own statistical model and a two-mode worker/finalize structure -- and is documented at length in the [appendix](#appendix----ao2d-bit-forensics-in-detail) below.

---

# Appendix -- AO2D bit forensics in detail

**File:** `zvtxBitForensics.cxx` -- **Step 7** of `run_all_wagons.sh`, once per wagon.

## 1. Why this exists

The consumer task uses `SameKindPair` event mixing to build a baseline in which the jet proxy of one collision is borrowed by another. That baseline is only meaningful if the borrowed proxy is genuinely uncorrelated with the target collision. Two failure modes would quietly break it:

1. **Self-correlation through indexing** -- a collision borrowing from itself, or an off-by-one in the relational index columns, so that jets/V0s are associated with the wrong collision.
2. **Self-correlation through duplicated rows** -- the same physical collision written twice into one dataframe under two different `globalIndex` values. The consumer's self-index guard cannot see this: the indices genuinely differ, the delta-index histogram reports a healthy separation, and the borrowed proxy is nonetheless from the same event.

Both are invisible from inside the task. This tool checks them from raw ROOT, with no O2 headers and no framework assumptions, so its verdict is independent of the code it is validating.

The investigation started from a different observation: **bit-identical `fZvtx` values appear in the AO2Ds**, and also bit-identical `fJetPt`. The primary vertex z is a fit output over many trajectories and has no detector-level discretisation, so the naive expectation is that exact repeats should essentially never occur. **That expectation is wrong, for two compounding reasons.**

## 2. Why bit-identical values are expected (i.e, have you heard of the "Birthday Paradox"?)

### 2.1 The stored grid is far coarser than float32

Between the vertex fit and the derived table sits the AO2D writer, which applies deliberate lossy compression: `truncateFloatFraction(value, mask)` performs a bitwise AND on the IEEE-754 word, zeroing the low mantissa bits. Different columns get different masks.

> **Unverified constant.** The collision position is probably using the `0xFFFFFFF0` mask (4 low bits cleared, 19 mantissa bits kept) and track `1/pt` something far more aggressive such as
> `0xFFFFFC00` (10 cleared, 13 kept). These specific hex values have **not** been confirmed against the O2 source yet, but you can find similar masking for them under "AliceO2/Detectors/AOD/src/AODMcProducerHelpers.cxx", "AliceO2/Common/MathUtils/include/MathUtils/detail/TypeTruncation.h" and "AliceO2/Detectors/PHOS/workflow/src/StandaloneAODProducerSpec.cxx". 

I myself couldn't really find much more info on those, so me and Claude just built a tool that measures the surviving bit count empirically and does not rely on them (thanks, Claudinho!).

With 19 surviving mantissa bits the relative step would be `2^-19 ~ 1.9e-6`, so when around `|z| = 5 cm` the absolute grid step is roughly `7.6e-6 cm`. That is a real, finite grid: about `10^6` usable cells across a +-10 cm acceptance, rather than the `~10^9` untruncated float32 would provide.

### 2.2 The birthday paradox is stronger than intuition suggests

The relevant quantity is not the probability that a *given* pair collides, but the number of pairs, which grows as `N^2`. For a smooth density `p(z)` sampled on a grid of local step `q(z)` (which is essentially the case of the FT0M centrality distributions, so this applies directly to them):

```
E[colliding pairs] = C(N,2) * INTEGRAL p(z)^2 q(z) dz
                   ~ ((N-1)/2) * SUM_i p(z_i) q(z_i)
```

For a Gaussian of `sigma ~ 5 cm` (`p ~ 0.08 cm^-1` at the peak) and `q ~ 7.6e-6 cm`, the per-pair collision probability is `~6e-7`. With `N = 1e5` collisions that is `~5e9` pairs and therefore
**~3000 expected duplicate pairs**; at `N = 1e6` it is `~3e5`.

**Seeing many repeated `fZvtx` values is the _null_ hypothesis, not the anomaly!** The formula above is what turns "this looks weird" into a quantitative test, and it is exactly what the `Counters/` directory implements.

### 2.3 What this does *not* test

Repeated `fZvtx` values are nearly orthogonal to mixing safety. Mixing bins on `axisPVz`, whose bin width is `1/3 cm`; two collisions agreeing to the last bit land in the same bin as two collisions differing by `1e-6 cm`. Bit-identity carries no extra information about whether mixing is correct. The checks that *do* bear on mixing are the `Fingerprints/`, `Integrity/` and `MixingPool/` modules.

In other words, this is just a measurement of how much the birthday paradox can affect us, not a QA of the actual event mixing engine by ALICE in O2Physics code (which we obviously assume is correct, up to an error on my dev hands, not the O2Group's dev team!).

## 3. Running it

Automatically, as Step 7 of the coordinator, is probably the most guaranteed form:

```bash
./run_all_wagons.sh                    # includes forensics
./run_all_wagons.sh --skip-forensics   # skips Step 7
```

Manually, in worker + finalize form (I myself never ran it like this, but this be working as is described below!):

```bash
# One manifest per batch: one absolute AO2D path per line.
ls -d $PWD/AO2Ds/AO2D*.root > batch_0.txt
zvtxBitForensics.exe batch_0.txt zvtxForensics_batch_0.root
# ... repeat per batch, in parallel ...
hadd -f zvtxBitForensics.root zvtxForensics_batch_*.root
zvtxBitForensics.exe --finalize zvtxBitForensics.root
```

Output lands in `results_consumer/zvtxBitForensics.root`, logs in `results_consumer/logs/zvtxBitForensics.log`, matching the `auxiliarySummaryPlots` convention. The manifests for each batch and the partial `.root` files live in a temporary `results_consumer/.forensics_batches_<PID>/` folder that is removed on completion **and** on interrupt.

> **This requires `hadd`**, i.e. the O2Physics environment (`alienv enter O2Physics/latest`).

### 3.1 Parallelism

`FORENSICS_JOBS` in `run_all_wagons.sh` defaults to **24**. Only a handful of branches are read per file (`SetBranchStatus("*", 0)`), so per-file cost is basket decompression rather than streaming whole AO2Ds. 24 saturates a SATA SSD and leaves headroom on NVMe. Raising it further only helps if the step proves CPU-bound rather than I/O-bound. Files are distributed **round-robin**, not in contiguous chunks, so that variation in AO2D size does not load a single worker.

### 3.2 Why two modes, and why `--finalize` is mandatory

**A ratio is not additive under `hadd`.** Merging two files each holding `obs/exp` would sum the ratios, which is meaningless. Every ratio in this tool is therefore stored by the workers as a separate `(numerator, denominator)` pair; `hadd` sums each independently; `--finalize` reopens the merged file and performs the division, writing the result back in-place.

The birthday expectation is quadratic in `N` **within a dataframe**, so it is computed per-DF and accumulated. Since each worker owns whole files (hence whole dataframes), summing per-DF contributions across workers is exactly correct.

> **If you add a histogram to this tool, it must be additive.** Never write a ratio, mean, or
> fraction from a worker. Write its numerator and denominator, and divide in `runFinalize()`.

### 3.3 Why not `TChain`

Chaining `DF_x/O2ringcollision` with `DF_y/O2ringcollision` erases the dataframe boundary -- and the dataframe boundary *is* the scope of event mixing, since `SameKindPair` never crosses it. A duplicate pair spanning two dataframes is harmless; the identical pair inside one dataframe is the one that can contaminate the baseline. Chaining would merge those two populations and report a number answering neither question. The tool iterates `DF_*` directories explicitly and treats each as its own unit; parallelism comes from distributing whole *files*, which is the natural I/O shard.

## 4. Output structure

### `Truncation/`

| Histogram | What it measures | How to read it |
|---|---|---|
| `hTrailingZeros_<col>` | Distribution of trailing zero mantissa bits | The **lowest populated bin** is the truncation depth: that many low bits were cleared by the writer. `23 - that` is the surviving precision. |
| `hTrailingZerosVsBinade_<col>` | Truncation depth vs `floor(log2\|value\|)` | A **horizontal band** means relative (mantissa-mask) truncation: constant significant bits, absolute step scaling with magnitude -- so a single "precision in cm" number would be misleading. A band **sloping with the binade** means fixed absolute rounding. |

### `Duplicates/`

| Histogram | What it measures | How to read it |
|---|---|---|
| `hMultiplicity_<col>` | How many rows share each distinct stored value | Long tail is expected wherever the density is high. You can compare its integral against `Counters/` rather than eyeballing it. |
| `hLogSpacing_<col>` | `log10` of the gap between adjacent distinct values | Spans many decades because sparse tails leave many empty grid cells between populated ones. Mostly superseded by the next row. |
| `hSpacingOverGridStep_<col>` | Gap divided by the **local** grid step | The strong evidence. Dividing out the magnitude dependence makes a true quantisation grid appear as a **picket fence of peaks at 1, 2, 3, ...**. If that fence is there, the grid is confirmed. A smooth continuum with nothing at small integers falsifies the grid hypothesis. |
| `hDuplicateRowGap_<col>` | `log10` of the row separation between **consecutive** value-sharing entries, in table order | Broad and featureless for chance coincidences. A **spike at small gaps** points at split vertices or double-written rows. See [section 6.1](#61-why-row-gaps-are-consecutive-rather-than-all-pairs) for why this is not an all-pairs enumeration. |
| `pMultiplicityVsValue_<col>` | Mean multiplicity vs value | Should track `p(z) * q(z)`: highest where the distribution is dense. |

### `Counters/` (additive; ratios appear after `--finalize`)

| Histogram | Meaning |
|---|---|
| `hEntriesPerColumn`, `hDistinctPerColumn` | Rows read and distinct values, per column |
| `hObservedDuplicatePairs`, `hExpectedDuplicatePairs` | Observed vs birthday-model pair counts, per column |
| `hNonFinitePerColumn` | NaN/inf guard. **Should be empty.** |
| `hConstantColumnDataframes` | Dataframes in which a column took a single value -- such a column contributes nothing to the fingerprint test |
| `hObservedPairsVsDfSize`, `hExpectedPairsVsDfSize` | The same counters for `Zvtx`, binned in `log10(collisions per dataframe)` |
| `hDuplicatePairRatio` *(finalize)* | **The headline number.** Ratio near 1 = pure numerical coincidence, expected and harmless. Large excess = suspect duplicated rows. Large deficit = the density model or the measured truncation is wrong. |
| `hDistinctFraction` *(finalize)* | Fraction of rows carrying a distinct value, per column |
| `hDuplicatePairRatioVsDfSize` *(finalize)* | **Should be flat at 1** across dataframe sizes. This is a nontrivial confirmation: observed and expected each vary by orders of magnitude across the range while their ratio does not. A **slope** means the quadratic birthday law is wrong and something other than chance generates the repeats. |

### `Fingerprints/`

A "fingerprint" is the concatenated bit patterns of the stored collision columns, in the order
`Zvtx, CentFT0M, CentFT0C, CentFV0A, InteractionRate`. The test is **cumulative**: bin `k` counts
collision pairs agreeing bit-for-bit on the *first k* columns.

| Histogram | How to read it |
|---|---|
| `hObservedFingerprintPairs` | **Read the shape, not the value.** Chance coincidences fall steeply with each added column. A genuinely duplicated row matches on *everything*, so real duplicates make the curve **plateau at a nonzero floor**. Steady falloff = healthy. Plateau = duplicated collision rows. |
| `hExpectedFingerprintPairs` | Chance expectation assuming column independence. The three centrality columns are all multiplicity-derived and therefore **mutually correlated**, so this product is a **lower bound** on the true chance rate, never an upper bound. Treat it as a reference curve, not a threshold. |
| `hFingerprintGroupSize` | Size of each full-fingerprint group. **Should be empty** on healthy data. |
| `hFingerprintRowGap` | Row separation between consecutive members of a full-fingerprint group, linear from 0. A spike at `\|dRow\| = 1` is consecutive duplicated rows -- the signature of a split vertex or a double write. |

The cumulative construction is what removes the need for an independence model: the plateau-versus-
falloff reading is assumption-free, which matters precisely because the centrality columns are
correlated in an uncontrolled way.

**Scope:** computed **per dataframe only**. Cross-dataframe repeats cannot reach the mixing pool,
and restricting to per-DF also removes the problem of merging a global hash map across workers.

### `Integrity/`

| Histogram | How to read it |
|---|---|
| `hIntegrityViolations` | Labelled counter. **Every bin should be zero.** Bins: index `< 0`; index `>= nCollisions`; index not sorted; more than one leadP row per collision; missing collision tree; missing indexed tree; non-finite value. |
| `hIndexRangeExcess` | Signed distance outside `[0, nCollisions)` per table (`0=jet, 1=leadP, 2=V0`). A stripe at `+1` or `-1` is the **off-by-one**, localised to a specific table. |
| `hRowsPerCollisionLeadP` | The producer writes **at most one** leading-particle row per collision, and the consumer's `getMixLeadPPt` relies on that by taking `rows.begin()`. Nothing in the data model enforces it. **Bin 2 and above must be empty.** |
| `hRowsPerCollisionJet`, `hRowsPerCollisionV0` | Characterisation of table occupancy |
| `hCollisionContentPattern` | Which tables reference each collision (`jet` / `leadP` / `V0` bitmask) |

> **The "empty" bin of `hCollisionContentPattern` will be populated on healthy data.** The producer
> writes the collision row *before* the jet and leading-particle logic, which has early returns
> (`fjParticles.size() < 1`, `jets.empty()`) and a `minLeadParticlePt` gate. Orphan collisions
> therefore exist by construction. This is characterisation, not an alarm.

### `MixingPool/`

Bins collisions on the mixing grid `(Zvtx, proxy pT, centrality)` using the consumer's own axes, and reports what each collision could possibly mix with. (this needs to be updated if you really want "consumer Vs this script" matching, if the axes are ever updated)

| Histogram | How to read it |
|---|---|
| `hCollisionsPerDataframe` | Whether dataframes are large enough for mixing at all |
| `hMixBinOccupancyLeadP`, `hMixBinOccupancyLeadJet` | Collisions per mixing bin, per dataframe |
| `hMixPoolOutcomeLeadP`, `hMixPoolOutcomeLeadJet` | Three-way per-collision verdict: *no proxy*; *proxy, alone in bin*; *proxy, partner available*. **The middle bin is the useful number**: those collisions can never be mixed regardless of `mixedEventWindowSize`, so it is an upper bound on mixing efficiency computed without running the consumer. |
| `hZvtxAcceptance` | Collisions surviving the `\|Zvtx\| < 10 cm` filter |

## 5. Reading a run: checklist

1. `Counters/hDuplicatePairRatio` near 1 for every column -> repeats are chance, as predicted.
2. `Counters/hDuplicatePairRatioVsDfSize` flat -> the birthday model itself is validated.
3. `Fingerprints/hObservedFingerprintPairs` falling steeply with no plateau -> no duplicated rows.
4. `Integrity/hIntegrityViolations` entirely empty -> indices and the one-leadP-per-collision
   invariant both hold.
5. `Duplicates/hSpacingOverGridStep_Zvtx` showing the integer picket fence -> the storage grid is
   confirmed and the measured truncation depth is trustworthy.
6. `MixingPool/hMixPoolOutcome*` -> how much of the sample is mixable at all.

If 1-4 all pass, the AO2Ds are clean and any residual mixing problem lives in the consumer, not in
its input.

### 5.1 Progress and logging

Each worker writes a timestamped progress log to
`results_consumer/logs/forensics_logs/batch_<N>.log`, flushed line by line so the
file is always current and `tail -f` works during a run:

```
[0:00:00] Worker started.
[0:00:00]   files    : 12
[0:00:00] [1/12] AO2D_7.root
[0:00:00]     DF 1/8 DF_2300474821381920 ...
[0:00:09]     DF 1 done: 9812 coll, 41233 jets, 8901 leadP, 22104 V0 | read 0.31s cols 7.94s fprint 0.28s integ 0.44s proxy 0.19s pool 0.03s | total 9.19s
[0:01:14] [1/12] done in 74.2s (78411 collisions) | elapsed 0:01:14 | ETA 0:13:36
```

The `DF n/m ...` line is printed **before** that dataframe is processed, so if a
worker hangs, the last line of its log names what it hung on. The per-stage
timings double as a profiler: `cols` covers the truncation and duplicate
analysis, `fprint` the fingerprint scan, `integ` the index validation, `proxy`
the jet and leading-particle columns, `pool` the mixing-bin occupancy.

Per-file lines carry a running ETA extrapolated from the mean time per file so
far. AO2D sizes vary, so treat it as an order of magnitude rather than a
deadline.

## 6. Design constraints worth preserving

### 6.1 Why row gaps are consecutive rather than all-pairs

The row-separation histograms record the gap between **consecutive** members of
a value-sharing group after sorting that group by row, not the gap of every
pair inside it. There are two independent reasons, and both matter.

**Cost.** All-pairs enumeration is `O(m^2)` in the group size. `fInteractionRate`
is a timeframe-level quantity and is effectively constant within a dataframe, so
it forms a single group containing every collision. At 10k collisions per
dataframe that is `5e7` pairs, each with a `log10` and a histogram fill, for one
column in one dataframe. Measured against the consecutive-gap alternative:

| group size | all pairs | consecutive gaps | speedup |
| --- | --- | --- | --- |
| 1 000 | 3.7 ms | 0.018 ms | ~200x |
| 5 000 | 86 ms | 0.086 ms | ~1000x |
| 10 000 | 339 ms | 0.184 ms | ~1800x |

(Benchmarked with a plain array increment in place of `TH1D::Fill`, so the
all-pairs column is a lower bound on its real cost.)

**Statistics.** The question the histogram asks is whether value-sharing entries
sit *adjacent* in the table. Consecutive gaps answer that directly. All-pairs is
dominated by the trivial combinatorics of large groups -- a group of `m` entries
contributes `m(m-1)/2` mostly-large separations -- which buries the small-gap
signal the plot exists to expose. The linear version is therefore the sharper
statistic as well as the affordable one.

> **This changed the contents of `hDuplicateRowGap_*` and `hFingerprintRowGap`.**
> Their normalisation and shape are not comparable with output produced before
> this change.

### 6.2 Other invariants

- **Additivity.** Workers write only summable quantities. Ratios exist only after `--finalize`.
- **Fixed axes.** Every histogram axis is fixed at book time. Auto-ranging (`TH1(name, title, n,
  0.0, 0.0)`) would give each worker a different axis and make `hadd` unmergeable *silently* --
  hence the explicit `axisLow`/`axisHigh` fields in `ColumnSpec`.
- **Per-dataframe scope.** Duplicate and fingerprint statistics are computed within a dataframe,
  matching the scope of `SameKindPair`.
- **No error bars on ratios.** `writeRatio()` deliberately sets zero errors: numerator and
  denominator are not counts of the same trials, so binomial errors would be wrong. Absent beats
  misleading.
- **Truncation measured per dataframe.** The grid step is needed during the same pass that uses it.
  With thousands of entries per dataframe the minimum saturates; the pooled distribution remains
  recoverable from `hTrailingZeros` after merging.
- **Each tree is read once.** `scanIndexedTree()` returns the full proxy branch
  alongside the per-collision counts, so `O2ringjet` and `O2ringleadp` are not
  reopened for the column-wise analysis. Do not reintroduce a second read.
- **`__builtin_ctz` for the trailing-zero count.** It compiles to a single TZCNT
  and is bit-identical to the shift loop it replaced, which needed up to 23
  iterations.
- **Fingerprint keys stay exact.** The cumulative scan keys on the concatenated
  bit patterns themselves (`std::map<std::vector<uint32_t>, ...>`), not on a
  hash. A 64-bit hash would be roughly twice as fast but would admit a small
  false-match probability, and a tool whose entire purpose is exactness about
  bit patterns should not carry one.

## 7. Known limitations and TODOs

- **`TODO:`** The mixing axes (`kAxisPVzLow/High/Bins`, `kAxisJetPtEdges`, `kAxisCentralityEdges`)
  are hardcoded and **must be kept in sync with `axisConfigurations`** in
  `lambdaJetPolarizationIonsDerived.cxx`. They should instead be read from the dpl-config JSON.
- **`TODO:`** `MixingPool/` currently reports bin occupancy only. Extending it to reproduce the
  `SameKindPair` sliding window would let the predicted candidate distribution be compared bin by
  bin against the consumer's own
  `EventMixingQA/CollLoopOutcome/hMixedEventLeadPCandidates` -- two independent measurements of the
  same quantity, one from raw ROOT and one from inside the task. Disagreement would localise a fault
  in the `FlexibleBinningPolicy` setup without reading framework source.
- A sliding-window scan of *value equality* was considered and rejected: whether two collisions share
  a stored value does not depend on their separation in the table, so sweeping window sizes would
  re-measure the same rate with worse statistics. The window sweep is only informative when applied
  to **bin occupancy**, which is what the TODO above proposes.
- The O2 truncation masks quoted in section 2.1 remain **unverified** against O2 source.