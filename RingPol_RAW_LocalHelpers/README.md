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
4. [Macro reference](#macro-reference)
5. [Appendix -- AO2D bit forensics in detail](#appendix----ao2d-bit-forensics-in-detail)

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
| `signalExtractionRing.cxx` | Sideband signal extraction of $\langle R\rangle_S$ from the $\Lambda$ invariant-mass spectrum |
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

The two per-wagon macros differ: `auxiliarySummaryPlots.exe` takes `<consumerDir> [mcRefDir] [ppRefDir] [toyModelPath] [cutFolder]`, and `zvtxBitForensics.exe` takes a manifest and an output file (see its [appendix](#appendix----ao2d-bit-forensics-in-detail)).

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

# Macro reference

## `signalExtractionRing.cxx`

**Step 3**, once per config. Output: `results_SigExtract/SignalExtraction_<SUFFIX>.root`.

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

Per bin of the differential observable ($\Delta\phi$, $\Delta\theta$, and 3D slices in
$p_{\rm T}^{\Lambda}$ and leading-jet $p_{\rm T}$):

1. Project the mass spectrum and fit the peak for $\mu$ and $\sigma$.
2. **Peak window** $[\mu - 4\sigma,\ \mu + 4\sigma]$, chosen to match the sigma0 analysis from Gianni, because it seemed very reasonable an analysis!
3. **Sidebands** fitted beyond $\mu \pm 6\sigma$, leaving a $2\sigma$ guard band on each side, to guarantee we are not reaching for signal (still could be tuned over more complete analyses).
4. Fit a second-order polynomial to the sidebands -- separately for the counts and for the ring numerator as a function of mass -- and integrate each over the peak window. In other words, **bin-counting**, then background subtraction via integrated background function (much safer than just fitting everything together! I tried it and the fits are honestly really bad -- it is all in a QA folder under the signal extraction's results, if you want to compare them).

Two quantities are extrapolated independently: the background **yield** $N_B$, and the background **ring numerator** $\sum_i r_i$, whose ratio is $\langle R\rangle_B$.

### Design constraints worth preserving

**The sideband fit uses `TGraphErrors`, not `TH1`.** The fit domain is *discontinuous* -- two
intervals with a gap between them. A `TH1` cannot express that: fitting over the full range would
include the peak, and a fit sub-range can only ever be one contiguous interval. A `TGraphErrors`
holding only the sideband points states the domain exactly.

**Sideband points are densities, not raw counts.** The mass axis is not uniformly binned, so raw
per-bin counts would hand the fit a spurious shape that follows the binning. Each point is divided by
its bin width; integrating over the peak window then recovers counts.

**Background uncertainty comes from `TF1::IntegralError`.** The three polynomial coefficients are
strongly correlated, so adding their individual `GetParError` values in quadrature would be wrong.
`IntegralError` uses the full covariance matrix from the fit result, which is the only correct way to
propagate an extrapolated integral.

**A `TProfile3D` must be converted before projecting.** `Project3D("yx")` on a `TProfile3D` returns a
`TProfile2D` whose every cell holds the *grand mean* over all $z$ bins in range -- not what a sliced
analysis wants, and with errors that are wrong for this purpose. `BuildNumFromProfile3D` therefore
converts to a `TH3D` first, in which

$$\text{bin content} = \bar R \cdot N = \sum_i r_i, \qquad
\text{bin error} = \frac{\sigma_R}{\sqrt N}\cdot N = \sigma_R\sqrt{N} ,$$

the absolute sum and its correct error, which *are* additive under projection. Do not reintroduce a
direct `Project3D` on the profile.

### Assumptions this rests on

1. Background polarization varies smoothly in invariant mass.
2. No strong mass--observable correlation exists (**must be validated**, not assumed).
3. Sidebands are representative of the background under the peak.

If these fail, the escape route is a simultaneous mass--polarization fit,

$$\text{Numerator}(m) = S(m)\,R_S + B(m)\,R_B, \qquad \text{Denominator}(m) = S(m) + B(m),$$

which avoids sideband extrapolation entirely. The `ResultsCombinedFit/` output directory is
groundwork for exactly that.

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

### Section 1 -- polarization vector fields

$\langle \vec{P}^{*}\rangle$ drawn over the $\Lambda$ momentum plane, in the same visual language as the Helicity Efficiency Toy Model's plotter ["PhD_codes/ToyModels/plotHelicityEfficiency.cxx"](../ToyModels/plotHelicityEfficiency.cxx), so that real data and toy model can be compared by eye.

| Panel | COLZ background | Arrows |
|---|---|---|
| X--Y plane | $\langle P^{*}_z\rangle$ (out of plane) | $(\langle P^{*}_x\rangle,\ \langle P^{*}_y\rangle)$ |
| Z--X plane | $\langle P^{*}_y\rangle$ (out of plane) | $(\langle P^{*}_z\rangle,\ \langle P^{*}_x\rangle)$ |

Arrows are **block-averaged** over `arrowBlockSize` x `arrowBlockSize` tiles of the underlying `TProfile2D`, because one arrow per histogram bin is unreadable (and too unstable in error!). Bins below `minEntries` are excluded
from the tile average rather than diluting it.

Arrow *length* is normalised to a **percentile** of the tile-magnitude distribution (95th by default) rather than to the maximum: a single noisy tile would otherwise set the scale and shrink every real arrow to nothing. Lengths are then capped at that reference, so an outlier shows up as a full-length arrow instead of one running off the pad, and the reference magnitude is printed on the panel so the scale never has to be guessed.

The COLZ range is symmetric about zero, so the sign of the out-of-plane component reads directly from the colour.

### Section 2 -- ring-observable 2D maps

$\langle R\rangle$ over the same two momentum planes, for both jet proxies (FastJet and leading particle). No arrow overlay: the ring observable is already one scalar per candidate, so the COLZ background *is* the whole result. No arrows here!

### On the duplicated `DrawVectorFieldPanel()`

It is copied almost verbatim from `plotHelicityEfficiency.cxx`, deliberately. The toy-model plotter and this post-processor are independent workflows, separately compiled and separately run, and the function is small, self-contained and stable. Sharing it through a common header would couple two otherwise unrelated build targets for very little. If it ever needs to diverge meaningfully, or grows non-trivially, that is the trigger to revisit though.

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

Beyond the overlays it produces subtracted curves against the data reference, folded versions where symmetry makes that meaningful, $\tanh$ fits of the $\eta$ dependence, mass signal-versus-background splits, and per-family integrated summaries with a cross-family canvas above them.

### Design constraints worth preserving

**Files are opened once.** Every fetch used to reopen and reclose the ROOT file in some previous versions of the code, so the number of opens scaled as (families x observables x variations) -- pure overhead (!), and a real cost once the observable list grew. `FileCache` opens each file once, keeps it open, and closes everything at the end. Failed opens are cached too, so a missing systematic warns **once** instead of once per observable.

**Objects handed out are detached clones.** `FetchClone` calls `SetDirectory(nullptr)` immediately, because `Clone()` otherwise attaches to `gDirectory` -- which at that point is the *output* file. That detachment is what lets the clones stay valid after `CloseAll()`.

**It reads one kinematic-cut family**, `Ring` by default, overridable with the fifth argument. `Ring` is the family the consumer always books, which is why it is the default.

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
| `hTrailingZerosVsBinade_<col>` | Truncation depth vs `floor(log2|value|)` | A **horizontal band** means relative (mantissa-mask) truncation: constant significant bits, absolute step scaling with magnitude -- so a single "precision in cm" number would be misleading. A band **sloping with the binade** means fixed absolute rounding. |

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
| `hFingerprintRowGap` | Row separation between consecutive members of a full-fingerprint group, linear from 0. A spike at `|dRow| = 1` is consecutive duplicated rows -- the signature of a split vertex or a double write. |

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
| `hZvtxAcceptance` | Collisions surviving the `|Zvtx| < 10 cm` filter |

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