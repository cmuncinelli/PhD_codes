# RingPol RAW Local Helpers

Local analysis framework for the **ring polarization observable**, operating on both raw AOD data
and Hyperloop-derived data outside of the Hyperloop environment. This folder collects the tools
needed to download, process, and analyze data locally -- complementing the O2Physics tasks that run
inside Hyperloop.

## Subfolders

### [`Local_framework/`](Local_framework)
Scripts and configuration for running the O2 DPL analysis workflow locally on raw AOD files.
Includes download helpers, path generators, and multi-threaded execution scripts.

### [`DerivedDataHY/`](DerivedDataHY)
Framework for downloading and analyzing Hyperloop-produced derived data locally.
See [`DerivedDataHY/README.md`](DerivedDataHY/README.md) for full documentation.

### [`JsonExamples/`](JsonExamples)
Example DPL workflow configuration files (`.json`) for various analysis scenarios: both-hyperon
mode, Lambda-only, AntiLambda-only, and full QA with permissive pT cuts.

## C++ Analysis Codes

(most of these are called by a coordinator in [`DerivedDataHY/run_all_wagons.sh`](DerivedDataHY/run_all_wagons.sh))

| File | Description |
|---|---|
| `signalExtractionRing.cxx` | Main signal extraction code for the ring polarization observable. Extracts the ring signal from invariant mass spectra across the (pT, centrality, Δφ) analysis bins. |
| `signalExtractionRingTest.cxx` | Test and development version of the signal extraction code. |
| `extractDeltaErrors.cxx` | Extracts and propagates statistical uncertainties using the Delta Method, which more accurately preserves the covariance terms to avoid SEM mathematical pitfalls. (Did not improve by much though) |
| `makeCumulativeDCAdauProfile.cxx` | Evaluates the robustness of the ring observable against varying Distance of Closest Approach (DCA) cuts, specially the DCA of the daughter particles wrt the PV. Condenses 2D/3D differential QA profiles into unified cumulative distributions. |
| `auxiliarySummaryPlots.cxx` | Cross-configuration aggregation macro. Combines baseline, MC reference, Helicity Toy Model inputs, and systematic variations (e.g., Data-like Jet, Rand Jet) of ring observables into appropriately scaled, unified canvases for easy comparison. |
| `auxiliaryPerConfigPlots.cxx` | Per-configuration post-processing macro. Goal is to analyze all configurations for a given macro and output useful plots which do not fit into O2Physics' principles. For instance, this code has the necessary functions to plot the polarization vector field graph inspired by the Helicity Efficiency Toy Model's plotting function. |
| `zvtxBitForensics.cxx` | Bit-level and integrity QA on the **raw derived AO2Ds**, before any O2 abstraction. Validates the inputs the event-mixing machinery depends on: storage truncation, repeated stored values, duplicated collision rows, index-column integrity, and mixing-bin occupancy. See [AO2D bit forensics](#ao2d-bit-forensics-and-event-mixing-input-qa) below. |
| `testVersionOfHeaders.h` | Development version of shared header definitions. |

## Utilities

| File | Description |
|---|---|
| `mergeAODDerived.sh` | Merges AOD derived data files prior to local analysis. |

---

# AO2D bit forensics and event-mixing input QA

**File:** `zvtxBitForensics.cxx` -- **Step 7** of `run_all_wagons.sh`, once per wagon.

## 1. Why this exists

The consumer task uses `SameKindPair` event mixing to build a baseline in which the jet proxy of
one collision is borrowed by another. That baseline is only meaningful if the borrowed proxy is
genuinely uncorrelated with the target collision. Two failure modes would quietly break it:

1. **Self-correlation through indexing** -- a collision borrowing from itself, or an off-by-one in
   the relational index columns, so that jets/V0s are associated with the wrong collision.
2. **Self-correlation through duplicated rows** -- the same physical collision written twice into
   one dataframe under two different `globalIndex` values. The consumer's self-index guard cannot
   see this: the indices genuinely differ, the delta-index histogram reports a healthy separation,
   and the borrowed proxy is nonetheless from the same event.

Both are invisible from inside the task. This tool checks them from raw ROOT, with no O2 headers
and no framework assumptions, so its verdict is independent of the code it is validating.

The investigation started from a different observation: **bit-identical `fZvtx` values appear in the
AO2Ds**, and also (more often) bit-identical `fJetPt`. The primary vertex z is a fit output over
many trajectories and has no detector-level discretisation, so the naive expectation is that exact
repeats should essentially never occur. That expectation is wrong, for two compounding reasons.

## 2. Why bit-identical values are expected

### 2.1 The stored grid is far coarser than float32

Between the vertex fit and the derived table sits the AO2D writer, which applies deliberate lossy
compression: `truncateFloatFraction(value, mask)` performs a bitwise AND on the IEEE-754 word,
zeroing the low mantissa bits. Different columns get different masks.

> **Unverified constant.** The collision position is believed to use `0xFFFFFFF0` (4 low bits
> cleared, 19 mantissa bits kept) and track `1/pt` something far more aggressive such as
> `0xFFFFFC00` (10 cleared, 13 kept). These specific hex values have **not** been confirmed against
> the O2 source. Grep your O2 checkout for `mCollisionPosition` and `truncateFloatFraction` in
> `AODProducerWorkflowSpec` before quoting them anywhere. The tool measures the surviving bit count
> empirically and does not rely on them.

With 19 surviving mantissa bits the relative step is `2^-19 ~ 1.9e-6`, so around `|z| = 5 cm` the
absolute grid step is roughly `7.6e-6 cm`. That is a real, finite grid: about `10^6` usable cells
across a +-10 cm acceptance, rather than the `~10^9` untruncated float32 would provide.

The heavier truncation applied to track momenta explains why `fJetPt` repeats even more often than
`fZvtx`: a single-constituent or single-hard-track jet inherits the coarse stored track momentum
almost directly.

### 2.2 The birthday paradox is stronger than intuition suggests

The relevant quantity is not the probability that a *given* pair collides, but the number of pairs,
which grows as `N^2`. For a smooth density `p(z)` sampled on a grid of local step `q(z)`:

```
E[colliding pairs] = C(N,2) * INTEGRAL p(z)^2 q(z) dz
                   ~ ((N-1)/2) * SUM_i p(z_i) q(z_i)
```

For a Gaussian of `sigma ~ 5 cm` (`p ~ 0.08 cm^-1` at the peak) and `q ~ 7.6e-6 cm`, the per-pair
collision probability is `~6e-7`. With `N = 1e5` collisions that is `~5e9` pairs and therefore
**~3000 expected duplicate pairs**; at `N = 1e6` it is `~3e5`.

**Seeing many repeated `fZvtx` values is the null hypothesis, not the anomaly.** The formula above
is what turns "this looks weird" into a quantitative test, and it is exactly what the `Counters/`
directory implements.

### 2.3 What this does *not* test

Repeated `fZvtx` values are nearly orthogonal to mixing safety. Mixing bins on `axisPVz`, whose bin
width is `1/3 cm`; two collisions agreeing to the last bit land in the same bin as two collisions
differing by `1e-6 cm`. Bit-identity carries no extra information about whether mixing is correct.
The checks that *do* bear on mixing are the `Fingerprints/`, `Integrity/` and `MixingPool/` modules.

## 3. Running it

Automatically, as Step 7 of the coordinator:

```bash
./run_all_wagons.sh                    # includes forensics
./run_all_wagons.sh --skip-forensics   # skips Step 7
```

Manually, in worker + finalize form:

```bash
# One manifest per batch: one absolute AO2D path per line.
ls -d $PWD/AO2Ds/AO2D*.root > batch_0.txt
zvtxBitForensics.exe batch_0.txt zvtxForensics_batch_0.root
# ... repeat per batch, in parallel ...
hadd -f zvtxBitForensics.root zvtxForensics_batch_*.root
zvtxBitForensics.exe --finalize zvtxBitForensics.root
```

Output lands in `results_consumer/zvtxBitForensics.root`, logs in
`results_consumer/logs/zvtxBitForensics.log`, matching the `auxiliarySummaryPlots` convention.
The batch manifests and partial `.root` files live in a temporary
`results_consumer/.forensics_batches_<PID>/` folder that is removed on completion **and** on
interrupt.

> **Requires `hadd`**, i.e. the O2Physics environment (`alienv enter O2Physics/latest`). This is the
> only step of `run_all_wagons.sh` that needs it.

### 3.1 Parallelism

`FORENSICS_JOBS` in `run_all_wagons.sh` defaults to **24**. Only a handful of branches are read per
file (`SetBranchStatus("*", 0)`), so per-file cost is basket decompression rather than streaming
whole AO2Ds. 24 saturates a SATA SSD and leaves headroom on NVMe. Raising it further only helps if
the step proves CPU-bound rather than I/O-bound. Files are distributed **round-robin**, not in
contiguous chunks, so that variation in AO2D size does not load a single worker.

### 3.2 Why two modes, and why `--finalize` is mandatory

**A ratio is not additive under `hadd`.** Merging two files each holding `obs/exp` would sum the
ratios, which is meaningless. Every ratio in this tool is therefore stored by the workers as a
separate `(numerator, denominator)` pair; `hadd` sums each independently; `--finalize` reopens the
merged file and performs the division, writing the result back in place.

The birthday expectation is quadratic in `N` **within a dataframe**, so it is computed per-DF and
accumulated. Since each worker owns whole files (hence whole dataframes), summing per-DF
contributions across workers is exactly correct.

> **If you add a histogram to this tool, it must be additive.** Never write a ratio, mean, or
> fraction from a worker. Write its numerator and denominator, and divide in `runFinalize()`.

### 3.3 Why not `TChain`

Chaining `DF_x/O2ringcollision` with `DF_y/O2ringcollision` erases the dataframe boundary -- and the
dataframe boundary *is* the scope of event mixing, since `SameKindPair` never crosses it. A
duplicate pair spanning two dataframes is harmless; the identical pair inside one dataframe is the
one that can contaminate the baseline. Chaining would merge those two populations and report a
number answering neither question. The tool iterates `DF_*` directories explicitly and treats each
as its own unit; parallelism comes from distributing whole *files*, which is the natural I/O shard.

## 4. Output structure

### `Truncation/`

| Histogram | What it measures | How to read it |
|---|---|---|
| `hTrailingZeros_<col>` | Distribution of trailing zero mantissa bits | The **lowest populated bin** is the truncation depth: that many low bits were cleared by the writer. `23 - that` is the surviving precision. |
| `hTrailingZerosVsBinade_<col>` | Truncation depth vs `floor(log2|value|)` | A **horizontal band** means relative (mantissa-mask) truncation: constant significant bits, absolute step scaling with magnitude -- so a single "precision in cm" number would be misleading. A band **sloping with the binade** means fixed absolute rounding. |

### `Duplicates/`

| Histogram | What it measures | How to read it |
|---|---|---|
| `hMultiplicity_<col>` | How many rows share each distinct stored value | Long tail is expected wherever the density is high. Compare its integral against `Counters/` rather than eyeballing it. |
| `hLogSpacing_<col>` | `log10` of the gap between adjacent distinct values | Spans many decades because sparse tails leave many empty grid cells between populated ones. Mostly superseded by the next row. |
| `hSpacingOverGridStep_<col>` | Gap divided by the **local** grid step | The strong evidence. Dividing out the magnitude dependence makes a true quantisation grid appear as a **picket fence of peaks at 1, 2, 3, ...**. If that fence is there, the grid is confirmed. A smooth continuum with nothing at small integers falsifies the grid hypothesis. |
| `hDuplicateRowGap_<col>` | `log10` of the row separation of value-sharing pairs | Broad and featureless for chance coincidences. A **spike at small gaps** points at split vertices or double-written rows. |
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
| `hFingerprintRowGap` | Row separation of full-fingerprint pairs, linear from 0. A spike at `|dRow| = 1` is consecutive duplicated rows -- the signature of a split vertex or a double write. |

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

Bins collisions on the mixing grid `(Zvtx, proxy pT, centrality)` using the consumer's own axes, and
reports what each collision could possibly mix with.

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

## 6. Design constraints worth preserving

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