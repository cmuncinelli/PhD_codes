# hydroTruthLab

Truth-level laboratory for the ring observable in hydrodynamics. The producer reads the per-cell Lambda polarization that iSS computes on the MUSIC freeze-out surface, one file per event, and evaluates the full 3D ring observable together with its longitudinal projection $R_z$ and its companions -- with no sampling and no decays. Hydro gives perfect PID, momentum, jet direction and polarization, so everything here is the truth that the data analysis (`useRingZ` in `lambdaJetPolarizationIonsDerived.cxx`) tries to reconstruct. The companion macro, PerConfig, combines the bullet and no_bullet samples of one configuration, computes the nonlinear quantities with their jackknife errors and builds the analogues of the figures of PhysRevC.109.014905 at the centrality of the dataset.

| File | Description |
|---|---|
| `hydroTruthLab.cxx` | Producer. Hydro `.dat` files of one sample --> `hydroTruthLab_<run>_<sample>.root`. |
| `hydroTruthLabPerConfig.cxx` | Downstream (PerConfig, formerly the summary). Both producer files of one run --> `hydroTruthLabPerConfig_<run>.root`. |
| `runHydroTruthLab.sh` | Bash coordinator: compiles both macros, runs them for every run of the registries, in parallel, and holds the defaults of every configurable (including the output folder). |
| `hydroTruthLab_datasets.conf` | Registry of the hydro datasets: samples, jet file and centrality of each. |
| `hydroTruthLab_runs.conf` | Registry of the runs (a dataset and the options it changes) and of the comparisons across runs. |
| `surveyHydroDatasets.sh` | Read-only survey of the hydro folders: samples, events, grid, MUSIC inputs and the bullet - no_bullet pairing, in a few lines per dataset. |
| `README.md` | This file: design and "how to run me". Implementation details live in the code comments. |

A run is one dataset with one configuration, labelled `<dataset>_<variation>` (`lumpy4050_default`, `lumpy4050_smudge0.1`). Comparisons across runs (centrality, smudge width, ...) will be drawn by `hydroTruthLabSummaryPlots.cxx`, not written yet (see Known issues and TODO).

Both macros are single-threaded (no OpenMP, no NUMA binding) and use no fits.

## IMPORTANT: the `S^mu` columns hold `P^mu`, not `S^mu`

> **WARNING -- there *intentionally* is no factor of 2 anywhere in this code.** The columns labelled `S^t S^x S^y S^z` and `S^t_LRF S^x_LRF S^y_LRF S^z_LRF` in the iSS files already hold the polarization $P^\mu = S^\mu/(1/2)$, not the spin vector $S^\mu$. Multiplying them by 2 doubles every result. This is a known labelling inconsistency in Chun's outputs (the physics itself is correct though! Only the output file has outdated labels).

How we know:

- The per-cell file is written by the `SpinPolarization` class of iSS (repository `chunshen1987/iSS`, branch `Spin_polarization`). Commit `a8e1b80` (2021-03-27, "fixed a missing factor of 2 for the averaged spin polarization") changed the normalization so that the arrays named `Smu` hold $P^\mu$, while the file header kept its `S^mu` labels. Before that commit, the same arrays held $S^\mu$.
- Since that commit, the code multiplies the three independent terms of the Levi-Civita contraction by $1/(2m)$, using the thermal vorticity from MUSIC (which already carries the $-1/2$ of its definition). That is exactly twice Becattini's spin vector, $S^\mu = -\frac{1}{8m}\epsilon^{\mu\rho\sigma\tau}p_\tau \int d\Sigma_\lambda p^\lambda n_F(1-n_F)\varpi_{\rho\sigma} / \int d\Sigma_\lambda p^\lambda n_F$, i.e. $P^\mu = S^\mu/(1/2)$ as in Eqs. (2)-(3) of the paper.
- The build that produced our files, `/home/vribeiro/vorticity/iSS-pol` in jarvis4, has the post-commit normalization:

```
186:    const double prefactor = 1./(2.*mass);
610:            (omega_yz*pmu[1] - omega_xz*pmu[2] + omega_xy*pmu[3])
```

Consequences:

- Both column blocks are used as they are. The column-to-polarization factor is 1, defined once in the producer and recorded in `Config/` of every output file.
- Vitor's plots (100 x column, labelled P) and the paper are consistent with this.
- The old `lambda_pol_toy_model.cxx` uses `Pol = 2S`, so its truth layer is twice too large (see Known issues and TODO).
- If the hydro outputs are ever regenerated, repeat the check on the iSS version used. `1./(2.*mass)` with no `0.5*(` in front of the vorticity terms means P. `-1./(8.*mass)` or a `0.5*(` means S.

```
grep -n "prefactor = \|omega_yz\*pmu\[1\]" <iSS>/src/spin_polarization.cpp
```

## Inputs

### Data layout

```
/storage1/vribeiro/<dataset>/<sample>/<i>/results/Smu_dpTdphidy_Thermal_rapidity_3122.dat
```

Vitor's events in jarvis4: fluctuating (lumpy) initial conditions for Pb-Pb collisions at 2.76 TeV, Lambda (PDG code 3122), thermal vorticity only, rapidity grid, events numbered from 0. The reference dataset is `lumpy_events_40_50` (40-50%, 250 events with iSS output in each sample, of the 300 folders of no_bullet). `hydroTruthLab_datasets.conf` lists every runnable dataset, with its samples, jet file and centrality, and `surveyHydroDatasets.sh` is how they were found:

```
./surveyHydroDatasets.sh                                       # every dataset under /storage1/vribeiro, a few lines each
./surveyHydroDatasets.sh -p break_symmetry lumpy_events_40_50/no_bullet   # pairing of a bullet sample with a no_bullet sample elsewhere
```

It reports, per sample, the events and those without iSS output, the grid, the MUSIC parameters that vary or stay constant, the format of the jet files, and whether bullet and no_bullet event i start from the same initial condition (comparing the tails of the initial-condition paths, since some datasets moved from `/hadrex1/storage2/vribeiro/vorticity_data` to `/storage1/vribeiro`).

Some events have no iSS output (an iSS run that never finished: 113 of `lumpy_events_0_5_fine`, half of `bullet_E30`, most of `bullet_P20_E*`). The producer skips them with a warning, records them in `Config/skippedEvents`, and PerConfig pairs the samples on the events both have.

### File format and grid

One header line plus 30000 lines:

```
# y  pT[GeV]  phi  dN/(pTdpTdphidy)[GeV^-2]  S^t  S^x  S^y  S^z  S^t_LRF  S^x_LRF  S^y_LRF  S^z_LRF
```

Columns 5-8 are the lab polarization $P^\mu$ and columns 9-12 the polarization in the Lambda rest frame $P^{*\mu}$ (despite the labels, see the warning above). Values are printed with 7 significant digits, which sets the floor of every closure check.

The grid has 10 y nodes (-0.45 to 0.45, step 0.1), 30 pT nodes (0.05 to 2.95 GeV/c, step 0.1) and 100 phi nodes ($2\pi k/100$). Every node is the centre of its cell (midpoint rule). The pT nodes are stock iSS, the y nodes come from a local edit in iSS-pol (stock iSS places them on [-5, 5] including the endpoints), and phi is periodic, so each phi node is the centre of its own cell. The producer derives the widths from the grid itself (no hardcoded 0.1 or $2\pi/100$), checks that the spacing is uniform and checks that every event has exactly the same (y, pT, phi) columns as event 0.

Each cell carries the weight

$$w = \frac{dN}{p_T\,dp_T\,d\phi\,dy}\;p_T\,\Delta p_T\,\Delta y\,\Delta\phi,$$

with $p_T$ the node value. The cell measure $p_T\,\Delta p_T\,\Delta y\,\Delta\phi$ is the same in every event, so it is computed once. Kinematic cuts are applied on cell centres and should sit on cell edges, as the defaults 0.5 < pT < 1.5 GeV/c and |y| < 0.5 do. The producer warns otherwise.

> The fourth column is $\int d\Sigma_\lambda p^\lambda f_0/(\hbar c)^3$ in this iSS version, without the degeneracy factor $g/(2\pi)^3$. The sum of $w$ is therefore only proportional to the Lambda yield and is never used as an absolute count. Every quantity here is a ratio, in which the constant cancels.

### The LRF block

iSS fills it with `transformSmuToLRF`, the pure boost along $\vec p$ of a four-vector orthogonal to $p$:

$$P^{*0} = \frac{E P^0 - \vec p\cdot\vec P}{m}, \qquad \vec P^{*} = \vec P - \frac{(\vec p\cdot\vec P)\,\vec p}{E(E+m)}.$$

We take $P^*$ directly from this block. Vitor never used it. Our own boost of the lab block (ROOT Math, as in the consumer) is kept **only as a closure check** (I trust Chun's code way more than my own!).

### Jet file

One per dataset (the registry names it), one line per event, for example `lumpy_events_40_50/random_parameters.dat`, columns `#event n_random phi(rad) momentum_x momentum_y`. The jet of a no_bullet event is that of the paired bullet event, so one file serves both samples. Rows have five fields, or four without `n_random`: the files of `lumpy_events_0_5_fine` (`phi_parameters.dat`, 180 directions on a regular grid in $[0, 2\pi)$) and `break_symmetry` keep the five-column header over rows of four, so the producer takes the layout from each row. It checks that `atan2(momentum_y, momentum_x)` agrees with `phi(rad)` modulo $2\pi$.

### MUSIC bullet settings

The bullet block of the MUSIC configuration of event 118, taken as representative:

```
total_current_energy 59.6
total_current_momentum_x -42.2447 #Gev/c
total_current_momentum_y 8.02405 #Gev/c
total_current_momentum_eta_s 0 #Gev/c
bullet_size_x 1 # fm
bullet_size_y 1 # fm
bullet_size_eta_s 0.4
bullet_position_x -0 # fm
bullet_position_y 0 # fm
bullet_position_eta_s 0 # fm
bullet_initial_time 1.001
bullet_final_time 1.002
```

The energy-momentum current is purely transverse (`total_current_momentum_eta_s 0`), so the trigger axis always has $t_z = 0$. It is deposited at the centre of the fireball, at $\tau \simeq 1$ fm/c. Other datasets change this block only: the deposited energy and momentum (`lumpy_events_60_70/bullet_E30`, `bullet_P20_E30`, `bullet_P20_E60`), or the position, $R$ fm along the jet direction for `position_lumpy_events_*/R+<R>` and against it for `R-<R>`, or fixed at $(2, 0)$ fm for `break_symmetry`. The comments of `hydroTruthLab_datasets.conf` give the values.

### KEY HYPOTHESIS: paired initial conditions

**Bullet event i and no_bullet event i evolve the same initial condition.** The MUSIC configurations of event 118 read, in both samples,

```
Initial_Distribution_input_filename /storage1/vribeiro/lumpy_bg_40_50/118/PbPb276_N_0p6553.root
```

and differ only by the bullet block above. Neither MUSIC nor the iSS polarization integral samples anything, so for a fixed initial condition the difference between the two samples is entirely due to the bullet.

`surveyHydroDatasets.sh` checks the same, event by event, and found it true for every pairing of `hydroTruthLab_datasets.conf`, including the no_bullet sample shared by several bullet samples (`lumpy_events_40_50/no_bullet` also serves `position_lumpy_events_40_50/*` and `break_symmetry`).

The whole `Difference/` folder of PerConfig rests on this hypothesis:

- Both samples must run with the same trigger mode and the same event range: the coordinator gives both the range of the bullet sample. Pairs are formed over the event IDs both samples have, so an event without iSS output in either sample enters neither curve of the difference. When only one sample misses it, the curves of the other sample keep it, and the drawn difference is not exactly the drawn bullet curve minus the drawn no_bullet curve. Every legend gives its event or pair count.
- Differences are taken pair by pair, and their errors are computed over pairs (see Errors).
- The trigger variants stay paired: random streams are seeded from (seed, eventId) instead of being drawn from one running sequence, both samples get the same seed (the coordinator sets one per run), and mixing partners are chosen by event index.
- For every shared `eventId`, PerConfig checks that `phiTrigger` agrees in the two event trees and stops otherwise.

If the hypothesis ever fails (another dataset, a different list of initial conditions), run PerConfig with `--paired=0`. It then takes differences of means and adds their errors in quadrature.

### Provenance

The analysis conventions mirrored here come from Vitor's codes under `/storage1/vribeiro/work_station/analysis/Polarization` in jarvis4: the grid nodes used as cell centres, the cell weights, the cuts 0.5 < pT < 1.5 GeV/c and |y| < 0.5, and the trigger axis $\hat t = (\cos\phi_J, \sin\phi_J, 0)$ with $\phi_J$ from the `phi(rad)` column of the jet file. No code was copied from there.

## Observables

All vectors are in lab coordinates (x along the impact parameter, z along the beam) unless stated otherwise. For a cell with momentum $\vec p$ and a trigger axis $\hat t$, the ring direction is

$$\hat n = \frac{\hat t\times\vec p}{|\hat t\times\vec p|}.$$

For $\hat t = (\cos\phi_t, \sin\phi_t, 0)$ and $\Delta\phi = \phi - \phi_t$, its lab components are

$$\hat n = \frac{\left(\sin\phi_t\,p_z,\; -\cos\phi_t\,p_z,\; p_T\sin\Delta\phi\right)}{\sqrt{p_z^2 + p_T^2\sin^2\Delta\phi}}.$$

The grid has no y = 0 node, so $p_z \neq 0$ and $\hat n$ is always defined.

Per-cell observables, named as in the histograms:

| Name | Definition | Remarks |
|---|---|---|
| `Ring` | $R = \vec P\cdot\hat n = \vec P^*\cdot\hat n$ | Full 3D ring observable, paper Eq. (4). Both forms agree because $\hat n$ is orthogonal to $\vec p$. |
| `RingZ` | $R_z = P^*_z\,n_z$ | Truth of the consumer's `useRingZ`. |
| `RingZT` | $R_z^T = (P^*_z - P^*_\parallel\hat p_z)\,n_z$ | Truth of `useRingZT`, with $P^*_\parallel = \vec P^*\cdot\hat p$. |
| `HelicityTerm` | $H = R_z - R_z^T = P^*_\parallel\,\hat p_z\,n_z$ | |
| `RingZLab` | $R_z^{lab} = P_z\,n_z = R_z^T + \gamma H$ | Projection of the lab polarization. |
| `RingPerp` | $R_\perp = R - R_z$ | |
| `NzSq` | $n_z^2$ | Denominator of $A_z$, precision input for $R_z$. |
| `NzSqSinSqTheta` | $n_z^2\sin^2\theta$ | Precision input for $R_z^T$, with $\theta$ the polar angle of $\vec p$. |

QA quantities: the decomposition of $\vec P^*$ along the ring basis, $P_R = \vec P^*\cdot\hat n$, $P_e = \vec P^*\cdot(\hat p\times\hat n)$ and $P_\parallel = \vec P^*\cdot\hat p$, the alignment $P_R/|\vec P^*|$ (which replaces the old toy's `NormPol`), and the consumer's angle chi, computed as in `computeProxyRing` of the consumer.

$A_z(\Delta\phi) = \langle R_z\rangle/\langle n_z^2\rangle$ per bin is the amplitude that the z projection would assign under the ansatz $\vec P^* = A\,\hat n$. The ansatz does not hold in this dataset ($R(\Delta\phi)$ looks like two summed Gaussians), so $A_z$ is a diagnostic, not a measurement. The theory of $R_z$, $R_z^T$ and the helicity term is in note 20 (LongitudinalRingProjection_AdditionalNotes) and is not repeated here.

**$R$ and $A_z$ need not agree, for any number of events.** Take a ring field $\vec P^* = P_R\,\hat n$ whose amplitude $P_R$ varies from cell to cell. Each cell then gives $R = P_R$ and $R_z = n_z^2 P_R$, so $\langle R\rangle$ averages $P_R$ with the yield weights $w$, while $A_z$ averages the same $P_R$ with the weights $w\,n_z^2$. ($A_z$ is also the least-squares fit of a constant ring amplitude to the beam-axis components of the decay directions, which is why note 20 calls it the natural amplitude of $R_z$.) Two weighted means of one quantity agree only if the quantity is uncorrelated with the extra weight:

$$A_z - \langle R\rangle = \frac{\mathrm{Cov}(n_z^2,\,P_R)}{\langle n_z^2\rangle}$$

(note 20, "The longitudinal and the full amplitude"). Here $n_z^2 = p_T^2\sin^2\Delta\phi/(p_z^2 + p_T^2\sin^2\Delta\phi)$ vanishes at $\Delta\phi = 0$ and $\pi$ and is largest near $\pm\pi/2$, while $P_R$ has its lobes wherever the jet put them. The covariance is a property of the polarization pattern, not a fluctuation that averages out with more events. A caricature: if only the Lambdas near $\Delta\phi = \pm\pi/2$ were polarized, with $P_R = P_0$, then $\langle R\rangle$ would be $P_0$ times their yield fraction, while $A_z \simeq P_0$, since $n_z^2$ already selects them. Beyond the ansatz, the non-ring components of hydro ($P_e$ and $P_\parallel$, nonzero here) enter $R_z$ but never $R$. So $R_z/R$ is not expected to approach $\langle n_z^2\rangle$, nor $A_z$ to approach $R$, and PerConfig reports both as diagnostics. $R_z/R$ is integrated only (per bin it divides by an $R$ that crosses zero), and is left undefined where $|R| < 10^{-12}$, i.e. rounding noise (the delta-scan averaged over $\delta$, for instance). $A_z(\Delta\phi)$ is drawn against $R(\Delta\phi)$: within one bin $n_z^2$ still varies with $p_T$ and $y$, so the two differ there too, by less.

### Averages

For an observable $X$, event $e$ and bin $b$, define $A_{e,b} = \sum_{c\in b} w_c X_c$ and $B_{e,b} = \sum_{c\in b} w_c$. Two averages are kept:

- **Pooled** (default name): $\bar X_b = \sum_e A_{e,b}/\sum_e B_{e,b}$, the average over all Lambdas of all events. Central value of every differential curve (Fig. 4 style).
- **Equal-weight event mean** (suffix `_EvtMean`): $\frac{1}{N}\sum_e A_{e,b}/B_{e,b}$. Central value of the per-event integrated quantities (Fig. 11 style).

Two integrated values are computed per event:

- **`D3N`** (d$^3N$-weighted): $\sum_b A_{e,b}/\sum_b B_{e,b}$, the multiplicity-weighted average of paper Eq. (4), and what an experiment measures.
- **`Flat`**: $\frac{1}{N_{bins}}\sum_b A_{e,b}/B_{e,b} = \frac{1}{2\pi}\int R(\Delta\phi)\,d\Delta\phi$, a flat average over $\Delta\phi$ of the differential ratio. This is what Vitor computed for Fig. 11.

They differ whenever the yield depends on $\Delta\phi$ within an event, which is exactly what the jet-induced flow causes.

### Binning

pT and y use the native nodes (30 and 10 bins). phi and $\Delta\phi$ use 100 bins centred on the phi lattice, spanning $[-\pi - \delta_{cell}/2,\ \pi - \delta_{cell}/2)$ with $\delta_{cell} = 2\pi/100$, and $\Delta\phi$ is wrapped into exactly that half-open range. The wrapped nodes of one (pT, y) row form a lattice whose spacing is the bin width, so every event puts exactly one node in every bin, whatever the trigger angle. For the `rp` trigger and for lattice-multiple deltas, the nodes sit on bin centres. Finer $\Delta\phi$ bins are impossible per event.

Vitor wrapped into $[-\pi, \pi]$ but binned over $(-\pi + \delta_{cell}/2,\ \pi + \delta_{cell}/2]$, so the cells in $[-\pi, -\pi + \delta_{cell}/2)$ fell into the underflow and the bin at $\Delta\phi \simeq \pi$ was empty in about half of the events (see Errors).

## Trigger modes

- **`jet`** (default for both samples): $\phi_t = \phi_J$ from the jet file. In no_bullet events, this is the direction the jet has in the paired bullet event. That makes the no_bullet sample a closure (no jet-induced ring can be present, so the integrated ring observable must be consistent with zero) and the subtraction partner of the paired differences.
- **`rp`**: $\phi_t = 0$, i.e. $\hat t = \hat x$, along the impact parameter. The reaction plane is fixed at $\Psi_{RP} = 0$ and the events are never rotated, so $\phi - \Psi_{RP} = \phi$. At 40-50%, this lines the trigger up with the elliptic-flow pattern, and the quadrupole of $P_z$ (Fig. 3 left) adds up coherently into $R(\Delta\phi)$. This is the reaction-plane picture of the smooth no-jet curve of Fig. 4, i.e. the behaviour one expects from events without a bullet. It is not a valid subtraction baseline for the bullet sample, whose $\phi_J$ is uniform relative to $\Psi_{RP}$.

The `jet` mode stops if an event has no jet entry, so the coordinator caps the range at the last event of the jet file. Mixing requires the `jet` mode, since all `rp` triggers coincide, and the coordinator turns it off for `rp` with a note in the log. `Paper/p2dPzLabVsPhiVsPt` (Fig. 3 left) does not depend on the trigger.

## Errors

The only randomness in this study is the event-to-event fluctuation of the initial conditions (and of $\phi_J$). Within an event everything is a deterministic integral, so errors computed from cell-to-cell spreads, as a TProfile filled per cell would give, depend on the grid and not on statistics. They are meaningless.

Every output is therefore filled at event level. A scratch profile collects the cells of one event with weight $w$, and each of its bins is then filled once into the master with value $A_{e,b}/B_{e,b}$:

- **Pooled master**, weight $B_{e,b}$: its mean is exactly $\sum_e A_{e,b}/\sum_e B_{e,b}$, and its error is the weighted spread over events divided by $\sqrt{N_{eff}}$, the standard approximation for a ratio of sums.
- **`_EvtMean` master**, weight 1: its mean and error are exactly the event mean and its standard error (SEM).

The jackknife is used for nonlinear quantities only, from the event tree, in PerConfig:

$$\sigma^2_{JK} = \frac{N-1}{N}\sum_i\left(\theta_{(i)} - \bar\theta_{(\cdot)}\right)^2,$$

with $\theta_{(i)}$ the estimate without event $i$ and $\bar\theta_{(\cdot)}$ the mean of the $\theta_{(i)}$. The value reported is the full-sample estimate, and the jackknife bias $(N-1)(\bar\theta_{(\cdot)} - \theta)$ is stored next to it as a check. For a plain mean of per-event values, the jackknife is identical to the SEM. A pooled ratio of sums is nonlinear, though, so in practice PerConfig jackknifes every pooled quantity: the curves per bin, their integrals, `Flat` of the pooled curve, $A_z$ (per bin and integrated), $R_z/R$, the lobe features, and the fold and asymmetry of `QA/`. The event means go through the same code, where the jackknife reduces to the SEM exactly. Each replica subtracts one event from the totals, so the whole jackknife costs one pass per event.

**Lobe features** (`Integrated/`), for Ring and RingZ, for each sample and for the difference. The lobes are the bins with $\Delta\phi < 0$ and with $\Delta\phi > 0$ (the bins at 0 and $\pi$ excluded) where the full-sample curve is positive. That set, and the bin of the maximum of each lobe, are chosen once from the full sample, and every replica reuses them, so the features are smooth functions of the data and the jackknife applies. The features are the centroid $\sum\Delta\phi\,R/\sum R$ of each lobe, the separation of the two centroids, the height at the fixed maximum bin, and the dip at $\Delta\phi = 0$. The tables list them from $\Delta\phi < 0$ to $\Delta\phi > 0$: [centroid, separation, centroid] and [height, dip, height]. The lobes of $R_z$ are also compared with those of $R$, replica by replica: the shifts $R_z - R$ of the two centroids and of the separation, and the height ratios $R_z/R$ of the two lobes. The features of the difference are those of the difference curve, not differences of the features of the two samples. The no_bullet features are stored for completeness: with smooth initial conditions no_bullet has lobes of its own (Fig. 4 of the paper), which fluctuating initial conditions largely wash out.

**No folding outside `QA/`.** $R(\Delta\phi) = R(-\Delta\phi)$ holds on average only ($\vec P$ and $\hat n$ are both axial), not event by event, and even less with a 3D smudged trigger. `QA/` holds the folded curves and the asymmetry $R(\Delta\phi) - R(-\Delta\phi)$, which should be consistent with zero. Every other curve keeps both signs of $\Delta\phi$.

**Differences** (`Difference/`, see the key hypothesis) are taken estimate by estimate, $\theta_{bullet} - \theta_{no\_bullet}$, with one exception: $R_z/R$ becomes the ratio of the differences, the fraction of the jet-induced ring that the z projection sees (the difference of the two ratios would divide by the no_bullet $R$, ~0 for the jet trigger). Their error is the leave-one-pair-out jackknife, where replica $k$ drops pair $k$ from both samples together. For event means it reduces to the SEM of $d_i = X_{bullet,i} - X_{no\_bullet,i}$. On the same pairs, PerConfig also computes the independent-samples jackknife (drop one bullet event with no_bullet whole, then one no_bullet event with bullet whole, and add the two variances) and stores the ratio of the two errors as the gain of the pairing. With `--paired=0` the differences use all events of each sample and carry only the independent-samples error.

No fits anywhere. The Fig. 11 point shows three estimators side by side: the event mean of `Flat` with its SEM (Vitor's estimator, with the plain mean instead of his Gaussian fit), and `D3N` as an event mean and pooled.

### Traps found in Vitor's error macros (I may want to recheck these, but they do seem wrong at first glance)

They sit in the same `Polarization` folder. We reuse none of them, and we do not want to fall into the same traps:

- `do_jackknife.C` divides the sum of squares by $N(N-1)$ instead of multiplying it by $(N-1)/N$. The result is SEM/(N-1), 249 times too small.
- It fills event n into bin n+1 but reads bin n back: the empty underflow enters every replica, the last event never enters, and replica 0 removes nothing.
- Its $\Delta\phi$ axis has 100 bins of width $1.01\,\delta_{cell}$ instead of the per-event axis: the first point gets zero error, points 2-49 get the spread of their left neighbour, and the point at $\Delta\phi \simeq 0$ gets the spread of two bins summed.
- The band (jackknife of the event mean) is drawn around the pooled curve, which is a different estimator.
- `do_jackknife2.C` never clears its replica buffer, so every replica reads replica 0 and the error becomes |mean|/sqrt(N-1).
- `draw_bar_sorted_x.C` adds per-bin errors in quadrature through `IntegralAndError`, although the bins of one event move together.
- The per-event $\Delta\phi$ histograms leave the bin at $\Delta\phi \simeq \pi$ empty in about half of the events, and `Divide` turns those 0/0 into zeros that enter event means, the jackknife and the integrals.

The bands of the published Fig. 4 are far wider than SEM/(N-1) could give, so the paper was probably made with an earlier version of these macros.

## Delta-scan, mixing and smudging

All three are trigger variants computed in the same pass as the true trigger. They fill a reduced set, with $\Delta\phi$ measured from the variant trigger: `Ring`, `RingZ` and `NzSq` vs $\Delta\phi$, the `D3N` and `Flat` integrals of `Ring` and `RingZ`, their per-event distributions, and the ring overlap with the true trigger (below).

- **Delta-scan** (`--deltaScanSteps=n`): the trigger rotated by $\delta_k = 2\pi k/n$, with $n$ even and a divisor of the number of phi nodes (default 100, one step per cell). This is the deterministic version of event mixing. It fills the integrals of `Ring` and `RingZ` (`D3N` and `Flat`) and of `NzSq` (`D3N`) vs $\delta$, and maps vs $(\Delta\phi', \delta)$. For $t_z = 0$, $\hat n(-\hat t) = -\hat n(\hat t)$ (the half-turn identity of note 13), so the integrals obey $R(\delta + \pi) = -R(\delta)$ exactly in every event, and their uniform average over $\delta$ vanishes for any polarization field.
- **Uniform mixing** (`--mixDepth=K`): event $i$ takes the $\phi_J$ of events $i+1, \dots, i+K$ (mod $N$, over the range of the run). It is a random sample of the delta-scan with $\delta_{eff} = \phi_J(\mathrm{partner}) - \phi_J(i)$ and is kept as a null test. A partner only lends its $\phi_J$, read from the jet file, so an event skipped for lack of iSS output still serves as one. The integrals are also filled vs $\delta_{eff}$, with the delta-scan binning, and PerConfig compares them with the delta-scan vs $\delta$ through pulls (below).
- **Smudging** (`--smudgeSigma`, `--smudgeMode`, `--smudgeDraws`): `cone` rotates the trigger by a half-normal angle around a uniformly random perpendicular axis (Rodrigues rotation, as in the consumer's `forceJetDirectionSmudge`). This gives a 3D trigger with $t_z \neq 0$, although the source is purely transverse. `phi` smears the azimuth only and keeps $t_z = 0$. Smudging also fills the smudge angle, the trigger pseudorapidity and the integrals vs smudge angle (the cone analogue of the delta-scan).

**Pulls of mixing against the delta-scan.** In $\delta_{eff}$ bin $k$, the mixed integral is a weighted mean over a subsample of the events the scan uses at $\delta_k$. With $v_e$ the integral of event $e$ at $\delta_k$, $W_e$ its scan weight (the in-cut yield, total $T$), $\bar v$ the pooled scan value, $s^2 = \sum_e W_e(v_e - \bar v)^2/T$ the spread over all events, and $c_e$ the summed weight of the mixed entries of event $e$ in the bin (total $C$),

$$\mathrm{pull}_k = \frac{\mathrm{mix}_k - \bar v}{s\sqrt{\sum_e\left(c_e/C - W_e/T\right)^2}}.$$

The sum is the exact variance factor of $\mathrm{mix} - \bar v$ for independent events with a common spread. For a random subsample it reduces to $1/n_{eff} - 1/N_{eff}$: the $-1/N_{eff}$ removes the events the two share, which a quadrature sum of the two errors counts twice (for $N = 250$ and $n \approx 25$, that sum makes the pulls 10% too small). The spread comes from all events, not from the ~25 in the bin, whose own spread would inflate the RMS through the tails of the per-event values. Events, not entries, are the units: two partners of one event that land in the same $\delta_{eff}$ bin carry nearly the same value ($K(K-1)/2$ partner pairs per event, so about one such event per bin for $K = 10$), and counting them as independent entries makes the RMS about 4% too large. The producer stores the $\delta_{eff}$ bin of every partner in the tree (`mixDeltaEffBin`), and PerConfig checks bin by bin that the weights of those entries add up to the producer's profile. Toy studies with these settings, Gaussian or heavy-tailed per-event values and spread weights, give an RMS of $0.99 \pm 0.01$. A clear excess means the mixed entries are not a random subsample of the events, a genuine failure of mixing.

### Ring overlap

For every variant trigger $\hat t'$, each in-cut cell also gets the overlap of its ring direction $\hat n' = \hat t'\times\vec p/|\hat t'\times\vec p|$ with the true one $\hat n$, as `computeRingOverlap` in the consumer (`EventMixingQA/RingOverlap/`). The angle $\zeta$ turns $\hat n$ into $\hat n'$, right-handed about $\hat p$:

$$\cos\zeta = \hat n\cdot\hat n', \qquad \sin\zeta = (\hat n\times\hat n')\cdot\hat p.$$

Both ring directions are orthogonal to $\vec p$, so $\hat n' = \cos\zeta\,\hat n + \sin\zeta\,\hat e$ with $\hat e = \hat p\times\hat n$, and the ring that the variant measures is, cell by cell and exactly,

$$R' = \vec P^*\cdot\hat n' = P_R\cos\zeta + P_e\sin\zeta,$$

with $P_R$ and $P_e$ of the true trigger (`QA/Decomposition/`). The leakage ansatz of the consumer, $\vec P^* = A\,\hat n + B\,\hat e$ with constant $A$ and $B$, predicts $\langle R'\rangle = A\langle\cos\zeta\rangle + B\langle\sin\zeta\rangle$. Hydro tests it directly: with $A = \langle P_R\rangle$ (the true $R$) and $B = \langle P_e\rangle$ (`polED3N` in the tree), the prediction misses exactly $\mathrm{Cov}(P_R, \cos\zeta) + \mathrm{Cov}(P_e, \sin\zeta)$.

What to expect. For the delta-scan, $\langle\cos\zeta\rangle(\delta)$ is set by the geometry alone and changes sign under $\delta \to \delta + \pi$ (the half-turn identity again), so both sides of the ansatz vanish on the uniform average and the test is informative at fixed $\delta$. Uniform mixing samples that: $\langle\cos\zeta\rangle \simeq 0$, with a broad distribution. Smudging keeps $\cos\zeta \simeq 1 - O(\sigma^2)$.

The producer stores the per-bin sums of $w\cos\zeta$ and $w\sin\zeta$ in the tree (mixing and smudging, binned in the variant's $\Delta\phi'$) and their integrals per $\delta$ (delta-scan), so PerConfig can jackknife and pair them. It also fills the pooled curves and integrals, and the d$^3N$-weighted distributions `hCosZeta` and `hSinZeta` (mixing and smudging), and `h2dCosZetaVsDelta` and `h2dSinZetaVsDelta` (delta-scan).

PerConfig (`RingOverlap/`) jackknifes $\langle\cos\zeta\rangle$ and $\langle\sin\zeta\rangle$ per bin and integrated, and tests the ansatz: the prediction $L = R_{true}\langle\cos\zeta\rangle + \langle P_e\rangle\langle\sin\zeta\rangle$ against the measured $R'$, and the residual $R' - L$, which is exactly the covariance term above. Replica $k$ of the true trigger and of the variant drop the same event, so $L$ and $R' - L$ are jackknifed jointly, and for the difference they are formed pair by pair (`--paired=1` only). The test runs vs $\delta$ for the delta-scan and on the integrals for mixing and smudging: with constant $A$ and $B$ the ansatz cannot follow $R(\Delta\phi')$ bin by bin. At $\delta = 0$ the scan trigger is the true one, so $\langle\cos\zeta\rangle = 1$ and $\langle\sin\zeta\rangle = 0$ there, event by event (a closure of `DeltaScan/`).

## Polarization maps

`PolarizationMaps/` mirrors `RawPolarization/PolMaps/` of the consumer for the true trigger: the components of $\vec P^*$ over momentum planes, the d$^3N$ yields of the same planes and, on the PrimeJet planes, the ring scalars $R$ and $R_z$ (the consumer's `RingMaps`). Full grid, no cuts. Object names follow the consumer, with `h2dCountsVs*` written `h2dYieldVs*` (sums of weights, not counts). Three coordinate systems:

- **Lab/**: $(p_x, p_y)$, $(p_z, p_x)$ and $(p_y, p_z)$, with $P^*_x$, $P^*_y$ and $P^*_z$.
- **PrimeJet/**: the consumer's jet coordinate system, $\hat z' = \hat t$, $\hat x' = (\hat z - t_z\hat t)/t_T$ (the beam projected perpendicular to the jet) and $\hat y' = \hat t\times\hat z/t_T$. For our transverse trigger $\hat x' = \hat z$. The $x'$-$y'$ plane is transverse to the jet, and a ring shows up there as arrows circling the origin.
- **PrimeV0/**: the production plane $(p_z, p_T)$, with $\hat x'_{V0} = \hat p_T$ and $\hat y'_{V0} = \hat z\times\hat x'_{V0}$, the azimuthal direction. $P^*_{y'V0}$ is the ring observable with the beam as the trigger.

There is no `Aee/`: the azimuthal efficiency effect needs the decay proton, and hydro has no decays.

**Smoothing.** The grid is polar in $(y, p_T, \phi)$ and the maps are Cartesian. At $p_T = 3$ GeV/c, neighbouring nodes are 0.19 GeV/c apart in azimuth and 0.32 GeV/c apart in $p_z$, wider than the 0.1 GeV/c bins, so a map filled at the nodes has empty bins (244 inside the populated region of the $(p_x, p_y)$ map). With `--mapSubdivisions=n` (default 3), each cell is spread over the centres of its $n\times n\times n$ sub-cells in $(y, p_T, \phi)$. Each point carries the field of the node, as the midpoint rule of every other output already assumes, and the d$^3N$ share $(p_{T,\mathrm{point}}/p_{T,\mathrm{node}})/n^3$, which keeps the $p_T\,dp_T$ Jacobian and adds up to the cell weight exactly. With $n = 3$ no bin inside the populated region is empty. The titles say "smoothed map" whenever $n > 1$. Only `PolarizationMaps/` is smoothed: every other output uses the nodes.

> This is essentially just a visualization trick to make sure we don't get a bunch of empty bins in the middle of the pxpypz grids of `PolMaps/`, to get an easier comparison between Hydro and my O2Physics consumer output!

**Axes** come from the grid: bins of the $p_T$ cell width, transverse components up to the last $p_T$ edge ($\pm 3.0$ GeV/c), $|p_z|$ up to $m_T\sinh|y|$ at the outer $p_T$ and $y$ edges ($\pm 1.7$ GeV/c), and $p_T$ on its native nodes. They hold the whole momentum volume of the grid with less than one bin of margin.

**Vector fields** (PerConfig `PolarizationMaps/`). For each coordinate system and each sample, one canvas per set of planes, as in the consumer's per-config plots: the out-of-plane component of $\vec P^*$ as colour, and arrows of the in-plane components averaged over tiles of $4\times4$ bins (0.4 GeV/c), with the true angle on screen whatever the aspect of the panel. The drawing is ported from `DrawVectorFieldPanel()` of `auxiliaryPerConfigPlots.cxx`, with three changes for hydro. The maps are TH2Ds in percent, since a difference of two samples is no TProfile2D. A bin enters where its yield is positive (hydro has no counting statistics to cut on), and empty bins are left blank. The arrow of a tile is the yield-weighted mean of its bins, i.e. the d$^3N$-weighted mean polarization of the Lambdas in the tile. The scale box gives the 95th percentile of the tile magnitudes, to which the arrows are scaled. Its error propagates the pooled-master errors of the bins as if they were independent, as in the original, so it is indicative only. The bullet - no_bullet canvases show the difference of the two pooled maps, with errors in quadrature and the bullet yields for the populated bins and the tile weights. PrimeJet also gets the $R$ and $R_z$ maps, and every system the yields.

`PolVsAngles/` (formerly `PolMaps/`) keeps the angular dependences: the lab $\vec P$ and $\vec P^*$ components vs $\phi - \Psi_{RP}$ and $(\phi - \Psi_{RP}, p_T)$, and the components in trigger-aligned coordinates ($\hat x'$ along $\hat t$, $\hat y' = \hat z\times\hat x'$, $\hat z$) vs $\Delta\phi$ and $(\Delta\phi, p_T)$. Its former momentum-space ring maps moved to `PolarizationMaps/PrimeJet/`, in the consumer's convention.

## Expected precision

For an ideal detector and a small polarization, each Lambda contributes the estimate $\frac{3}{\alpha}\,\hat k^*\cdot\vec e$, with $\hat k^*$ the proton direction in the Lambda rest frame. Its variance is $3|\vec e|^2/\alpha^2$, with $\vec e = \hat n$ for $R$, $n_z\hat z$ for $R_z$ and $n_z(\hat z - \hat p_z\hat p)$ for $R_z^T$. Hence

$$\sigma_b = \frac{\sqrt{3\,\langle|\vec e|^2\rangle_b}}{|\alpha|\sqrt{N_b}}, \qquad |\vec e|^2 = 1,\ n_z^2,\ n_z^2\sin^2\theta,$$

with $N_b = N f_b$ and $f_b$ the hydro yield fraction of bin $b$. Hydro only provides the shapes $\langle|\vec e|^2\rangle_b$ and $f_b$, both taken from the bullet sample. PerConfig always fills `Precision/` with $\sigma_b\sqrt N$ per bin and integrated ($f = 1$), with $\alpha$ from `--alphaLambda`. With `--nLambdaPrecision=N` it also stores $\sigma_b$ and the expected significance of the jet-induced signal: the paired difference over $\sigma_b$, the error of one measurement with $N$ Lambdas (a baseline measured with the same statistics would add a factor $\sqrt2$). $R_z^T$ has no per-bin difference (no per-bin sums in the tree), so its significance is integrated only. The coordinator passes $N = 288574517$, the Lambdas left after every cut of the O-O analysis: the significance of a signal as clean as hydro allows, a ceiling rather than a prediction (Pb-Pb hydro, ideal detector).

> A small comment on that value: totally preliminary! Also, that counts $\Lambda$'s in a larger rapidity range than what is available in Hydro, so we should probably rescale this by some factor of around 3 ($|y|_{data}<1.5$, whilst $|y|_{hydro}<0.5$ after the cuts). It is only a rough estimate for reference though. There isn't much problem with this number being slightly off -- it offers only a better estimate of the required statistics (or of the uncertainty bars) in this specific O-O scenario.

## Closure checks

| Check | Where | Expected |
|---|---|---|
| Uniform grid, identical (y, pT, phi) columns in every event | producer, every event | hard stop otherwise |
| Cuts on cell edges | producer, start | warning otherwise |
| `atan2(momentum_y, momentum_x)` vs `phi(rad)`, modulo $2\pi$ | producer, jet file | below 1e-4 rad, hard stop otherwise |
| $P\cdot p/m$ (Minkowski) | `QA/Closure/hPdotP` | print precision. Tests the column order and the kinematics. |
| $P^{*0}$ | `QA/Closure/hPStarTime` | print precision |
| Our boost of the lab block vs the LRF block | `QA/Closure/hBoostClosure` | print precision |
| $\vec P\cdot\hat n - \vec P^*\cdot\hat n$ | `QA/Closure/hRingLabMinusStar` | print precision |
| Delta-scan at $\delta = 0$ vs the true trigger, event by event (`Ring`, `RingZ`, `NzSq`, and $\langle\cos\zeta\rangle - 1$, $\langle\sin\zeta\rangle$) | PerConfig `DeltaScan/hClosure` | float noise |
| $\max_k \lvert R(\delta_k + \pi) + R(\delta_k)\rvert$ per event | producer `DeltaScan/hClosureHalfTurn`, summed over the samples in PerConfig `DeltaScan/` | float noise |
| Half-turn and uniform average over $\delta$ of the pooled $R(\delta)$ | PerConfig `DeltaScan/hClosure` | float noise |
| Mixing vs $\delta_{eff}$ against the delta-scan vs $\delta$ | PerConfig `Mixing/hPullMean`, `Mixing/hPullRms`, `Config/mixingPulls` | pull mean ~0, RMS ~1 (statistical scatter about 0.07 for 100 bins) |
| Mixed entries of the tree (`mixDeltaEffBin`) vs the producer's profiles | PerConfig, per $\delta_{eff}$ bin | equal weights, warning and no pulls otherwise |
| no_bullet in `jet` mode: integrated `Ring` and `RingZ` | PerConfig log, `Config/noBulletClosure`, `Integrated/` | consistent with zero |
| Pooled-master error vs jackknife, per bin | PerConfig `Jackknife/` | close agreement |
| `_EvtMean` profile error vs jackknife | PerConfig `Jackknife/hErrorRatioIntegrated` | $\sqrt{(N-1)/N}$ exactly (TProfile uses the population spread) |
| $R(\Delta\phi) - R(-\Delta\phi)$ | PerConfig `QA/` | consistent with zero |
| Vitor's pooled bullet curve (`data/bullet_lumpy_40_50.root` in his `Polarization` folder, `h_RP_delta_phi`, in percent) vs `Paper/pRingVsDeltaPhi` | by eye | equal away from $\Delta\phi = \pm\pi$ (same columns, weights and trigger) |
| Paper Fig. 2 (right) and Fig. 3 (left) | by eye | same patterns |

## How to run me

### Coordinator

```
./runHydroTruthLab.sh [OPTIONS]
```

| Option | Effect |
|---|---|
| `-h`, `--help` | Show the help message. |
| `-p`, `--post-process-only` | Skip the producers and rerun PerConfig on the existing producer files. |
| `-t`, `--tag NAME` | Name of the campaign folder (default: `default`). |
| `-o`, `--only RUN` | Run only this run, e.g. `lumpy4050_default`. Repeatable. |
| `-j`, `--jobs N` | Runs at a time (default: `MAX_PARALLEL_RUNS`, half the cores). 1 runs them one after the other. |
| `--runs FILE`, `--datasets FILE` | Other registries (default: the two next to the script). |
| `--cleanup` | Remove the compiled executables and exit. |

**Configuration.** The CONFIGURATION block at the top of the script holds the defaults of every configurable, the registries only what changes:

- `PRODUCER_OPTS` and `PERCONFIG_OPTS`: the options of each macro, each commented with its source (the consumer for the Lambda mass and the decay parameter, the PRC paper for the cuts, the O-O analysis for `--nLambdaPrecision`);
- `HYDRO_BASE` and `INPUT_NAME`, `FIRST_EVENT` and `MAX_EVENTS` (0 = no cap), `MAX_PARALLEL_RUNS`, and the output base folder `OUT_BASE` (set only there).

`hydroTruthLab_datasets.conf` has one line per dataset, `label centrality bullet noBullet jetFile`, with paths relative to `HYDRO_BASE` and `none` for a dataset without no_bullet (PerConfig then runs on the bullet sample alone). `hydroTruthLab_runs.conf` has one line per run, `run <dataset> <variation> [option=value ...]`. An override replaces the default of that option in whichever macro has it, so a run lists only what it changes:

```
run lumpy4050 default
run lumpy4050 rp          trigger=rp
run lumpy4050 smudge0.1   smudgeSigma=0.1 seed=smudgeScan4050
```

An option missing from the two blocks is added there first, with its value, so that no run depends on a default hidden in a macro. `--jetFile`, `--inputName`, `--firstEvent`, `--nEvents` and `--seed` of the producer, and `--centralityLabel` and `--runLabel` of PerConfig, come from the registries. `trigger=rp` also turns mixing off, with a note in the log, since all `rp` triggers coincide. The registries are checked before anything runs: an unknown dataset, a duplicate run, an override that names no option of either macro (a typo), a seed above 32 bits or a comparison that names an unknown run stops the script there.

**Event range.** Both samples of a run get the same range, which the mixing partners and the pairing need: from `FIRST_EVENT` to the last event of the bullet sample with iSS output, capped at the last event of the jet file for the jet trigger or mixing, and at `MAX_EVENTS`. Events of the range without iSS output in either sample are listed in the log ("which the producer skips"), and the producer skips them with a warning.

**Seeds.** A run draws its own seed from the clock, shared by its two samples (paired smudging). `seed=<integer>` fixes it, to repeat a run exactly. `seed=<name>` gives every run of the invocation with that name one clock-drawn seed: the smudge draws are then identical, each rescaled by its width, so a width scan differs only through the width (common random numbers), and the first draws of a larger `smudgeDraws` repeat those of a smaller one. A group shares its seed within one invocation only, and rerunning part of it with `--only` gives that part a new seed (the log says so). `cone` and `phi` consume their draws differently, so they never share them. The seed of every run is printed in the plan and stored in each producer's `Config/`.

**Execution.** The coordinator compiles both macros ahead of time, plans every run (range, seed, skipped events, all printed before anything runs), then starts up to `--jobs` runs at once. Each run starts its two producers together, single-threaded, then PerConfig. Progress lines carry the run label, and a table of failed or skipped steps closes the log. Every output is written as `<file>.part` and renamed only on success, after removing the file of an earlier invocation, so a failed or interrupted run never leaves a partial or stale file under the final name. Ctrl+C stops every run, with its producers. The compilation uses the flags of the main analysis chain, without `-DNDEBUG`, so `assert` stays active:

```
g++ -Wall -Wextra -Wpedantic -Wshadow -O3 -march=native -fno-math-errno -flto=auto -fno-trapping-math $(root-config --cflags) -o hydroTruthLab.exe hydroTruthLab.cxx $(root-config --glibs) -lGenVector
```

**Output layout**, one flat folder per run, every file named after it:

```
${OUT_BASE}/${RUN_TAG}/
  logHydroTruthLab.log                                      appended, one block per invocation
  backupOf_runHydroTruthLab.sh                              the script and the registries as last run: the complete configuration
  backupOf_hydroTruthLab_runs.conf
  backupOf_hydroTruthLab_datasets.conf
  lumpy4050_default/
    hydroTruthLab_lumpy4050_default_bullet.root     (.log)
    hydroTruthLab_lumpy4050_default_noBullet.root   (.log)
    hydroTruthLabPerConfig_lumpy4050_default.root   (.log)
  lumpy4050_smudge0.1/
    ...
```

Downstream code never parses these names: everything that identifies a file (run label, centrality, sample folders, options, seed, skipped events) is in its `Config/`.

**Comparisons.** `hydroTruthLab_runs.conf` also lists the comparisons for `hydroTruthLabSummaryPlots.cxx`, `compare <name> <x> <run> <run> ...`, with `x` an option of either macro (its value in each run), `centralityLabel`, or `runLabel` (one point per run, in the order listed). The coordinator checks them and runs nothing for them.

### Running the executables by hand

```
./hydroTruthLab.exe <sampleDir> <outputFile> [--option=value ...]
```

| Option | Default | Meaning |
|---|---|---|
| `--jetFile` | none | Jet file, rows of 5 or 4 fields (see Jet file). Required for `--trigger=jet` and for mixing. |
| `--trigger` | `jet` | `jet` or `rp` (see Trigger modes). |
| `--firstEvent`, `--nEvents` | 0, 250 | Event range. Events without the input file are skipped with a warning. |
| `--inputName` | `Smu_dpTdphidy_Thermal_rapidity_3122.dat` | File name inside `<sampleDir>/<i>/results/`. |
| `--ptMin`, `--ptMax`, `--absYMax` | 0.5, 1.5, 0.5 | Kinematic cuts on cell centres (GeV/c for pT). |
| `--lambdaMass` | 1.115683 | GeV/c^2, used to build $E$ and $p_z$ from (y, pT). |
| `--deltaScanSteps` | 100 | Number of deltas, 0 = off. Even, and a divisor of the number of phi nodes. |
| `--mixDepth` | 0 | Mixing partners per event, 0 = off. |
| `--smudgeSigma` | 0 | Smudge width in rad, 0 = off. The consumer uses 0.05 jetR, i.e. 0.02 for jetR = 0.4. |
| `--smudgeMode` | `cone` | `cone` or `phi`. |
| `--smudgeDraws` | 1 | Smudged triggers per event. |
| `--seed` | current time | Base seed. Each event uses its own stream, seeded from (seed, eventId). Paired smudging needs the same seed for both samples: the coordinator passes one per run. |
| `--mapSubdivisions` | 3 | `PolarizationMaps/`: each cell spread over $n^3$ points of its momentum volume. 1 = nodes only (see Polarization maps). |

```
./hydroTruthLabPerConfig.exe <bulletFile> <noBulletFile|none> <outputFile> [--option=value ...]
```

| Option | Default | Meaning |
|---|---|---|
| `--paired` | 1 | Paired differences (key hypothesis). 0 = differences of means. |
| `--alphaLambda` | 0.749 | Decay parameter, used only in `Precision/`. |
| `--nLambdaPrecision` | 0 | Number of Lambdas for $\sigma$ and the significance in `Precision/`. 0 = $\sigma\sqrt N$ only. |
| `--centralityLabel` | `40-50%` | Label of the Fig. 11 point. |
| `--runLabel` | `none` | Run of the coordinator, `<dataset>_<variation>`, stored in `Config/` for the cross-run plots. |

With `none` as the no_bullet file, PerConfig skips everything that needs it: the difference curves, `Difference/`, Fig. 3 (left) and the significances.

## Output reference

### Producer: `hydroTruthLab_<run>_<sample>.root`

```
Config/            one TNamed per option (name = option, title = value), input paths, grid summary, skippedEvents (ranges, or none), column-to-polarization
                   factor and its provenance note
Inputs/            yield vs phi, pT, y and (dphi, pT), trigger angles, per-event yield in the cuts
Paper/
  p2dPzLabVsDeltaPhiVsPt                    Fig. 2 (right): lab P_z vs (dphi, pT), |y| < 0.5, all pT
  p2dPzLabVsPhiVsPt                         Fig. 3 (left): lab P_z vs (phi - Psi_RP, pT), Psi_RP = 0
  pRingVsDeltaPhi                           Fig. 4 analogue at the centrality of the run (+ _EvtMean)
  pRingFlat, pRingD3N                       Fig. 11 point, single-bin profiles (+ _EvtMean)
  hRingFlatPerEvent, hRingD3NPerEvent       per-event distributions (Fig. 11 band)
Ring/                         true trigger, inside the cuts
  p<X>VsDeltaPhi              X = Ring, RingZ, RingZT, HelicityTerm, RingZLab, RingPerp, NzSq, NzSqSinSqTheta
  p<X>D3N                     integrated, same X
  pRingVsDeltaPhi_EvtMean, pRingZVsDeltaPhi_EvtMean
  pRingFlat, pRingZFlat       (+ _EvtMean)
  p2dRingVsDeltaPhiVsPt, p2dRingVsDeltaPhiVsY, pRingVsPt, pRingVsY      same for RingZ
  hRingZFlatPerEvent, hRingZD3NPerEvent
  NoCuts/                     pRingVsDeltaPhi, pRingZVsDeltaPhi on the full grid (old toy "TrueValue")
PolVsAngles/       lab P and P* components vs phi and (phi, pT), components in trigger-aligned coordinates (x' along t) vs dphi and (dphi, pT)
PolarizationMaps/  full grid, smoothed (see Polarization maps)
  Lab/             p2dP{x,y,z}Star_vs{PxPy,PzPx,PyPz}, h2dYieldVs{PxPy,PzPx,PyPz}
  PrimeJet/        p2dP{x,y,z}StarPrimeJet_vs{PxPy,PzPx,PyPz}PrimeJet, h2dYieldVs{PxPy,PzPx,PyPz}PrimeJet, p2dRingVs{...}PrimeJet, p2dRingZVs{...}PrimeJet
  PrimeV0/         p2dPxStarPrimeV0_vsPzPt, p2dPyStarPrimeV0_vsPzPt, p2dPzStar_vsPzPt, h2dYieldVsPzPt
QA/
  Decomposition/   P_e, P_par and P_R/|P*| vs dphi (P_R is Ring/pRingVsDeltaPhi itself)
  Chi/             Ring and RingZ vs the consumer's chi, chi distribution
  Closure/         hPdotP, hPStarTime, hBoostClosure, hRingLabMinusStar
DeltaScan/         Ring, RingZ (D3N, Flat) and NzSq (D3N) vs delta (pooled + _EvtMean), p2dRingVsDeltaPhiVsDelta, p2dRingZVsDeltaPhiVsDelta, hClosureHalfTurn,
                   ring overlap: pCosZetaD3NVsDelta, pSinZetaD3NVsDelta, h2dCosZetaVsDelta, h2dSinZetaVsDelta
Mixing/            reduced set (below), Ring and RingZ integrals vs delta_eff
Smudge/            reduced set (below), hSmudgeAngle, hTriggerEta, Ring and RingZ integrals vs smudge angle
                   reduced set: pRingVsDeltaPhi, pRingZVsDeltaPhi, pNzSqVsDeltaPhi, pRingD3N, pRingFlat, pRingZD3N, pRingZFlat (+ _EvtMean), their per-event
                   distributions, and the ring overlap: pCosZetaVsDeltaPhi, pSinZetaVsDeltaPhi, pCosZetaD3N, pSinZetaD3N, hCosZeta, hSinZeta
events             TTree, one entry per event (below)
```

### Event tree `events`

```
eventId, phiJ, phiTrigger, nCellsInCuts, sumWInCuts
<x>D3N                        x = ring, ringZ, ringZT, helicityTerm, ringZLab, ringPerp, nzSq, nzSqSinSqTheta
ringFlat, ringZFlat
polED3N                                                                   P_e of the true trigger (ring overlap)
sumWDeltaPhi[100], sumARingDeltaPhi[100], sumARingZDeltaPhi[100], sumANzSqDeltaPhi[100]
nDelta, scanRingD3N[nDelta], scanRingFlat[nDelta], scanRingZD3N[nDelta], scanRingZFlat[nDelta], scanNzSqD3N[nDelta], scanCosZetaD3N[nDelta], scanSinZetaD3N[nDelta]
mixRingD3N, mixRingFlat, mixRingZD3N, mixRingZFlat                        pooled over the K partners
nMixPartners, mixDeltaEffBin[nMixPartners]                               K, and the delta_eff bin of each partner (from 0, -1 outside the axis)
nMixDeltaPhiBins, mixSumWDeltaPhi[100], mixSumARingDeltaPhi[100], mixSumARingZDeltaPhi[100], mixSumANzSqDeltaPhi[100],
                  mixSumACosZetaDeltaPhi[100], mixSumASinZetaDeltaPhi[100]                              the same, per bin (no bins when mixing is off)
smudgeRingD3N, smudgeRingFlat, smudgeRingZD3N, smudgeRingZFlat            pooled over the draws
nSmudgeDeltaPhiBins, smudgeSumWDeltaPhi[100], smudgeSumARingDeltaPhi[100], smudgeSumARingZDeltaPhi[100], smudgeSumANzSqDeltaPhi[100],
                     smudgeSumACosZetaDeltaPhi[100], smudgeSumASinZetaDeltaPhi[100]                     the same, per bin (no bins when smudging is off)
```

The per-bin sums let PerConfig jackknife and pair any pooled quantity, the true, mixed and smudged curves included. The delta-scan stores integrals only, which PerConfig turns back into sums with `sumWInCuts`: exact for transverse triggers, which see every in-cut cell. `mixDeltaEffBin` places each event's partners in the $\delta_{eff}$ bins of `Mixing/` for the pulls, so PerConfig never needs the partner rule.

### PerConfig: `hydroTruthLabPerConfig_<run>.root`

Every folder holds its canvases at the top, for a first skim, and every object they draw in its `Parts/` subfolder. "Samples" below means bullet, no_bullet and the bullet - no_bullet difference.

```
Config/            options (runLabel included), inputs, event and pair counts, the no_bullet and delta-scan closures, the mixing pulls, the leakage test,
                   every warning. ProducerBullet/ and ProducerNoBullet/: copies of the Config/ of each producer file (seed, event range, skipped events, ...)
Paper/             in percent: Fig. 2 (right), Fig. 3 (left), Fig. 4 analogue (R(dphi) of the samples, jackknife bands), Fig. 11 point (three estimators)
Integrated/        per sample: every integral of the true trigger (in units of P, grouped by observable, and dimensionless), lobe positions and heights
                   of Ring and RingZ, R against R_z (lobe positions overlaid, centroid shifts, height ratios)
DifferentialRing/  per sample: R, R_z, R_z^T, H and R_z^lab vs dphi (true trigger), <n_z^2>(dphi), and cAmplitudeZ: R and A_z of both samples,
                   their differences on a lower pad
Difference/        paired against independent-samples errors, per bin and for the integrals, for every trigger set, and the gain of the pairing
DeltaScan/         per sample: R(delta) and R_z(delta). Closures (delta = 0, half-turn, uniform average) and the (dphi', delta) maps
Mixing/            per sample: mixed R(dphi) and R_z(dphi), true against mixed curves and integrals, mixing vs delta_eff against the scan, pulls
Smudge/            per sample: smudged R(dphi) and R_z(dphi), true against smudged curves and integrals, smudged/true integrals, integrals vs smudge angle
RingOverlap/       <cos zeta> and <sin zeta> of every variant (per bin and integrated), their distributions, and the leakage ansatz (README, "Ring overlap")
PolarizationMaps/  Lab/, PrimeJet/ and PrimeV0/: vector fields of P* per sample, in percent, R and R_z on the PrimeJet planes, yields (README, "Polarization maps")
Precision/         sigma x sqrt(N) per dphi bin and integrated for Ring, RingZ and RingZT, plus sigma and the significance for --nLambdaPrecision > 0
Jackknife/         validation: pooled-master over jackknife errors (per bin for every set, and for the single-bin profiles), jackknife bias over error
QA/                folded curves and asymmetries R(dphi) - R(-dphi) of Ring and RingZ, the only place where curves are folded
```

Drawing conventions:

- Every curve with errors is a band with its central line, all bands drawn before all lines. The bands are TGraphErrors drawn with option "3" in opaque light colours: JSROOT (the vscode viewer) does not render the E3 bands of a TH1, and transparent colours are session-dependent indices that come out wrong when a canvas is reopened. Where bands overlap, the one drawn later hides the others, but every central line stays visible.
- `Difference/` nests the bands of one central line: the wide independent-samples band behind, the paired band in front.
- Curves without errors (pulls, error ratios, gains, precision) are lines, tables of labelled values are markers side by side, distributions are histograms. Tables of integrals are grouped by observable, with dashed separators: $R$ by its four estimators, $R_z$ by the same four, then the other projections.
- Labels: `d^{3}N weight` (pooled, d$^3N$-weighted), `flat` (flat average over $\Delta\phi$) and `evt. mean` (equal-weight event mean). The log uses the ASCII forms (`R_d3N`, `R_flat (evt. mean)`, ...).

Raw units everywhere but `Paper/` and `PolarizationMaps/`. Every curve with per-bin sums in the tree carries jackknife errors. The few without them (R_z^T, H and R_z^lab vs dphi, the (dphi', delta) maps, mixing vs delta_eff, the integrals vs smudge angle, the polarization maps) are the producer's pooled curves with their pooled-master errors, and their titles or legends say so. An undefined estimate (for instance $A_z$ at $\Delta\phi = 0$ for the rp trigger, where $n_z = 0$) is left as an empty bin in `Parts/`, out of every drawn band or line, and logged.

## Known issues and TODO

- **Old toy model.** `lambda_pol_toy_model.cxx` multiplies the iSS columns by 2 (`Pol = 2S`), so its truth-level polarization is twice the true value. Not corrected yet.
- **Paper scope.** In `lumpy_events_40_50` the jet is deposited at the centre. Fig. 2 (right) and Fig. 3 (left) are reproduced directly. Fig. 4 is published for 0-5% only, so we build its analogue at the centrality of every run (`lumpy_events_0_5_fine` included). Fig. 11 gets one point per run, and the centrality comparison of the runs registry collects them. Whether the `position_lumpy_events_*` datasets are those of Fig. 10 (jets at random positions) is not checked yet.
- **`hydroTruthLabSummaryPlots.cxx`.** The cross-run plots of the `compare` lines (Fig. 11 vs centrality, integrals vs smudge width, deposition, injection position, seed spread) are not written yet. It will read the `Parts/` and `Config/` of the PerConfig files.
- **Rotated dataset.** `lumpy_events_40_50_rotated` (iSS output in events 0-99 of both samples) is in the datasets registry, commented out: what "rotated" means is not known yet.
- **Config/ and hadd.** The histograms and the tree merge correctly, but `Config/` strings cannot be merged, so a merged file would not carry a correct event range. Not needed while each sample runs in one process.
- **$A_z$.** Diagnostic only, since the ansatz $\vec P^*(\phi_\Lambda-\phi_{Jet}) = A\,\hat n$ (with A constant in $\phi_\Lambda-\phi_{Jet}$) does not hold here. It was something I tested by hand, analytically solvable, based on Figure 7.13 from Vitor's dissertation. See the logs close to 775 or earlier -- or just look into internal note 26, as it has most of the "Event Mixing Validity" considerations.
- **Expected precision.** Ideal detector: no efficiency, no acceptance holes and no fake polarization.
- **Per-bin RingZT and RingZLab.** The tree has no per-bin sums of them, so their curves keep the producer's pooled-master errors and have no paired difference. Their integrals are fully jackknifed and paired.