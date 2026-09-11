# RingGeometryToy

> WARNING! This was essentially AI generated and is a simple tool to visualize the distorted rings and some geometrical effects. This is NOT meant as an actual physics check, just a qualitative tool for visualization and some intuition. If any of these results is ever used in a presentation, then that result will be throughly QAed beforehand.

> This code may not be ready as-is!

A geometric toy for the response of the ring observable to idealized and
distorted vortex rings. It accompanies the note *Ring Design Discussion* and
answers geometric questions only: how the shape of a polarized ring, its
position relative to the vertex, and a simple elliptic flow change what the
ring observable measures.

What it is **not**: there are no $\Lambda$ hyperons, no decays, no statistics,
and no hydrodynamics. Each ring element carries a prescribed polarization
vector and a momentum direction, so every observable is an exact weighted mean
over ring elements. Comparisons with the model papers (Lisa et al. 2021,
Serenone et al. 2021, Ribeiro et al. 2024, and V. H. Ribeiro's dissertation)
are qualitative by construction.

## Dependencies:

To run this code, you will need to install the following packages:
- `numpy`;
- `scipy`;
- `matplotlib`;
- `plotly`;
- `pytest`;
- `anywidget` (for dynamic ipynb widgets);
- `pyvista` optional for presentation renders.

## Layout

Every file of the toy, relative to `ToyModels/RingGeometryToy/`:

```
.gitignore             ignores untrackedOutput/ and Python caches
RingGeometryToy.ipynb  interactive explorer and explorations (Jupyter)
run_ring_toy.py        command-line driver (presets, explorations, overrides)
ringtoy/
  __init__.py          package marker (required for `import ringtoy`)
  geometry.py          jet frame, ring families, placement
  flow.py              elliptic flow: shape map, momentum map, v2 calibration
  almond.py            QGP-like almond consistent with v2 (context only)
  regimes.py           polarized ring states for regimes R0-R3
  observable.py        ring observable, angular variables, closed forms
  scenarios.py         configuration dataclasses, presets, overrides
  scans.py             regime evaluation and parameter scans
  views.py             everything a scene shows, without plotting code
  viz3d.py             plotly scenes, widget updates and movies from one trace builder
  plots2d.py           matplotlib figures
  explorations.py      movies and checks answering the toy's questions
  notebook_ui.py       ipywidgets explorer
  output.py            output location and provenance
tests/
  conftest.py          puts ToyModels/RingGeometryToy/ on sys.path for pytest
  test_geometry.py     frames, sign conventions, ring families
  test_flow.py         v2 calibration, affine tangent transport
  test_almond.py       almond eccentricity, contours, inversion
  test_regimes.py      regime reductions, rotation covariance, bookkeeping
  test_closed_forms.py closed forms of the note and of the ellipse average
  test_views.py        sweeps, trace bookkeeping of widgets and movies
```

`ringtoy/__init__.py` and `tests/conftest.py` must stay in those folders:
without the first, `ringtoy` becomes a namespace package; without the second,
pytest cannot import `ringtoy`.

**One implementation of the physics.** `views.compute_view` evaluates a
configuration; `viz3d.trace_data` turns the result into trace data, always in
the same order. The static scenes, the in-place updates of the notebook
widget, and every movie frame are built from these two functions.

## Quick start

```
cd ToyModels/RingGeometryToy
pytest -q tests
python run_ring_toy.py --list
python run_ring_toy.py --preset regimes --tube-radius 0.12
python run_ring_toy.py --preset position --set shape.radius=1.5 --no-3d
python run_ring_toy.py --explore morph
jupyter lab --no-browser --port 8890        # then open RingGeometryToy.ipynb
```

Overrides use dotted keys and Python literals, for example
`--set flow.v2=0.1`, `--set shape.kind=ellipse`,
`--set "shape.fourier_modes=((2,0.2,0.0),)"`,
`--set "placement.center_lab=(2.0,0.0,0.0)"`, `--set almond.kappa2=0.25`.

Requirements: numpy and matplotlib; plotly for scenes and movies (skip them
with `--no-3d`); scipy for the ellipticity check and two tests (skipped
without it); `jupyterlab`, `ipywidgets` and `anywidget` for the notebook.

### Where outputs go

Each run creates `<output root>/<label>_<UTC timestamp>/`. The output root is
`--output-dir` when given, otherwise `DEFAULT_OUTPUT_DIR` in
`ringtoy/output.py`, currently `/home/users/cicerodm/RingPol/RingGeometryToy`.
The notebook passes `OUTPUT_ROOT` (default `None`, meaning the same default).
Missing directories are created. The git commit recorded in `provenance.json`
is always the commit of this code, wherever the outputs are written.

### Running the notebook on jarvis15 over ssh (just use vscode, duh...)

```
# on jarvis15, from ToyModels/RingGeometryToy/
jupyter lab --no-browser --port 8890
# on the laptop
ssh -N -L 8890:localhost:8890 cicerodm@jarvis15
```

and open the `http://localhost:8890/lab?token=...` link printed by Jupyter.

## Tests

The tests check the toy against results that are known independently of the
code: analytic geometry, the derivations in this README, and the closed forms
of the *Ring Design Discussion* note. They use generic, non-symmetric
parameter values, because symmetric choices can hide sign and factor errors.

### Running them and reading the output

```
pytest -q tests                    # compact: one character per test
pytest -v tests                    # one line per test, with its name
pytest -q tests -k elliptic        # only tests whose name matches
```

In the compact output `.` is a pass, `F` a failed assertion, `E` an error
while running or collecting a test, and `s` a skipped test (the scipy-based
tests without scipy, `test_views.py` without plotly). A healthy run ends with
`41 passed`. On a failure, pytest prints the assertion with the actual and
expected values; the test name and the table below say which part of the toy
is affected.

The tests write no result files. The only by-products are `.pytest_cache/` and
`__pycache__/`, both git-ignored and safe to delete. There is nothing to keep:
rerun the suite after any change to `ringtoy/`.

### What each test protects

`tests/test_geometry.py` -- frames, conventions and ring families

| Test | Checks | Protects against |
|---|---|---|
| `test_jet_frame_is_right_handed_orthonormal` | $(e_1,e_2,\hat t)$ is orthonormal and $e_1\times e_2=\hat t$ | a left-handed frame, which flips the apparent circulation in every plot |
| `test_sign_convention_right_handed_ring` | jet along $+\hat x$ and $\hat p=+\hat y$ give $\hat n=+\hat z$; the default ring gives $r=+f$ | a global sign error of the observable |
| `test_coaxial_circle_gives_full_ring_value` (3 cases) | a circle centred on the jet gives $r=f$ at every element for $L>0$, $L<0$, $L=0$; weights add up to $2\pi a$ | errors in placement, in $\hat n$, or in the weights |
| `test_ellipse_ring_value_is_rho_dpsi_ds` | $r=\rho\,d\psi/ds=ab/(\rho\,ds/d\lambda)$ element by element | wrong tangents or wrong in-plane rotation of non-circular rings |
| `test_reversed_circulation_flips_ring_value` | `circulation=-1` gives exactly $-r$ | a circulation flag that moves points instead of reversing tangents |
| `test_fourier_ring_rejects_non_positive_radius` | a modulation with $\rho\le0$ raises | silently self-intersecting rings |

`tests/test_flow.py` -- elliptic flow

| Test | Checks | Protects against |
|---|---|---|
| `test_v2_calibration_is_exact` | emitters uniform in azimuth acquire $\langle\cos2(\phi_p-\Psi_2)\rangle=v_2$ to $10^{-12}$, with a longitudinal component present | a wrong $r(v_2)$ relation, a wrong event-plane rotation, or longitudinal leakage |
| `test_isotropic_flow_is_identity` | $v_2=0$, $h_z=1$ give $H=1$ for any $\Psi_2$ | spurious anisotropy from the rotation matrices |
| `test_affine_tangent_transport_matches_finite_differences` | $(1+sH)\,d\vec x/d\lambda$ equals the numerical derivative of the deformed curve | tangents that do not follow the deformed ring |
| `test_shape_map_rejects_folding` | strengths that fold the ring through the vertex raise | unphysical deformations passing silently |

`tests/test_almond.py` -- the almond drawn for context

| Test | Checks | Protects against |
|---|---|---|
| `test_toy_relation_gives_flow_anisotropy` | at $\kappa_2=1$, $a^2/b^2=r$, $\varepsilon_2=v_2$, and the rms radius is kept | an almond inconsistent with the flow |
| `test_kappa2_scales_eccentricity_and_limits_it` | $\varepsilon_2=v_2/\kappa_2$; $|\varepsilon_2|\ge1$ raises | silently impossible almonds |
| `test_contour_is_rotated_ellipse` | contour points satisfy the ellipse equation in the event-plane frame | a wrong orientation or wrong axes |
| `test_evolving_almond_becomes_round_at_inversion_strength` | the evolved almond is round at $s=(b-a)/(a-br)$ and in-plane beyond | a wrong shape map on the almond |

`tests/test_views.py` -- sweeps and the bookkeeping of widgets and movies (skipped without plotly)

| Test | Checks | Protects against |
|---|---|---|
| `test_parameter_round_trip` (4 cases) | `configure_parameter` and `get_parameter` are inverse, including `delta` | sweep markers at the wrong position |
| `test_sweep_does_not_depend_on_current_value_of_swept_parameter` | the cached sweep ignores where the swept parameter currently is | a sweep curve that changes during its own movie |
| `test_trace_layout_is_invariant` | every configuration and display option gives the same traces in the same order | widget updates writing data into the wrong traces |
| `test_update_figure_matches_fresh_figure` | an updated figure equals a freshly built one | stale data left behind by in-place updates |
| `test_animation_has_one_frame_per_value_with_all_traces` | movies have one complete frame per value | broken or partial movies |

`tests/test_regimes.py` -- regimes and bookkeeping

| Test | Checks | Protects against |
|---|---|---|
| `test_regimes_reduce_to_reference` (R1, R2, R3) | with trivial flow each regime reproduces R0 exactly | regimes that deform or deflect when they should not |
| `test_rotation_covariance` (R0-R3) | rotating the jet and $\Psi_2$ together about the beam leaves $r$, $\Delta\phi$, $\Delta\theta_{{\rm Jet},\Lambda}$ and the integrated results unchanged | hidden dependence on absolute lab angles |
| `test_signal_density_integrates_to_total_signal` | $\sum_{\rm bins}$ density $\times$ width $=\sum_{\rm acc}w\,r/W_{\rm tot}$ | normalization errors of the differential curves |
| `test_stretch_scaling_follows_line_stretching` | with stretch scaling, $|\vec P|=f\,|F\tau|/|\tau|$ | a stretch factor applied with the wrong ratio |

`tests/test_closed_forms.py` -- results of the *Ring Design Discussion* note

| Test | Checks | Protects against |
|---|---|---|
| `test_cos_c_sin_c_closed_form_matches_vector_algebra` | the closed forms of $\cos C$ and $\sin C$ against vector algebra at 200 random configurations | errors in the note's proposition or in `ring_directions` |
| `test_elliptic_average` (4 values of $\eta_{\rm Jet}$) | $\frac{2}{\pi}\tanh\eta_{\rm Jet}K(\mathrm{sech}\,\eta_{\rm Jet})$ against numerical integration | the $\tanh$-like fake-ring shape of the note |
| `test_ellipse_ring_average_matches_toy` (2 cases) | the toy ring average of a coaxial ellipse equals $(b/a)K(k)/E(k)$, with $a>b$ and $b>a$ | errors in weights or tangents of non-circular rings |
| `test_quadrupole_leakage_selection_rule` | $\langle n_z\sin2(u+\delta)\rangle$ vanishes when $\eta_{\rm Jet}=0$ or $\eta_\Lambda=0$ and scales as $\cos2\delta$ otherwise | the Test 4 prediction of the note ($v_2$-polarization leakage) |

## Conventions

* Lab frame: beam along $\hat z$, vertex at the origin, lengths in fm, doubles
  throughout, Cartesian vectors.
* Jet: a direction $\hat t(\eta_{\rm Jet},\phi_{\rm Jet})$ only. The jet frame
  $(e_1,e_2,\hat t)$ is right-handed, with $e_2$ the projection of the beam
  axis onto the plane transverse to the jet ("up") and $e_1=e_2\times\hat t$
  ("right"). A jet pointing at the viewer therefore sees right-handed
  circulation as counterclockwise in the $(e_1,e_2)$ plane.
* Ring direction and observable, per element:
  $\hat n=\hat t\times\hat p/|\hat t\times\hat p|$ and $r=\vec P\cdot\hat n$.
  No $3/\alpha$ factor appears because $\vec P$ is prescribed, not
  reconstructed from decays. An ideal coaxial ring gives $r=f$ exactly.
* Opening angle: $\Delta\theta_{{\rm Jet},\Lambda}$ between $\hat t$ and
  $\hat p$ (code: `delta_theta`), computed with `atan2` for full precision
  near $0$ and $\pi$.
* Pseudorapidity of $\hat p$ stands in for rapidity: the toy has directions
  only, so there is no $p_T$ axis.

## The ring

The core curve is built in the jet frame, centred at $\vec D$, with
right-handed circulation about $\hat t$ (`circulation=-1` reverses it):

* `circle` (radius $a$), `ellipse` (semi-axes and in-plane angle), and
  `fourier`, with $\rho(\lambda)=a\,[1+\sum_m A_m\cos m(\lambda-\lambda_m)]$;
* `tilt` rotates the ring plane out of the plane transverse to the jet;
* placement: `distance_along_jet` $L$ ($L>0$ outward, $L<0$ inward),
  transverse offsets, or a fixed lab centre `center_lab` (jet--flow
  misalignment, where the centre stays put while the jet direction varies).

Material points sit on a uniform grid in the curve parameter with analytic
tangents. The curves are smooth and periodic, so ring averages on this grid
are spectrally accurate.

Why emission from the vertex matters: with $\hat p=\hat x$ and
$\vec D=L\hat t$, the ring is seen from the vertex as a cone of half-angle
$\arctan(a/|L|)$ around $+\hat t$ for $L>0$ and around $-\hat t$ for $L<0$,
and in both cases $\hat t\times\vec x=a\,\hat t\times\hat\rho=a\,\hat\psi$ is the
ring tangent. Outward rings therefore give narrow double peaks around
$\Delta\phi=0$, inward rings give peaks near $\pm\pi$, and both are positive,
which is the pattern of the insertion-position studies of the model papers.

## Elliptic flow

The flow is the linear velocity field $\vec u(\vec x)=H\vec x$ about the
vertex, an anisotropic version of Hubble's law:

$$H=R_z(\Psi_2)\,{\rm diag}(1,\;r,\;h_z)\,R_z(-\Psi_2),\qquad r=\frac{1-v_2}{1+v_2}.$$

Along the event plane the expansion rate is $1$, perpendicular to it $r$, and
along the beam $h_z$. It is used in two independent ways.

**Shape map.** Each ring point moves along the flow for a "time" $s$
(`strength`): $\vec x\,'=(1+sH)\,\vec x$. The map is affine, so tangents
transform exactly with the same matrix, $d\vec x\,'/d\lambda=(1+sH)\,d\vec x/d\lambda$.
A circle around the vertex becomes an ellipse elongated along the event
plane; an off-centre ring is stretched, rotated and pushed outward. Strengths
for which $1+sH$ has a non-positive eigenvalue are rejected, since they would
fold the ring through the vertex.

**Momentum map.** Each element emits along its local flow velocity,
$\hat p={\rm normalize}(H\vec x)$. For $H\propto 1$ this is $\hat p=\hat x$.

**Calibration of $v_2$.** For emitters uniform in azimuth about the vertex,
the transverse momentum azimuth in the event-plane frame satisfies
$\tan\phi_p'=r\tan\phi_x'$. The resulting density is
$dN/d\phi_p'=\frac{1}{2\pi}\,r/(r^2\cos^2\phi_p'+\sin^2\phi_p')$. Writing the
denominator as $A+B\cos2\phi_p'$ with $A=(r^2+1)/2$ and $B=(r^2-1)/2$, and using
$\int_0^{2\pi}\cos2\phi\,/(A+B\cos2\phi)\,d\phi=(2\pi/B)(1-A/\sqrt{A^2-B^2})$
with $\sqrt{A^2-B^2}=r$, gives

$$\langle\cos2(\phi_p-\Psi_2)\rangle=\frac{r-A}{B}=\frac{1-r}{1+r}=v_2 .$$

The longitudinal component does not enter, because $H$ does not mix it with the
transverse plane. `test_v2_calibration_is_exact` checks this at generic values.

## The almond

`almond.py` draws a QGP-like almond for context: it does not enter any
observable. A Gaussian almond with rms widths $a$ along the event plane and
$b$ perpendicular to it has pressure gradients $\propto(x/a^2,\,y/b^2)$, so the
fluid is accelerated with the anisotropy ratio $r=a^2/b^2$ of the elliptic
flow. With $r=(1-v_2)/(1+v_2)$,

$$\varepsilon_2=\frac{b^2-a^2}{b^2+a^2}=\frac{(1+v_2)-(1-v_2)}{(1+v_2)+(1-v_2)}=v_2 ,$$

that is, a response coefficient $\kappa_2=1$ in this toy. Realistic
hydrodynamics gives $\kappa_2\sim0.2$-$0.3$, so `kappa2` is a parameter and
$\varepsilon_2=v_2/\kappa_2$. Given the rms radius $R=\sqrt{(a^2+b^2)/2}$,
$a=R\sqrt{1-\varepsilon_2}$ and $b=R\sqrt{1+\varepsilon_2}$; for $v_2>0$ the long axis
is perpendicular to $\Psi_2$. Scenes show the $1\sigma$ (filled) and $2\sigma$
contours at $z=0$.

With `evolve_with_flow`, the almond follows the same shape map as the ring:
its axes scale as $a(1+s)$ and $b(1+sr)$, it becomes round at
$s=(b-a)/(a-br)$, and it is elongated in plane beyond that (eccentricity
inversion). The current eccentricity is quoted in scene titles and in the
explorer.

## Regimes

| Regime | Ring shape | Momentum direction |
|---|---|---|
| R0 | $\vec x$ | $\hat x$ |
| R1 | $(1+sH)\,\vec x$ | $\hat x\,'$ |
| R2 | $\vec x$ | ${\rm normalize}(H\vec x)$ |
| R3 | $(1+sH)\,\vec x$ | ${\rm normalize}(H\vec x\,')$ |

R1 isolates the deformation of the ring, R2 the divergence between the ring
and the particle directions, and R3 combines them (the ring is first carried
by the flow and then emits along the flow at its new position). Jet--flow
misalignment is a placement choice and combines with any regime.

## Polarization and weights

Each element carries $\vec P=f\,\hat\tau$ with $f\in[0,1]$ and $\hat\tau$ the
unit tangent of the (possibly deformed) ring. Two choices are open and are
kept as options, both flagged with `TODO:` in `regimes.py`:

* `stretch_scaling` (default off): scale $|\vec P|$ by the line stretching
  $|(1+sH)\,\tau|/|\tau|$, as vortex stretching would.
* `weight_mode` (default `material`): `material` keeps the weight of each
  element fixed under the shape map (emitters conserved);
  `deformed_arclength` re-weights by the deformed length. This is where ring
  rarefaction will be revisited.

## Observables

* **Integrated** (`RingSummary`): the ring average $\langle\vec P\cdot\hat n\rangle$,
  the parity-forbidden components $\langle\vec P\cdot\hat\vartheta\rangle$ and
  $\langle\vec P\cdot\hat p\rangle$ (with $(\hat p,\hat\vartheta,\hat n)$
  right-handed), the acceptance fraction, and
  $\sum_{\rm acc}w\,r/W_{\rm tot}$. Results are kept as additive sums
  (`IntegratedSums`) and ratios are formed last, so rings and jet orientations
  can be combined correctly.
* **Differential** in $\Delta\phi$ and $\Delta\theta_{{\rm Jet},\Lambda}$
  (`BinnedRing`), in two forms:
  * *ring average per bin*: how well the ring elements in that bin align with
    $\hat n$;
  * *signal density* $\sum w\,r/(W_{\rm tot}\,\Delta x)$. If a dominant
    background of unpolarized emitters is flat in the binned variable, the
    measured ring value per bin is proportional to this density, so it plays
    the role of the $\mathcal R(\Delta\phi)$ curves of the model papers without
    simulating any background. Caveat: under the momentum map a background
    emitted from the vertex is no longer flat in $\Delta\phi$, so for R2 and
    R3 the correspondence is approximate.
* **Divergence diagnostic**: the distribution of
  $\cos\gamma=\hat P\cdot\hat n$, which is $1$ for a perfectly aligned ring.

## Interactive explorer (notebook)

`RingExplorer` in `notebook_ui.py` shows the scene of `viz3d` as a plotly
FigureWidget next to controls grouped in tabs:

* **jet**: $\eta_{\rm Jet}$, $\phi_{\rm Jet}$ (with play button);
* **ring**: kind, radius, minor radius and angle, one Fourier mode, tilt and
  tilt axis, circulation;
* **placement**: along the jet ($L$ with play button, transverse offsets) or a
  fixed lab centre;
* **flow/almond**: $v_2$, $\Psi_2$, strength $s$ (with play button), $h_z$;
  almond visibility, $\kappa_2$, rms radius, evolution with the flow;
* **polarization**: $f$, stretch scaling, weight mode, $\eta$ window;
* **display**: regimes, side-plot quantity, cosmetic torus, sweep parameter,
  highlighted $\Delta\phi$ bin, export buttons.

The figure has the position space (vertex, jet axis, event plane, almond,
rings with polarization arrows coloured by $r$), the momentum sphere (ring
elements coloured by $r$, $\eta$ window), and three side plots: the chosen
quantity versus $\Delta\phi$ and $\Delta\theta_{{\rm Jet},\Lambda}$, and the
integrated value along the chosen sweep with the current value marked. A
table gives $\langle\vec P\cdot\hat n\rangle$, $\sum w\,r/W$, the acceptance and
$\langle\cos\gamma\rangle$ per regime. Clicking a $\Delta\phi$ curve (or the
*highlight bin* slider) marks the ring elements of that bin in both 3D panels.
Invalid combinations, such as $|v_2/\kappa_2|\ge1$, are reported below the
controls instead of raising. Each refresh takes about 0.05-0.15 s; sweeps are
cached, so playing the swept parameter does not recompute its curve.

## Explorations

`explorations.py` holds the studies below, runnable from the notebook or with
`python run_ring_toy.py --explore <name>`. Each writes a run directory with
provenance.

| Name | What it shows | Question |
|---|---|---|
| `morph` | movie over $L$ from $-4$ to $4$ fm | inward rings peak near $\Delta\phi=\pm\pi$ on a cone around $-\hat t$, outward rings near 0 around $+\hat t$ |
| `misalignment` | movie over $\phi_{\rm Jet}$ with the ring centre fixed at $x=2$ fm | how much signal jet--flow misalignment removes, compared with the model papers |
| `eta_jet_movie` | movie over $\eta_{\rm Jet}$ | how the $\eta$ window shapes the response |
| `evolution` | movie over $s$ with ring and almond evolving | ring deformation next to the eccentricity inversion |
| `delta_test` | integrated response versus $\delta=\phi_{\rm Jet}-\Psi_2$ and its relative variation | whether selecting jets relative to $\Psi_2$ is worth it |
| `ellipticity` | toy versus $(b/a)K(k)/E(k)$ for coaxial ellipses | validation of non-circular rings, and a clean figure |

Movies are self-contained HTML files with play/pause buttons and a slider.
Every frame stores its trace data, so movies draw rings with at most 240
elements and weigh roughly 10-25 MB (about 5 MB of which is the embedded
plotly.js).

## Presets

| Preset | Question |
|---|---|
| `regimes` | How R1-R3 change the ring at one configuration; all 2D figures and a movie over the strength. |
| `position` | Distance of the ring along the jet, outward and inward ($L=-4\ldots4$ fm). |
| `misaligned` | Ring centre fixed at $x=2$ fm while the jet azimuth is scanned, compared with a ring centred on the jet at the same distance. |
| `eta_jet` | Integrated response versus $\eta_{\rm Jet}$ inside a fixed $\eta$ window. |
| `delta` | Integrated response versus $\delta=\phi_{\rm Jet}-\Psi_2$ (signal shape only). |
| `strength` | Tolerance to the shape-map strength. |

## Figures

* `*_delta_phi`, `*_delta_theta`: ring average (top) and signal density
  (bottom), as step histograms.
* `*_jet_plane`: polarization vectors in the plane transverse to the jet
  through the (deformed) ring centre, $(\vec P\cdot e_1,\vec P\cdot e_2)$,
  coloured by $r$; $\langle|\vec P\cdot\hat t|\rangle$ is quoted per panel.
* `*_momentum_map`: ring elements in $(\Delta\phi,\eta)$ of their momentum,
  coloured by $r$, with the $\eta$ window.
* `*_cos_gamma`: the divergence diagnostic.
* `scan_*`: integrated ring average, diluted signal, and acceptance fraction
  versus the scanned parameter.
* `base_scene.html` (scan presets) and `regimes_strength_movie.html`: the
  scene described under *Interactive explorer*, static or as a movie. Legend
  entries toggle regimes in all panels.

Every run writes `provenance.json` with the git commit (with `-dirty`), user,
UTC timestamp, command line, full configuration, and the integrated results.

## Open items

* `TODO:` stretch scaling of $|\vec P|$ (see above).
* `TODO:` ring rarefaction and the choice of weights under deformation.
* The torus in the scenes is cosmetic; thickening does not enter any
  observable yet.
* The almond is context only; $\kappa_2=1$ is the toy relation, not a
  hydrodynamic prediction.
* Non-affine flows and thermal smearing between position and momentum are
  deferred.