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
- `pyvista` optional for presentation renders.

## Layout

Every file of the toy, relative to `ToyModels/RingGeometryToy/`:

```
.gitignore             ignores untrackedOutput/ and Python caches
run_ring_toy.py        command-line driver (presets, overrides, provenance)
ringtoy/
  __init__.py          package marker (required for `import ringtoy`)
  geometry.py          jet frame, ring families, placement
  flow.py              elliptic flow: shape map, momentum map, v2 calibration
  regimes.py           polarized ring states for regimes R0-R3
  observable.py        ring observable, angular variables, closed forms
  scenarios.py         configuration dataclasses, presets, overrides
  scans.py             regime evaluation and parameter scans
  plots2d.py           matplotlib figures
  viz3d.py             plotly scenes
tests/
  conftest.py          puts ToyModels/RingGeometryToy/ on sys.path for pytest
  test_geometry.py     frames, sign conventions, ring families
  test_flow.py         v2 calibration, affine tangent transport
  test_regimes.py      regime reductions, rotation covariance, bookkeeping
  test_closed_forms.py closed forms of the Ring Design Discussion note
untrackedOutput/       run outputs (created on first run, git-ignored)
```

`ringtoy/__init__.py` and `tests/conftest.py` must stay in those folders:
without the first, `ringtoy` becomes a namespace package; without the second,
pytest cannot import `ringtoy`.

## Quick start

```
cd ToyModels/RingGeometryToy
pytest -q tests
python run_ring_toy.py --list-presets
python run_ring_toy.py --preset regimes --tube-radius 0.12
python run_ring_toy.py --preset position --set shape.radius=1.5 --no-3d
python run_ring_toy.py --preset misaligned --set scan.n_jet_azimuths=90
```

Overrides use dotted keys and Python literals, for example
`--set flow.v2=0.1`, `--set shape.kind=ellipse`,
`--set "shape.fourier_modes=((2,0.2,0.0),)"`,
`--set "placement.center_lab=(2.0,0.0,0.0)"`. Requirements: numpy and
matplotlib; plotly for the scenes (skip them with `--no-3d`); scipy only for
one test, which is skipped without it.

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
(`strength`): $\vec x'=(1+sH)\vec x$. The map is affine, so tangents
transform exactly with the same matrix, $d\vec x'/d\lambda=(1+sH)\,d\vec x/d\lambda$.
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

## Regimes

| Regime | Ring shape | Momentum direction |
|---|---|---|
| R0 | $\vec x$ | $\hat x$ |
| R1 | $(1+sH) \vec x$ | $\hat x'$ |
| R2 | $\vec x$ | ${\rm normalize}(H\vec x)$ |
| R3 | $(1+sH) \vec x$ | ${\rm normalize}(H\vec x')$ |

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

## Presets

| Preset | Question |
|---|---|
| `regimes` | How R1-R3 change the ring at one configuration; all 2D figures, the 3D scene, and an R3 strength slider. |
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
* `*_scene.html`: position space (vertex, jet axis, event plane, rings,
  polarization arrows, optional cosmetic torus) next to the momentum sphere
  (directions coloured by $r$, $\eta$ window). Legend entries toggle regimes
  in both panels.
* `*_R3_strength_slider.html`: R3 with a slider over the shape-map strength.

Every run writes `provenance.json` with the git commit (with `-dirty`), user,
UTC timestamp, command line, full configuration, and the integrated results.

## Open items

* `TODO:` stretch scaling of $|\vec P|$ (see above).
* `TODO:` ring rarefaction and the choice of weights under deformation.
* The torus in the scenes is cosmetic; thickening does not enter any
  observable yet.
* Non-affine flows and thermal smearing between position and momentum are
  deferred.
