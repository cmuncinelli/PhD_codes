"""
notebook_ui.py -- Interactive explorer for Jupyter (ipywidgets + plotly FigureWidget).

Part of the ring geometry toy (README.md explains the controls). Every change
builds a configuration, calls views.compute_view, and refreshes the figure in
place through viz3d.update_figure: the widget adds controls, not physics.
"""

import dataclasses

import ipywidgets as widgets
import numpy as np
from IPython.display import display

from .geometry import RingPlacement, RingShape
from .observable import Acceptance
from .output import make_run_dir, write_provenance
from .regimes import REGIMES, PolarizationModel
from .scenarios import get_preset
from .views import SweepSpec, compute_view, readouts
from .viz3d import FigureOptions, animated_figure, make_figure, update_figure

# Range and number of points of each sweep offered in the explorer.
SWEEPS = {
    "jet_phi": (-np.pi, np.pi, 41),
    "jet_eta": (-1.5, 1.5, 31),
    "delta": (0.0, np.pi, 25),
    "placement.distance_along_jet": (-4.0, 4.0, 33),
    "flow.strength": (0.0, 2.5, 26),
    "flow.v2": (-0.4, 0.4, 33),
    "shape.radius": (0.3, 3.0, 28),
    "shape.tilt": (-1.2, 1.2, 25),
}


def _slider(description, value, low, high, step=0.01):
    return widgets.FloatSlider(value=value, min=low, max=high, step=step, description=description,
                               continuous_update=False, readout_format=".2f",
                               style=dict(description_width="120px"), layout=widgets.Layout(width="360px"))


def _play(slider, interval_ms=150):
    """Play button that steps a slider through its range and loops."""
    n_steps = int(round((slider.max - slider.min) / slider.step))
    play = widgets.Play(min=0, max=n_steps, step=max(1, n_steps // 60), interval=interval_ms,
                        layout=widgets.Layout(width="140px"))

    def on_step(change):
        slider.value = slider.min + change["new"] * slider.step

    play.observe(on_step, names="value")
    return play


class RingExplorer:
    """Sliders for every toy parameter, a live scene, and export buttons."""

    def __init__(self, config=None, sweep_parameter="jet_phi", output_root=None):
        self.base = config if config is not None else get_preset("regimes")
        self.output_root = output_root
        self._busy = False
        c = self.base

        # ---- jet ----
        self.jet_eta = _slider("jet eta", c.jet_eta, -1.5, 1.5)
        self.jet_phi = _slider("jet phi [rad]", c.jet_phi, -np.pi, np.pi)
        # ---- ring ----
        self.kind = widgets.Dropdown(options=["circle", "ellipse", "fourier"], value=c.shape.kind,
                                     description="ring kind", style=dict(description_width="120px"))
        self.radius = _slider("radius [fm]", c.shape.radius, 0.2, 4.0)
        self.radius_minor = _slider("minor radius [fm]", c.shape.radius_minor, 0.2, 4.0)
        self.ellipse_angle = _slider("ellipse angle", c.shape.ellipse_angle, -np.pi, np.pi)
        self.fourier_m = widgets.IntSlider(value=2, min=2, max=6, description="Fourier m",
                                           continuous_update=False, style=dict(description_width="120px"))
        self.fourier_amp = _slider("Fourier amplitude", 0.2, 0.0, 0.8)
        self.fourier_phase = _slider("Fourier phase", 0.0, -np.pi, np.pi)
        self.tilt = _slider("tilt [rad]", c.shape.tilt, -1.5, 1.5)
        self.tilt_axis = _slider("tilt axis angle", c.shape.tilt_axis_angle, -np.pi, np.pi)
        self.circulation = widgets.Dropdown(options=[("right-handed", 1), ("reversed", -1)],
                                            value=c.shape.circulation, description="circulation",
                                            style=dict(description_width="120px"))
        # ---- placement ----
        self.placement_mode = widgets.Dropdown(options=["along jet", "fixed lab centre"],
                                               value="along jet" if c.placement.center_lab is None
                                               else "fixed lab centre",
                                               description="placement", style=dict(description_width="120px"))
        self.distance = _slider("L along jet [fm]", c.placement.distance_along_jet, -5.0, 5.0)
        self.offset_e1 = _slider("offset e1 [fm]", c.placement.offset_e1, -3.0, 3.0)
        self.offset_e2 = _slider("offset e2 [fm]", c.placement.offset_e2, -3.0, 3.0)
        centre = c.placement.center_lab if c.placement.center_lab is not None else (2.0, 0.0, 0.0)
        self.centre_x = _slider("centre x [fm]", centre[0], -5.0, 5.0)
        self.centre_y = _slider("centre y [fm]", centre[1], -5.0, 5.0)
        self.centre_z = _slider("centre z [fm]", centre[2], -5.0, 5.0)
        # ---- flow and almond ----
        self.v2 = _slider("v2", c.flow.v2, -0.6, 0.6)
        self.psi2 = _slider("Psi2 [rad]", c.flow.psi2, -np.pi / 2, np.pi / 2)
        self.strength = _slider("strength s", c.flow.strength, 0.0, 3.0)
        self.h_z = _slider("h_z", c.flow.h_z, 0.2, 3.0)
        self.almond_show = widgets.Checkbox(value=c.almond.show, description="show almond")
        self.kappa2 = _slider("kappa2", c.almond.kappa2, 0.1, 1.0)
        self.almond_radius = _slider("almond rms [fm]", c.almond.rms_radius, 0.5, 4.0)
        self.almond_evolve = widgets.Checkbox(value=c.almond.evolve_with_flow, description="evolve with flow")
        # ---- polarization and acceptance ----
        self.magnitude = _slider("|P| = f", c.polarization.magnitude, 0.0, 1.0)
        self.stretch = widgets.Checkbox(value=c.polarization.stretch_scaling, description="stretch scaling")
        self.weight_mode = widgets.Dropdown(options=["material", "deformed_arclength"],
                                            value=c.polarization.weight_mode, description="weights",
                                            style=dict(description_width="120px"))
        self.use_window = widgets.Checkbox(value=c.acceptance.eta_max is not None, description="eta window")
        self.eta_max = _slider("eta max", c.acceptance.eta_max if c.acceptance.eta_max is not None else 0.5,
                               0.1, 2.0)
        # ---- display ----
        self.regimes = widgets.SelectMultiple(options=REGIMES, value=REGIMES, description="regimes",
                                              rows=4, style=dict(description_width="120px"))
        self.quantity = widgets.Dropdown(options=["signal_density", "ring_average"], value="signal_density",
                                         description="side plots", style=dict(description_width="120px"))
        self.tube = _slider("torus radius [fm]", 0.0, 0.0, 0.5)
        self.sweep_parameter = widgets.Dropdown(options=list(SWEEPS), value=sweep_parameter,
                                                description="sweep", style=dict(description_width="120px"))
        self.highlight = widgets.IntSlider(value=-1, min=-1, max=c.binning.n_delta_phi - 1,
                                           description="highlight bin", continuous_update=False,
                                           style=dict(description_width="120px"))
        # ---- actions ----
        self.export_scene = widgets.Button(description="export scene HTML")
        self.export_movie = widgets.Button(description="export movie of sweep")
        self.status = widgets.HTML()
        self.table = widgets.HTML()

        self.controls = [self.jet_eta, self.jet_phi, self.kind, self.radius, self.radius_minor,
                         self.ellipse_angle, self.fourier_m, self.fourier_amp, self.fourier_phase, self.tilt,
                         self.tilt_axis, self.circulation, self.placement_mode, self.distance, self.offset_e1,
                         self.offset_e2, self.centre_x, self.centre_y, self.centre_z, self.v2, self.psi2,
                         self.strength, self.h_z, self.almond_show, self.kappa2, self.almond_radius,
                         self.almond_evolve, self.magnitude, self.stretch, self.weight_mode, self.use_window,
                         self.eta_max, self.regimes, self.quantity, self.tube, self.sweep_parameter,
                         self.highlight]
        for control in self.controls:
            control.observe(self._on_change, names="value")
        self.export_scene.on_click(self._on_export_scene)
        self.export_movie.on_click(self._on_export_movie)

        view, options = self._current()
        self.figure = make_figure(view, options, widget=True)
        self._built_for = (self.sweep_parameter.value, self.quantity.value)
        self._layout = None
        self._connect_clicks()
        self._show_readouts(view)

    # ---- configuration from the controls ----
    def config(self):
        """The ToyConfig currently described by the controls."""
        if self.kind.value == "fourier":
            modes = ((int(self.fourier_m.value), float(self.fourier_amp.value), float(self.fourier_phase.value)),)
        else:
            modes = ()
        shape = RingShape(kind=self.kind.value, radius=self.radius.value, radius_minor=self.radius_minor.value,
                          ellipse_angle=self.ellipse_angle.value, fourier_modes=modes, tilt=self.tilt.value,
                          tilt_axis_angle=self.tilt_axis.value, circulation=int(self.circulation.value),
                          n_points=self.base.shape.n_points)
        if self.placement_mode.value == "along jet":
            placement = RingPlacement(distance_along_jet=self.distance.value, offset_e1=self.offset_e1.value,
                                      offset_e2=self.offset_e2.value)
        else:
            placement = RingPlacement(center_lab=(self.centre_x.value, self.centre_y.value, self.centre_z.value))
        return dataclasses.replace(
            self.base, jet_eta=self.jet_eta.value, jet_phi=self.jet_phi.value, shape=shape, placement=placement,
            flow=dataclasses.replace(self.base.flow, v2=self.v2.value, psi2=self.psi2.value,
                                     strength=self.strength.value, h_z=self.h_z.value),
            almond=dataclasses.replace(self.base.almond, show=self.almond_show.value, kappa2=self.kappa2.value,
                                       rms_radius=self.almond_radius.value,
                                       evolve_with_flow=self.almond_evolve.value),
            polarization=PolarizationModel(magnitude=self.magnitude.value, stretch_scaling=self.stretch.value,
                                           weight_mode=self.weight_mode.value),
            acceptance=Acceptance(eta_max=self.eta_max.value if self.use_window.value else None))

    def sweep(self):
        low, high, n = SWEEPS[self.sweep_parameter.value]
        return SweepSpec(parameter=self.sweep_parameter.value, values=tuple(np.linspace(low, high, n)))

    def options(self):
        return FigureOptions(regimes_shown=tuple(self.regimes.value), quantity=self.quantity.value,
                             tube_radius=self.tube.value,
                             highlight_bin=None if self.highlight.value < 0 else int(self.highlight.value))

    def _current(self):
        return compute_view(self.config(), sweep=self.sweep()), self.options()

    # ---- callbacks ----
    def _on_change(self, _change):
        # Play buttons can fire faster than a refresh; later events are
        # dropped while one refresh is running rather than queued.
        if self._busy:
            return
        self._busy = True
        try:
            view, options = self._current()
            if self._needs_rebuild():
                self.figure = make_figure(view, options, widget=True)
                self._connect_clicks()
                if self._layout is not None:
                    self._layout.children = self._layout.children[:-1] + (self.figure,)
            else:
                update_figure(self.figure, view, options)
            self._show_readouts(view)
            self.status.value = ""
        except ValueError as error:
            # Invalid combinations (e.g. |v2/kappa2| >= 1) are reported, not raised.
            self.status.value = f"<b style='color:#cc0000'>{error}</b>"
        finally:
            self._busy = False

    def _needs_rebuild(self):
        # Subplot titles and axis labels depend on the sweep and on the quantity.
        wanted = (self.sweep_parameter.value, self.quantity.value)
        rebuild = self._built_for != wanted
        self._built_for = wanted
        return rebuild

    def _connect_clicks(self):
        # Clicking a Delta phi curve highlights the ring elements of that bin.
        edges = np.linspace(-np.pi, np.pi, self.base.binning.n_delta_phi + 1)

        def on_click(_trace, points, _selector):
            if points.xs:
                index = int(np.clip(np.searchsorted(edges, points.xs[0], side="right") - 1, 0, len(edges) - 2))
                self.highlight.value = index

        for trace in self.figure.data:
            if trace.type == "scatter" and trace.xaxis == "x":
                trace.on_click(on_click)

    def _show_readouts(self, view):
        rows = "".join(f"<tr><td><b>{regime}</b></td><td>{r['ring_average']:.4f}</td>"
                       f"<td>{r['signal_per_total_weight']:.4f}</td><td>{r['acceptance_fraction']:.3f}</td>"
                       f"<td>{r['mean_cos_gamma']:.4f}</td></tr>"
                       for regime, r in readouts(view).items())
        almond = (f"almond eps2 = {view.almond_eccentricity:+.4f}" if view.almond_contours else "almond hidden")
        self.table.value = ("<table><tr><th></th><th>&lt;P.n&gt;</th><th>sum w r / W</th><th>acceptance</th>"
                            f"<th>&lt;cos gamma&gt;</th></tr>{rows}</table><p>{almond}</p>")

    def _on_export_scene(self, _button):
        view, options = self._current()
        run_dir = make_run_dir("explorer_scene", self.output_root)
        make_figure(view, options).write_html(f"{run_dir}/scene.html", include_plotlyjs=True)
        write_provenance(run_dir, "explorer_scene", view.config, readouts(view), "notebook_ui.export_scene")
        self.status.value = f"scene written to {run_dir}"

    def _on_export_movie(self, _button):
        self.status.value = "rendering movie..."
        config, options, sweep = self.config(), self.options(), self.sweep()
        run_dir = make_run_dir("explorer_movie", self.output_root)
        fig = animated_figure(config, sweep.parameter, sweep.values, options, sweep=sweep)
        fig.write_html(f"{run_dir}/movie.html", include_plotlyjs=True, auto_play=False)
        write_provenance(run_dir, "explorer_movie", config,
                         {"parameter": sweep.parameter, "values": list(sweep.values)},
                         "notebook_ui.export_movie")
        self.status.value = f"movie written to {run_dir}"

    # ---- layout ----
    def display(self):
        def box(title, children):
            return widgets.VBox([widgets.HTML(f"<b>{title}</b>")] + children)

        tabs = widgets.Tab(children=[
            box("jet", [self.jet_eta, widgets.HBox([self.jet_phi, _play(self.jet_phi)])]),
            box("ring", [self.kind, self.radius, self.radius_minor, self.ellipse_angle, self.fourier_m,
                         self.fourier_amp, self.fourier_phase, self.tilt, self.tilt_axis, self.circulation]),
            box("placement", [self.placement_mode, widgets.HBox([self.distance, _play(self.distance)]),
                              self.offset_e1, self.offset_e2, self.centre_x, self.centre_y, self.centre_z]),
            box("flow and almond", [self.v2, self.psi2, widgets.HBox([self.strength, _play(self.strength)]),
                                    self.h_z, self.almond_show, self.kappa2, self.almond_radius,
                                    self.almond_evolve]),
            box("polarization", [self.magnitude, self.stretch, self.weight_mode, self.use_window, self.eta_max]),
            box("display", [self.regimes, self.quantity, self.tube, self.sweep_parameter, self.highlight,
                            widgets.HBox([self.export_scene, self.export_movie])]),
        ])
        for index, title in enumerate(["jet", "ring", "placement", "flow/almond", "polarization", "display"]):
            tabs.set_title(index, title)
        self._layout = widgets.VBox([widgets.HBox([tabs, self.table]), self.status, self.figure])
        display(self._layout)
