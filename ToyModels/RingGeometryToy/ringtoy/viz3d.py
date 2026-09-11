"""
viz3d.py -- Plotly figures of a View: static scenes, live widgets, and movies.

Part of the ring geometry toy (README.md describes the scenes).

One function, trace_data, turns a View into the data of every trace, always
in the same order and number. The static figure, the in-place updates of the
notebook widget, and every movie frame are built from it, so what is drawn
is defined once.
"""

import dataclasses
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .regimes import REGIMES
from .views import compute_view, highlight_mask, readouts, with_parameter

REGIME_COLOR = {"R0": "#000000", "R1": "#0033cc", "R2": "#cc0000", "R3": "#008800"}
REGIME_DASH = {"R0": "solid", "R1": "dash", "R2": "dashdot", "R3": "dot"}
RING_COLORSCALE = "RdBu_r"
ALMOND_FILL = "#d9a066"
TRACE_CLASSES = {"scatter3d": go.Scatter3d, "mesh3d": go.Mesh3d,
                 "surface": go.Surface, "scatter": go.Scatter}

LABEL_DPHI = "Delta phi = phi_Lambda - phi_Jet [rad]"
LABEL_DTHETA = "Delta theta_(Jet,Lambda) [rad]"


@dataclass(frozen=True)
class FigureOptions:
    """How a View is drawn (nothing here changes the physics).

    regimes_shown : regimes drawn; the others keep empty traces.
    quantity      : "signal_density" or "ring_average" for the side plots.
    tube_radius   : cosmetic torus radius [fm]; 0 draws none.
    max_arrows    : polarization arrows per ring.
    highlight_bin : Delta phi bin whose ring elements are highlighted, or None.
    scene_extent  : half-width of the position scene [fm]; None fits the view.
                    Movies fix it over all frames so that the camera does not jump.
    height        : figure height [px].
    """
    regimes_shown: tuple = REGIMES
    quantity: str = "signal_density"
    tube_radius: float = 0.0
    max_arrows: int = 36
    highlight_bin: object = None
    scene_extent: object = None
    height: int = 760


def _xyz(points):
    return dict(x=points[:, 0], y=points[:, 1], z=points[:, 2])


def _segments(starts, vectors):
    """One trace of segments start --> start + vector, separated by gaps."""
    n = len(starts)
    points = np.full((3 * n, 3), np.nan)
    points[0::3] = starts
    points[1::3] = starts + vectors
    return _xyz(points)


def _steps(edges, values):
    """Step-line coordinates of a histogram; NaN bins become gaps."""
    return dict(x=np.repeat(edges, 2)[1:-1], y=np.repeat(values, 2))


def _closed(points):
    return np.vstack([points, points[:1]])


def tube_mesh(state, r_values, tube_radius, n_alpha=12):
    """Cosmetic torus around the core curve, coloured by the core ring value.

    The cross-section uses the local radial direction from the ring centre,
    which is well defined for every star-shaped ring family of the toy.
    """
    x, tangent = state.x, state.tangent
    rel = x - state.center_deformed
    radial = rel - np.einsum("ij,ij->i", rel, tangent)[:, None] * tangent
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    binormal = np.cross(tangent, radial)
    alpha = 2.0 * np.pi * np.arange(n_alpha) / n_alpha
    n = len(x)
    vertices = (x[:, None, :] + tube_radius * (np.cos(alpha)[None, :, None] * radial[:, None, :]
                                               + np.sin(alpha)[None, :, None] * binormal[:, None, :])).reshape(-1, 3)
    k, j = np.meshgrid(np.arange(n), np.arange(n_alpha), indexing="ij")
    a = (k * n_alpha + j).ravel()
    b = (((k + 1) % n) * n_alpha + j).ravel()
    c = (k * n_alpha + (j + 1) % n_alpha).ravel()
    d = (((k + 1) % n) * n_alpha + (j + 1) % n_alpha).ravel()
    return dict(**_xyz(vertices), i=np.concatenate([a, b]), j=np.concatenate([b, d]),
                k=np.concatenate([c, c]), intensity=np.repeat(np.nan_to_num(r_values, nan=0.0), n_alpha))


def _empty3d():
    return dict(x=[], y=[], z=[])


def colour_scale(view):
    """Symmetric colour range of r, following the largest |P| in the view."""
    magnitude = max(float(np.max(np.linalg.norm(res.state.P, axis=1))) for res in view.results.values())
    return magnitude if magnitude > 0.0 else 1.0


def scene_extent(view):
    """Half-width of a position scene that holds vertex, rings and almond."""
    reach = [float(np.max(np.abs(res.state.x))) for res in view.results.values()]
    reach += [float(np.max(np.abs(points))) for _, points in view.almond_contours]
    return 1.15 * max(reach + [1.0])


def trace_data(view, options):
    """Data of every trace, in a fixed order: [(kind, row, col, properties), ...]."""
    config = view.config
    frame = next(iter(view.results.values())).state.frame
    scale = colour_scale(view)
    extent = options.scene_extent if options.scene_extent is not None else scene_extent(view)
    traces = []

    # ---- position space: vertex, jet axis, event plane, almond ----
    traces.append(("scatter3d", 1, 1, dict(x=[0.0], y=[0.0], z=[0.0], mode="markers", name="vertex",
                                           marker=dict(size=4, color="#444444", symbol="diamond"))))
    jet_tip = 0.9 * extent * frame.t_hat
    traces.append(("scatter3d", 1, 1, dict(x=[0.0, jet_tip[0]], y=[0.0, jet_tip[1]], z=[0.0, jet_tip[2]],
                                           mode="lines", name="jet axis",
                                           line=dict(color="#aa6600", width=7))))
    plane = 0.9 * extent * np.array([np.cos(config.flow.psi2), np.sin(config.flow.psi2), 0.0])
    traces.append(("scatter3d", 1, 1, dict(x=[-plane[0], plane[0]], y=[-plane[1], plane[1]], z=[0.0, 0.0],
                                           mode="lines", name="event plane",
                                           line=dict(color="#888888", width=3, dash="dash"))))
    contours = dict(view.almond_contours)
    if 1.0 in contours:
        inner = contours[1.0]
        fan = np.vstack([np.zeros(3), inner])
        n = len(inner)
        traces.append(("mesh3d", 1, 1, dict(**_xyz(fan), i=np.zeros(n - 1, dtype=int),
                                            j=np.arange(1, n), k=np.arange(2, n + 1),
                                            color=ALMOND_FILL, opacity=0.35, name="almond (1 sigma)",
                                            showlegend=True, hoverinfo="skip")))
        outer = contours[2.0]
        traces.append(("scatter3d", 1, 1, dict(**_xyz(outer), mode="lines", name="almond (2 sigma)",
                                               line=dict(color=ALMOND_FILL, width=3))))
    else:
        traces.append(("mesh3d", 1, 1, dict(**_empty3d(), i=[], j=[], k=[], name="almond (1 sigma)")))
        traces.append(("scatter3d", 1, 1, dict(**_empty3d(), mode="lines", name="almond (2 sigma)")))

    # ---- position space: one fixed block of traces per regime ----
    for regime in REGIMES:
        shown = regime in view.results and regime in options.regimes_shown
        color = REGIME_COLOR[regime]
        if not shown:
            traces.append(("scatter3d", 1, 1, dict(**_empty3d(), mode="lines", name=regime,
                                                   legendgroup=regime, showlegend=False)))
            traces.append(("scatter3d", 1, 1, dict(**_empty3d(), mode="lines", legendgroup=regime,
                                                   showlegend=False)))
            traces.append(("scatter3d", 1, 1, dict(**_empty3d(), mode="markers", legendgroup=regime,
                                                   showlegend=False)))
            traces.append(("mesh3d", 1, 1, dict(**_empty3d(), i=[], j=[], k=[], legendgroup=regime)))
            traces.append(("scatter3d", 1, 1, dict(**_empty3d(), mode="markers", legendgroup=regime,
                                                   showlegend=False)))
            continue
        result = view.results[regime]
        state, r = result.state, result.measurement.r
        size = float(np.max(np.linalg.norm(state.x - state.center_deformed, axis=1)))
        arrow_scale = 0.4 * size / scale
        step = max(1, len(state.x) // options.max_arrows)
        starts = state.x[::step]
        arrows = arrow_scale * state.P[::step]
        traces.append(("scatter3d", 1, 1, dict(**_xyz(_closed(state.x)), mode="lines", name=regime,
                                               legendgroup=regime, showlegend=True,
                                               line=dict(color=color, width=5))))
        traces.append(("scatter3d", 1, 1, dict(**_segments(starts, arrows), mode="lines",
                                               legendgroup=regime, showlegend=False, hoverinfo="skip",
                                               line=dict(color=color, width=3))))
        traces.append(("scatter3d", 1, 1, dict(**_xyz(starts + arrows), mode="markers",
                                               legendgroup=regime, showlegend=False,
                                               marker=dict(size=3, color=np.nan_to_num(r[::step], nan=0.0),
                                                           colorscale=RING_COLORSCALE, cmin=-scale, cmax=scale))))
        if options.tube_radius > 0.0:
            traces.append(("mesh3d", 1, 1, dict(**tube_mesh(state, r, options.tube_radius),
                                                colorscale=RING_COLORSCALE, cmin=-scale, cmax=scale,
                                                opacity=0.5, showscale=False, legendgroup=regime,
                                                hoverinfo="skip")))
        else:
            traces.append(("mesh3d", 1, 1, dict(**_empty3d(), i=[], j=[], k=[], legendgroup=regime)))
        mask = highlight_mask(result, options.highlight_bin)
        traces.append(("scatter3d", 1, 1, dict(**_xyz(state.x[mask]), mode="markers", legendgroup=regime,
                                               showlegend=False,
                                               marker=dict(size=5, color="#ffcc00",
                                                           line=dict(color="#000000", width=1)))))

    # ---- momentum sphere ----
    theta, phi = np.meshgrid(np.linspace(0.0, np.pi, 25), np.linspace(0.0, 2.0 * np.pi, 49))
    traces.append(("surface", 1, 2, dict(x=np.sin(theta) * np.cos(phi), y=np.sin(theta) * np.sin(phi),
                                         z=np.cos(theta), opacity=0.12, showscale=False, hoverinfo="skip",
                                         colorscale=[[0.0, "#bbbbbb"], [1.0, "#bbbbbb"]])))
    tip = 1.3 * frame.t_hat
    traces.append(("scatter3d", 1, 2, dict(x=[0.0, tip[0]], y=[0.0, tip[1]], z=[0.0, tip[2]], mode="lines",
                                           showlegend=False, line=dict(color="#aa6600", width=7))))
    circle = np.linspace(0.0, 2.0 * np.pi, 97)
    eta_max = config.acceptance.eta_max
    for sign in (-1.0, 1.0):
        if eta_max is None:
            traces.append(("scatter3d", 1, 2, dict(**_empty3d(), mode="lines", showlegend=False)))
        else:
            traces.append(("scatter3d", 1, 2, dict(x=np.cos(circle) / np.cosh(eta_max),
                                                   y=np.sin(circle) / np.cosh(eta_max),
                                                   z=np.full(circle.shape, sign * np.tanh(eta_max)),
                                                   mode="lines", showlegend=False, hoverinfo="skip",
                                                   line=dict(color="#888888", width=3, dash="dash"))))
    for regime in REGIMES:
        if regime in view.results and regime in options.regimes_shown:
            result = view.results[regime]
            m, p_hat = result.measurement, result.state.p_hat
            finite = np.isfinite(m.r)
            traces.append(("scatter3d", 1, 2, dict(**_xyz(p_hat[finite]), mode="markers", legendgroup=regime,
                                                   showlegend=False,
                                                   marker=dict(size=2.5, color=m.r[finite],
                                                               colorscale=RING_COLORSCALE,
                                                               cmin=-scale, cmax=scale))))
            mask = highlight_mask(result, options.highlight_bin)
            traces.append(("scatter3d", 1, 2, dict(**_xyz(p_hat[mask]), mode="markers", legendgroup=regime,
                                                   showlegend=False,
                                                   marker=dict(size=4.5, color="#ffcc00",
                                                               line=dict(color="#000000", width=1)))))
        else:
            traces.append(("scatter3d", 1, 2, dict(**_empty3d(), mode="markers", showlegend=False)))
            traces.append(("scatter3d", 1, 2, dict(**_empty3d(), mode="markers", showlegend=False)))

    # ---- side plots: Delta phi, Delta theta, sweep ----
    for row, variable in ((1, "delta_phi"), (2, "delta_theta")):
        for regime in REGIMES:
            if regime in view.results and regime in options.regimes_shown:
                binned = getattr(view.results[regime], variable)
                values = binned.signal_density() if options.quantity == "signal_density" else binned.ring_average()
                data = _steps(binned.edges, values)
            else:
                data = dict(x=[], y=[])
            traces.append(("scatter", row, 3, dict(**data, mode="lines", legendgroup=regime, showlegend=False,
                                                   line=dict(color=REGIME_COLOR[regime], dash=REGIME_DASH[regime],
                                                             width=2))))
    for regime in REGIMES:
        if view.sweep is not None and regime in view.sweep_summaries and regime in options.regimes_shown:
            y = [s.ring_average if options.quantity == "ring_average" else s.signal_per_total_weight
                 for s in view.sweep_summaries[regime]]
            data = dict(x=np.asarray(view.sweep.values, dtype=float), y=np.asarray(y, dtype=float))
        else:
            data = dict(x=[], y=[])
        traces.append(("scatter", 3, 3, dict(**data, mode="lines", legendgroup=regime, showlegend=False,
                                             line=dict(color=REGIME_COLOR[regime], dash=REGIME_DASH[regime],
                                                       width=2))))
    return traces


def layout_updates(view, options):
    """Layout pieces that follow the View: title readouts, markers, scene ranges."""
    shapes = []
    if options.highlight_bin is not None:
        edges = next(iter(view.results.values())).delta_phi.edges
        shapes.append(dict(type="rect", xref="x", yref="y domain", x0=edges[options.highlight_bin],
                           x1=edges[options.highlight_bin + 1], y0=0.0, y1=1.0,
                           fillcolor="#ffcc00", opacity=0.3, line=dict(width=0)))
    if view.sweep is not None:
        shapes.append(dict(type="line", xref="x3", yref="y3 domain", x0=view.sweep_value,
                           x1=view.sweep_value, y0=0.0, y1=1.0, line=dict(color="#aa6600", width=2)))
    parts = []
    for regime, row in readouts(view).items():
        parts.append(f"{regime}: R={row['ring_average']:.3f} acc={row['acceptance_fraction']:.2f}")
    title = " | ".join(parts)
    if view.almond_contours:
        title += f" | almond eps2={view.almond_eccentricity:+.3f}"
    extent = options.scene_extent if options.scene_extent is not None else scene_extent(view)
    axis = dict(range=[-extent, extent], autorange=False)
    return dict(shapes=shapes, title=dict(text=title, font=dict(size=12)),
                scene=dict(xaxis=axis, yaxis=axis, zaxis=axis, aspectmode="cube"))


def make_figure(view, options=FigureOptions(), widget=False):
    """Build a figure (or a FigureWidget) with the traces of a View."""
    quantity_label = "signal density" if options.quantity == "signal_density" else "ring average per bin"
    sweep_label = f"integrated vs {view.sweep.parameter}" if view.sweep is not None else "integrated (no sweep)"
    fig = make_subplots(rows=3, cols=3, column_widths=[0.38, 0.30, 0.32],
                        specs=[[{"type": "scene", "rowspan": 3}, {"type": "scene", "rowspan": 3}, {"type": "xy"}],
                               [None, None, {"type": "xy"}],
                               [None, None, {"type": "xy"}]],
                        subplot_titles=("position space", "momentum directions",
                                        f"{quantity_label} vs Delta phi", f"{quantity_label} vs Delta theta",
                                        sweep_label),
                        horizontal_spacing=0.03, vertical_spacing=0.09)
    if widget:
        fig = go.FigureWidget(fig)
    for kind, row, col, properties in trace_data(view, options):
        fig.add_trace(TRACE_CLASSES[kind](**properties), row=row, col=col)
    fig.update_xaxes(title_text=LABEL_DPHI, range=[-np.pi, np.pi], row=1, col=3)
    fig.update_xaxes(title_text=LABEL_DTHETA, range=[0.0, np.pi], row=2, col=3)
    if view.sweep is not None:
        fig.update_xaxes(title_text=view.sweep.parameter, row=3, col=3)
    fig.update_scenes(aspectmode="cube", row=1, col=2)
    fig.update_layout(height=options.height, margin=dict(l=10, r=10, t=70, b=40),
                      legend=dict(orientation="h", y=-0.02), **layout_updates(view, options))
    return fig


def update_figure(fig, view, options):
    """Refresh an existing figure in place (same trace order as make_figure)."""
    data = trace_data(view, options)
    if len(data) != len(fig.data):
        raise RuntimeError("trace layout changed; rebuild the figure with make_figure")
    with fig.batch_update():
        for trace, (kind, _, _, properties) in zip(fig.data, data):
            trace.update(**properties)
        fig.update_layout(**layout_updates(view, options))


def animated_figure(config, parameter, values, options=FigureOptions(), sweep=None, regimes=REGIMES,
                    frame_ms=120, n_points=240):
    """A movie of the scene while one parameter runs through values.

    Every frame stores all trace data in the HTML file, so the ring is drawn
    with at most n_points elements (integrated results are unaffected at the
    precision shown); this keeps movies at a few MB instead of tens.
    """
    config = dataclasses.replace(config, shape=dataclasses.replace(
        config.shape, n_points=min(config.shape.n_points, n_points)))
    views = [compute_view(with_parameter(config, parameter, value), regimes, sweep) for value in values]
    extent = options.scene_extent if options.scene_extent is not None else max(scene_extent(v) for v in views)
    options = dataclasses.replace(options, scene_extent=extent)
    fig = make_figure(views[0], options)
    frames = []
    for index, (value, view) in enumerate(zip(values, views)):
        data = [TRACE_CLASSES[kind](**properties) for kind, _, _, properties in trace_data(view, options)]
        frames.append(go.Frame(name=str(index), data=data, layout=layout_updates(view, options)))
    fig.frames = frames
    play = dict(label="play", method="animate",
                args=[None, dict(frame=dict(duration=frame_ms, redraw=True), fromcurrent=True,
                                 transition=dict(duration=0))])
    pause = dict(label="pause", method="animate",
                 args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])
    steps = [dict(method="animate", label=f"{float(value):.2f}",
                  args=[[str(index)], dict(frame=dict(duration=0, redraw=True), mode="immediate")])
             for index, value in enumerate(values)]
    fig.update_layout(updatemenus=[dict(type="buttons", buttons=[play, pause], x=0.0, y=1.08, direction="left")],
                      sliders=[dict(steps=steps, currentvalue=dict(prefix=f"{parameter} = "), y=-0.06)])
    return fig
