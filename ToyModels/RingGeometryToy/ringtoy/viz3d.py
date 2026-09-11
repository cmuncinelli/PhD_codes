"""
viz3d.py -- Interactive plotly scenes of the ring geometry toy.

Part of the ring geometry toy (README.md describes the scenes). Scenes are
written as self-contained HTML files, which works on headless nodes.
"""

import dataclasses

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

REGIME_COLOR = {"R0": "#000000", "R1": "#0033cc", "R2": "#cc0000", "R3": "#008800"}
RING_COLORSCALE = "RdBu_r"


def _segments(starts, vectors):
    """Line segments start --> start + vector, separated by None for one trace."""
    n = len(starts)
    points = np.full((3 * n, 3), np.nan)
    points[0::3] = starts
    points[1::3] = starts + vectors
    # plotly breaks lines at None; NaN rows are converted explicitly.
    coords = [[None if np.isnan(value) else float(value) for value in points[:, axis]]
              for axis in range(3)]
    return coords


def _closed(array):
    return np.vstack([array, array[:1]])


def tube_mesh(state, r_values, tube_radius, n_alpha=16):
    """Cosmetic torus around the core curve, coloured by the core ring value.

    The tube does not enter any observable. Its cross-section uses the local
    radial direction from the ring centre, which is well defined for every
    star-shaped ring family of the toy.
    """
    x, tangent = state.x, state.tangent
    rel = x - state.center_deformed
    radial = rel - np.einsum("ij,ij->i", rel, tangent)[:, None] * tangent
    radial /= np.linalg.norm(radial, axis=1)[:, None]
    binormal = np.cross(tangent, radial)
    alpha = 2.0 * np.pi * np.arange(n_alpha) / n_alpha
    n = len(x)
    ring = (x[:, None, :]
            + tube_radius * (np.cos(alpha)[None, :, None] * radial[:, None, :]
                             + np.sin(alpha)[None, :, None] * binormal[:, None, :]))
    vertices = ring.reshape(-1, 3)
    k, j = np.meshgrid(np.arange(n), np.arange(n_alpha), indexing="ij")
    a = (k * n_alpha + j).ravel()
    b = (((k + 1) % n) * n_alpha + j).ravel()
    c = (k * n_alpha + (j + 1) % n_alpha).ravel()
    d = (((k + 1) % n) * n_alpha + (j + 1) % n_alpha).ravel()
    return dict(x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
                i=np.concatenate([a, b]), j=np.concatenate([b, d]), k=np.concatenate([c, c]),
                intensity=np.repeat(np.nan_to_num(r_values, nan=0.0), n_alpha))


def _position_traces(regime, result, arrow_scale, max_arrows, scale, tube_radius, show_legend):
    state, r = result.state, result.measurement.r
    color = REGIME_COLOR[regime]
    traces = []
    closed = _closed(state.x)
    traces.append(go.Scatter3d(x=closed[:, 0], y=closed[:, 1], z=closed[:, 2], mode="lines",
                               line=dict(color=color, width=4), name=regime,
                               legendgroup=regime, showlegend=show_legend))
    step = max(1, len(state.x) // max_arrows)
    starts = state.x[::step]
    tips = starts + arrow_scale * state.P[::step]
    sx, sy, sz = _segments(starts, arrow_scale * state.P[::step])
    traces.append(go.Scatter3d(x=sx, y=sy, z=sz, mode="lines", line=dict(color=color, width=3),
                               legendgroup=regime, showlegend=False, hoverinfo="skip"))
    traces.append(go.Scatter3d(x=tips[:, 0], y=tips[:, 1], z=tips[:, 2], mode="markers",
                               marker=dict(size=3, color=np.nan_to_num(r[::step], nan=0.0),
                                           colorscale=RING_COLORSCALE, cmin=-scale, cmax=scale),
                               legendgroup=regime, showlegend=False,
                               hovertemplate="r = %{marker.color:.3f}<extra>" + regime + "</extra>"))
    if tube_radius:
        traces.append(go.Mesh3d(**tube_mesh(state, r, tube_radius), colorscale=RING_COLORSCALE,
                                cmin=-scale, cmax=scale, opacity=0.55, showscale=False,
                                legendgroup=regime, showlegend=False, hoverinfo="skip"))
    return traces


def _momentum_traces(regime, result, scale):
    m = result.measurement
    p_hat = result.state.p_hat
    finite = np.isfinite(m.r)
    return [go.Scatter3d(x=p_hat[finite, 0], y=p_hat[finite, 1], z=p_hat[finite, 2], mode="markers",
                         marker=dict(size=2.5, color=m.r[finite], colorscale=RING_COLORSCALE,
                                     cmin=-scale, cmax=scale),
                         legendgroup=regime, showlegend=False,
                         hovertemplate="r = %{marker.color:.3f}<extra>" + regime + "</extra>")]


def _static_traces(config, frame, length):
    traces = [
        go.Scatter3d(x=[0.0], y=[0.0], z=[0.0], mode="markers",
                     marker=dict(size=5, color="#444444", symbol="diamond"), name="vertex"),
        go.Scatter3d(x=[0.0, length * frame.t_hat[0]], y=[0.0, length * frame.t_hat[1]],
                     z=[0.0, length * frame.t_hat[2]], mode="lines",
                     line=dict(color="#aa6600", width=6), name="jet axis"),
    ]
    psi2 = config.flow.psi2
    traces.append(go.Scatter3d(x=[-length * np.cos(psi2), length * np.cos(psi2)],
                               y=[-length * np.sin(psi2), length * np.sin(psi2)], z=[0.0, 0.0],
                               mode="lines", line=dict(color="#888888", width=2, dash="dash"),
                               name="event plane"))
    return traces


def _sphere_traces(frame, eta_max):
    theta, phi = np.meshgrid(np.linspace(0.0, np.pi, 30), np.linspace(0.0, 2.0 * np.pi, 60))
    traces = [go.Surface(x=np.sin(theta) * np.cos(phi), y=np.sin(theta) * np.sin(phi),
                         z=np.cos(theta), opacity=0.12, showscale=False,
                         colorscale=[[0.0, "#bbbbbb"], [1.0, "#bbbbbb"]], hoverinfo="skip"),
              go.Scatter3d(x=[0.0, 1.3 * frame.t_hat[0]], y=[0.0, 1.3 * frame.t_hat[1]],
                           z=[0.0, 1.3 * frame.t_hat[2]], mode="lines",
                           line=dict(color="#aa6600", width=6), showlegend=False)]
    if eta_max is not None:
        circle = np.linspace(0.0, 2.0 * np.pi, 120)
        for sign in (-1.0, 1.0):
            traces.append(go.Scatter3d(x=np.cos(circle) / np.cosh(eta_max),
                                       y=np.sin(circle) / np.cosh(eta_max),
                                       z=np.full(circle.shape, sign * np.tanh(eta_max)),
                                       mode="lines", line=dict(color="#888888", width=2, dash="dash"),
                                       showlegend=False, hoverinfo="skip"))
    return traces


def _geometry_scales(results):
    states = [res.state for res in results.values()]
    size = max(float(np.max(np.linalg.norm(s.x - s.center_deformed, axis=1))) for s in states)
    reach = max(float(np.max(np.linalg.norm(s.x, axis=1))) for s in states)
    magnitude = max(float(np.max(np.linalg.norm(s.P, axis=1))) for s in states)
    magnitude = magnitude if magnitude > 0.0 else 1.0
    return size, reach, magnitude


def scene_figure(config, results, tube_radius=None, max_arrows=48):
    """Position space (left) and momentum sphere (right) for several regimes."""
    first = next(iter(results.values()))
    frame = first.state.frame
    size, reach, magnitude = _geometry_scales(results)
    arrow_scale = 0.4 * size / magnitude
    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=("position space", "momentum directions"))
    for trace in _static_traces(config, frame, 1.2 * reach):
        fig.add_trace(trace, row=1, col=1)
    for trace in _sphere_traces(frame, config.acceptance.eta_max):
        fig.add_trace(trace, row=1, col=2)
    for regime, result in results.items():
        for trace in _position_traces(regime, result, arrow_scale, max_arrows, magnitude,
                                      tube_radius, show_legend=True):
            fig.add_trace(trace, row=1, col=1)
        for trace in _momentum_traces(regime, result, magnitude):
            fig.add_trace(trace, row=1, col=2)
    fig.update_scenes(aspectmode="data", row=1, col=1)
    fig.update_scenes(aspectmode="cube", row=1, col=2)
    fig.update_layout(height=720, legend=dict(itemsizing="constant"))
    return fig


def strength_slider_figure(config, strengths, evaluate, regime="R3", tube_radius=None):
    """One regime with a slider over the shape-map strength.

    evaluate is scans.evaluate, passed in to keep this module free of the
    scan machinery.
    """
    steps_results = [evaluate(dataclasses.replace(
        config, flow=dataclasses.replace(config.flow, strength=float(s))), regime) for s in strengths]
    all_results = {f"{regime}_{k}": res for k, res in enumerate(steps_results)}
    frame = steps_results[0].state.frame
    size, reach, magnitude = _geometry_scales(all_results)
    arrow_scale = 0.4 * size / magnitude

    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
                        subplot_titles=("position space", "momentum directions"))
    static = _static_traces(config, frame, 1.2 * reach) + _sphere_traces(frame, config.acceptance.eta_max)
    n_static_left = len(_static_traces(config, frame, 1.2 * reach))
    for index, trace in enumerate(static):
        fig.add_trace(trace, row=1, col=1 if index < n_static_left else 2)
    n_static = len(static)

    step_trace_ranges = []
    for k, result in enumerate(steps_results):
        start = len(fig.data)
        for trace in _position_traces(regime, result, arrow_scale, 48, magnitude, tube_radius,
                                      show_legend=False):
            trace.visible = (k == 0)
            fig.add_trace(trace, row=1, col=1)
        for trace in _momentum_traces(regime, result, magnitude):
            trace.visible = (k == 0)
            fig.add_trace(trace, row=1, col=2)
        step_trace_ranges.append((start, len(fig.data)))

    slider_steps = []
    for k, (start, stop) in enumerate(step_trace_ranges):
        visible = [True] * n_static + [False] * (len(fig.data) - n_static)
        for index in range(start, stop):
            visible[index] = True
        slider_steps.append(dict(method="update", args=[{"visible": visible}],
                                 label=f"{strengths[k]:.2f}"))
    fig.update_layout(sliders=[dict(active=0, steps=slider_steps,
                                    currentvalue=dict(prefix="shape-map strength s = "))],
                      height=760, title=f"{regime} versus shape-map strength")
    fig.update_scenes(aspectmode="data", row=1, col=1)
    fig.update_scenes(aspectmode="cube", row=1, col=2)
    return fig
