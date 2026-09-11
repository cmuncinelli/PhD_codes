"""Tests of the view layer and of the figure bookkeeping shared by widgets and movies."""

import dataclasses

import numpy as np
import pytest

from ringtoy.scans import configure_parameter, get_parameter
from ringtoy.scenarios import apply_override, get_preset
from ringtoy.views import SweepSpec, compute_view, sweep_summaries

plotly = pytest.importorskip("plotly")
from ringtoy.viz3d import FigureOptions, animated_figure, make_figure, trace_data, update_figure  # noqa: E402


@pytest.mark.parametrize("parameter, value", [("delta", 0.41), ("jet_eta", -0.33),
                                              ("flow.strength", 1.7), ("shape.radius", 1.25)])
def test_parameter_round_trip(parameter, value):
    config = configure_parameter(get_preset("regimes"), parameter, value)
    np.testing.assert_allclose(get_parameter(config, parameter), value, rtol=1e-14)


def test_sweep_does_not_depend_on_current_value_of_swept_parameter():
    base = get_preset("regimes")
    sweep = SweepSpec(parameter="jet_phi", values=tuple(np.linspace(-1.0, 1.0, 5)))
    first = sweep_summaries(configure_parameter(base, "jet_phi", 0.3), sweep)
    second = sweep_summaries(configure_parameter(base, "jet_phi", -2.2), sweep)
    for regime in first:
        assert [s.ring_average for s in first[regime]] == [s.ring_average for s in second[regime]]


def test_trace_layout_is_invariant():
    # The notebook updates traces by position, so every configuration and
    # display option must produce the same number, kind and subplot of traces.
    base = get_preset("regimes")
    sweep = SweepSpec(parameter="flow.strength", values=(0.0, 0.5, 1.0))
    variants = [
        (base, FigureOptions()),
        (apply_override(base, "shape.kind", "fourier"), FigureOptions(tube_radius=0.1, highlight_bin=20)),
        (apply_override(base, "almond.show", False), FigureOptions(regimes_shown=("R2",))),
        (apply_override(base, "acceptance.eta_max", 0.5), FigureOptions(quantity="ring_average")),
    ]
    reference = [(kind, row, col) for kind, row, col, _ in trace_data(compute_view(base, sweep=sweep), FigureOptions())]
    for config, options in variants:
        layout = [(kind, row, col) for kind, row, col, _ in trace_data(compute_view(config, sweep=sweep), options)]
        assert layout == reference


def test_update_figure_matches_fresh_figure():
    base = get_preset("regimes")
    options = FigureOptions(tube_radius=0.1, highlight_bin=33)
    changed = dataclasses.replace(base, jet_eta=0.44, flow=dataclasses.replace(base.flow, v2=-0.12))
    figure = make_figure(compute_view(base), options)
    update_figure(figure, compute_view(changed), options)
    fresh = make_figure(compute_view(changed), options)
    for updated_trace, fresh_trace in zip(figure.data, fresh.data):
        for axis in ("x", "y"):
            a, b = getattr(updated_trace, axis), getattr(fresh_trace, axis)
            if a is not None and len(a):
                np.testing.assert_allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float), equal_nan=True)


def test_animation_has_one_frame_per_value_with_all_traces():
    base = apply_override(get_preset("regimes"), "shape.n_points", 90)
    values = (0.0, 0.8, 1.6)
    figure = animated_figure(base, "flow.strength", values, FigureOptions())
    assert len(figure.frames) == len(values)
    assert all(len(frame.data) == len(figure.data) for frame in figure.frames)
