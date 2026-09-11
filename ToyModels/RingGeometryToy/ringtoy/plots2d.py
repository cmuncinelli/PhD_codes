"""
plots2d.py -- Matplotlib figures of the ring geometry toy.

Part of the ring geometry toy (README.md explains what each figure shows).
Every figure is written as both PDF and PNG.
"""

import matplotlib

# Headless backend: the toy runs over ssh on the compute nodes.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

from .observable import summary_from_sums

REGIME_STYLE = {
    "R0": dict(color="#000000", linestyle="-"),
    "R1": dict(color="#0033cc", linestyle="--"),
    "R2": dict(color="#cc0000", linestyle="-."),
    "R3": dict(color="#008800", linestyle=":"),
}

REGIME_TITLE = {
    "R0": "R0: rigid ring, radial momenta",
    "R1": "R1: deformed ring, radial momenta",
    "R2": "R2: rigid ring, deflected momenta",
    "R3": "R3: deformed ring, deflected momenta",
}

LABEL_DPHI = r"$\Delta\phi = \phi_{\Lambda} - \phi_{\mathrm{Jet}}$ [rad]"
LABEL_DTHETA = r"$\Delta\theta_{\mathrm{Jet},\Lambda}$ [rad]"
LABEL_RING = r"$\langle \vec{P} \cdot \hat{n} \rangle$"
LABEL_DENSITY = r"$\sum w\,r \,/\, (W_{\mathrm{tot}}\,\Delta x)$"


def apply_style():
    """Sober HEP-like style: ticks inside on all sides, minor ticks on."""
    plt.rcParams.update({
        "font.size": 11,
        "axes.linewidth": 1.0,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "legend.frameon": False,
        "savefig.bbox": "tight",
    })


def save_figure(fig, path_stem):
    """Write a figure as PDF and PNG and close it."""
    fig.savefig(f"{path_stem}.pdf")
    fig.savefig(f"{path_stem}.png", dpi=150)
    plt.close(fig)


def _magnitude_scale(results):
    # Colour and axis scales follow the polarization magnitude, so that
    # f < 1 does not leave the plots mostly empty.
    magnitudes = [np.nanmax(np.linalg.norm(res.state.P, axis=1)) for res in results.values()]
    scale = max(magnitudes) if magnitudes else 1.0
    return scale if scale > 0.0 else 1.0


def plot_binned_overlay(curves, variable, path_stem, title=""):
    """Ring average and signal density versus one binned variable.

    curves   : list of (label, BinnedRing, style dict).
    variable : "delta_phi" or "delta_theta".
    """
    xlabel = LABEL_DPHI if variable == "delta_phi" else LABEL_DTHETA
    fig, (ax_avg, ax_den) = plt.subplots(2, 1, figsize=(6.4, 7.2), sharex=True)
    for label, binned, style in curves:
        # Step plots: the histograms are bin sums, not samples of a function.
        ax_avg.stairs(binned.ring_average(), binned.edges, label=label, **style)
        ax_den.stairs(binned.signal_density(), binned.edges, label=label, **style)
    ax_avg.set_ylabel(LABEL_RING + " per bin")
    ax_den.set_ylabel(LABEL_DENSITY)
    ax_den.set_xlabel(xlabel)
    ax_avg.axhline(0.0, color="#888888", linewidth=0.6)
    ax_den.axhline(0.0, color="#888888", linewidth=0.6)
    ax_avg.legend(fontsize=9)
    if title:
        ax_avg.set_title(title, fontsize=10)
    save_figure(fig, path_stem)


def plot_regime_overlays(results, path_prefix):
    """Delta phi and Delta theta overlays of all evaluated regimes."""
    for variable in ("delta_phi", "delta_theta"):
        curves = [(REGIME_TITLE[regime], getattr(res, variable), REGIME_STYLE[regime])
                  for regime, res in results.items()]
        plot_binned_overlay(curves, variable, f"{path_prefix}_{variable}")


def plot_cos_gamma(results, path_stem):
    """Distribution of the angle between P and n_hat for accepted elements."""
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    edges = np.linspace(-1.0, 1.0, 51)
    for regime, res in results.items():
        m = res.measurement
        sel = m.accepted & np.isfinite(m.cos_gamma)
        counts, _ = np.histogram(m.cos_gamma[sel], bins=edges, weights=m.weight[sel])
        total = np.sum(counts)
        if total > 0.0:
            counts = counts / (total * np.diff(edges))
        ax.stairs(counts, edges, label=REGIME_TITLE[regime], **REGIME_STYLE[regime])
    ax.set_xlabel(r"$\cos\gamma = \hat{P} \cdot \hat{n}$")
    ax.set_ylabel("weighted density")
    ax.set_yscale("log")
    ax.legend(fontsize=9, loc="upper left")
    save_figure(fig, path_stem)


def plot_jet_plane(results, path_stem, max_arrows=48):
    """Polarization vectors projected on the plane transverse to the jet.

    Positions are measured from the (possibly deformed) ring centre in the jet
    frame (e1, e2); arrows are (P . e1, P . e2) and are coloured by r = P . n_hat.
    The out-of-plane component P . t_hat is quoted in each panel.
    """
    scale = _magnitude_scale(results)
    norm = Normalize(vmin=-scale, vmax=scale)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.6), squeeze=False)
    quiver = None
    for ax, (regime, res) in zip(axes.flat, results.items()):
        state = res.state
        frame = state.frame
        rel = state.x - state.center_deformed
        u, v = rel @ frame.e1, rel @ frame.e2
        pu, pv, pt = state.P @ frame.e1, state.P @ frame.e2, state.P @ frame.t_hat
        colour = np.nan_to_num(res.measurement.r, nan=0.0)
        step = max(1, len(u) // max_arrows)
        closed = np.r_[np.arange(len(u)), 0]
        ax.plot(u[closed], v[closed], color="#999999", linewidth=0.8)
        extent = np.max(np.hypot(u, v))
        quiver = ax.quiver(u[::step], v[::step], pu[::step], pv[::step], colour[::step],
                           cmap="RdBu_r", norm=norm, angles="xy", scale_units="xy",
                           scale=scale / (0.35 * extent), width=0.006)
        ax.plot(0.0, 0.0, marker="+", color="#000000")
        ax.set_aspect("equal")
        ax.set_xlim(-1.45 * extent, 1.45 * extent)
        ax.set_ylim(-1.45 * extent, 1.45 * extent)
        ax.set_title(REGIME_TITLE[regime], fontsize=10)
        ax.set_xlabel(r"$e_1$ [fm]")
        ax.set_ylabel(r"$e_2$ [fm]")
        weights = state.weight
        mean_abs_pt = np.sum(weights * np.abs(pt)) / np.sum(weights)
        ax.text(0.03, 0.03, r"$\langle|\vec{P}\cdot\hat{t}|\rangle$ = " + f"{mean_abs_pt:.3f}",
                transform=ax.transAxes, fontsize=9)
    for ax in axes.flat[len(results):]:
        ax.set_visible(False)
    if quiver is not None:
        fig.colorbar(quiver, ax=axes, shrink=0.8, label=r"$r = \vec{P} \cdot \hat{n}$")
    save_figure(fig, path_stem)


def plot_momentum_map(results, eta_max, path_stem):
    """Ring elements in (Delta phi, eta) of their momentum, coloured by r."""
    scale = _magnitude_scale(results)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.6), sharex=True, sharey=True, squeeze=False)
    scatter = None
    for ax, (regime, res) in zip(axes.flat, results.items()):
        m = res.measurement
        finite = np.isfinite(m.r)
        scatter = ax.scatter(m.delta_phi[finite], m.eta[finite], c=m.r[finite], s=4,
                             cmap="RdBu_r", vmin=-scale, vmax=scale)
        if eta_max is not None:
            for sign in (-1.0, 1.0):
                ax.axhline(sign * eta_max, color="#888888", linewidth=0.8, linestyle="--")
        ax.set_title(REGIME_TITLE[regime], fontsize=10)
        ax.set_xlim(-np.pi, np.pi)
    for ax in axes[-1, :]:
        ax.set_xlabel(LABEL_DPHI)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\eta$ of $\hat{p}$")
    for ax in axes.flat[len(results):]:
        ax.set_visible(False)
    if scatter is not None:
        fig.colorbar(scatter, ax=axes, shrink=0.8, label=r"$r = \vec{P} \cdot \hat{n}$")
    save_figure(fig, path_stem)


def plot_scan(values, sums_per_regime, xlabel, path_stem, title=""):
    """Integrated ring average, diluted signal, and acceptance versus a parameter."""
    fig, axes = plt.subplots(3, 1, figsize=(6.4, 8.6), sharex=True)
    for regime, sums_list in sums_per_regime.items():
        summaries = [summary_from_sums(s) for s in sums_list]
        style = REGIME_STYLE[regime]
        axes[0].plot(values, [s.ring_average for s in summaries], marker="o", markersize=3,
                     label=REGIME_TITLE[regime], **style)
        axes[1].plot(values, [s.signal_per_total_weight for s in summaries], marker="o",
                     markersize=3, **style)
        axes[2].plot(values, [s.acceptance_fraction for s in summaries], marker="o",
                     markersize=3, **style)
    axes[0].set_ylabel(LABEL_RING)
    axes[1].set_ylabel(r"$\sum w\,r \,/\, W_{\mathrm{tot}}$")
    axes[2].set_ylabel("acceptance fraction")
    axes[2].set_xlabel(xlabel)
    for ax in axes:
        ax.axhline(0.0, color="#888888", linewidth=0.6)
    axes[0].legend(fontsize=9)
    if title:
        axes[0].set_title(title, fontsize=10)
    save_figure(fig, path_stem)


def plot_ellipticity_check(ratios, toy_values, analytic_values, path_stem):
    """Toy ring average of a coaxial ellipse against its closed form."""
    fig, (ax, ax_diff) = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True,
                                      gridspec_kw=dict(height_ratios=[3, 1]))
    ax.plot(ratios, analytic_values, color="#000000", label=r"$(b/a)\,K(k)/E(k)$")
    ax.plot(ratios, toy_values, linestyle="none", marker="o", markersize=4, color="#cc0000", label="toy")
    ax.set_ylabel(LABEL_RING)
    ax.legend()
    # A difference, not a ratio: both curves vanish as b/a --> 0.
    ax_diff.plot(ratios, np.asarray(toy_values) - np.asarray(analytic_values), marker="o", markersize=3,
                 color="#cc0000")
    ax_diff.axhline(0.0, color="#888888", linewidth=0.6)
    ax_diff.set_ylabel("toy - analytic")
    ax_diff.set_xlabel("minor / major semi-axis $b/a$")
    save_figure(fig, path_stem)
