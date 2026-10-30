#!/usr/bin/env python3
"""
Figures for the scaling over the number of MAGs (random, nested MAG subsets of the
full dataset, each built with magrix build and measured with test_real_data.py).

  scaling_latency.pdf/.png  median query latency over the number of MAGs, first access
                            (from disk) and repeated (from the page cache), scaling_latency.csv
  scaling_build.pdf/.png    build time, peak memory, gene families and container size over
                            the number of MAGs, scaling_summary.tsv

The full dataset was not rebuilt in the scaling run; its build time and memory are
taken from a separate build of the same data on the same cluster (FULL_BUILD_*).

Usage: python3 plot_scaling.py
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
LATENCY_CSV = HERE / "scaling_latency.csv"
SUMMARY_TSV = HERE / "scaling_summary.tsv"

# Build of the full dataset in a separate run (alone on the node): total time and peak memory.
FULL_BUILD_S = 2590
FULL_BUILD_PEAK_MB = 888

# Chart chrome (reference palette, light mode) on white paper, as in plot_buffer_sweep.py.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"

# query in the CSV -> (label, colour, marker); fixed slot order of the categorical palette
QUERIES = [
    ("row (families of a MAG)", "Row", "#2a78d6", "o"),
    ("column (MAGs of a family)", "Column", "#eb6834", "s"),
    ("get (one cell)", "Single cell", "#1baf7a", "^"),
    ("members (genes of a cell)", "Member genes", "#eda100", "D"),
    ("members worst case (largest family, last cell)", "Member genes, largest family", "#e87ba4", "v"),
    ("submatrix (10 MAGs x 100 families)", "Submatrix 10 × 100", "#008300", "P"),
]

plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "text.color": INK,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "axes.axisbelow": True,
    "savefig.bbox": "tight",
    "savefig.dpi": 200,
    "pdf.fonttype": 42,
})


def time_label(us, _pos=None):
    if us >= 1000:
        return f"{us / 1000:g} ms"
    return f"{us:g} µs"


def mags_axis(ax, sizes, log):
    if log:
        # one tick per measured size; the label of the full dataset one line lower,
        # it would collide with 400k on the log axis
        ax.set_xscale("log")
        ax.set_xlim(sizes[0] / 1.4, sizes[-1] * 1.4)
        ax.xaxis.set_major_locator(FixedLocator(sizes))
        ax.xaxis.set_major_formatter(FuncFormatter(
            lambda n, _: ("\n" if n == sizes[-1] else "") + f"{n / 1000:.0f}k"))
    else:
        ax.set_xlim(0, sizes[-1] * 1.05)
        ax.xaxis.set_major_locator(FixedLocator([0, 100_000, 200_000, 300_000, 400_000, 500_000]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda n, _: f"{n / 1000:.0f}k" if n else "0"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("MAGs")


def load_latency():
    """(query, pass) -> list of (num_mags, median_us), sorted by num_mags"""
    out = {}
    for r in csv.DictReader(open(LATENCY_CSV)):
        out.setdefault((r["query"], r["pass"]), []).append((int(r["num_mags"]), float(r["median_us"])))
    return {k: sorted(v) for k, v in out.items()}


def load_summary():
    rows = list(csv.DictReader(open(SUMMARY_TSV), delimiter="\t"))
    return {k: [float(r[k]) if r[k] else None for r in rows] for k in rows[0]}


def latency_figure(data):
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.0), sharey=True)
    for ax, kind, title in [(axes[0], "cold", "First access (from disk)"),
                            (axes[1], "warm", "Repeated (from the page cache)")]:
        for query, _, color, marker in QUERIES:
            points = data.get((query, kind))
            if not points:  # the worst case was only measured warm
                continue
            ax.plot([n for n, _ in points], [us for _, us in points], color=color, linewidth=2,
                    marker=marker, markersize=6, markeredgecolor="white", markeredgewidth=1.0, zorder=3)
        sizes = [n for n, _ in data[(QUERIES[0][0], "warm")]]
        mags_axis(ax, sizes, log=True)
        ax.set_yscale("log")
        ax.set_title(title, loc="left")
    axes[0].set_ylim(2, 50_000)
    axes[0].yaxis.set_major_locator(FixedLocator([10, 100, 1_000, 10_000]))
    axes[0].yaxis.set_major_formatter(FuncFormatter(time_label))
    axes[0].yaxis.set_minor_locator(NullLocator())
    axes[0].set_ylabel("latency (median, log scale)")
    handles = [Line2D([], [], color=c, linewidth=2, marker=m, markersize=5, markeredgecolor="white", label=label)
               for _, label, c, m in QUERIES]
    fig.legend(handles=handles, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 1.06),
               handletextpad=0.4, columnspacing=1.0, fontsize=7)
    fig.tight_layout(w_pad=1.5)
    return fig


def scaling_panel(ax, mags, values, title, fmt, full_value=None):
    xs = [m for m, v in zip(mags, values) if v is not None]
    ys = [v for v in values if v is not None]
    ax.plot(xs, ys, color=BLUE, linewidth=2, marker="o", markersize=6,
            markeredgecolor="white", markeredgewidth=1.0, zorder=3)
    for x, y in zip(xs, ys):
        if x < 100_000:  # the small subsets sit too close together for labels
            continue
        ax.annotate(fmt(y), (x, y), xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom", fontsize=6.5, color=INK_2)
    if full_value is not None:  # full dataset from a separate build
        ax.plot([xs[-1], mags[-1]], [ys[-1], full_value], color=BLUE, linewidth=1.2, linestyle=":", zorder=2)
        ax.plot(mags[-1], full_value, marker="o", markersize=6, markerfacecolor="white",
                markeredgecolor=BLUE, markeredgewidth=1.4, linestyle="none", zorder=3)
        ax.annotate(fmt(full_value), (mags[-1], full_value), xytext=(-5, 5), textcoords="offset points",
                    ha="right", va="bottom", fontsize=6.5, color=INK_2)
    mags_axis(ax, mags, log=False)
    ax.set_ylim(0, max(ys + [full_value or 0]) * 1.2)
    ax.set_title(title, loc="left")


def build_figure(s):
    mags = [int(m) for m in s["num_mags"]]
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 4.6))
    scaling_panel(axes[0][0], mags, [v / 60 if v else None for v in s["build_s"]],
                  "Build time (min)", lambda v: f"{v:.1f}" if v < 10 else f"{v:.0f}", FULL_BUILD_S / 60)
    scaling_panel(axes[0][1], mags, s["build_peak_mb"], "Peak memory of the build (MB)",
                  lambda v: f"{v:.0f}", FULL_BUILD_PEAK_MB)
    scaling_panel(axes[1][0], mags, [v / 1e6 for v in s["num_families"]], "Gene families (millions)",
                  lambda v: f"{v:.1f}" if v < 10 else f"{v:.0f}")
    scaling_panel(axes[1][1], mags, [v / 1e9 for v in s["container_bytes"]], "Container size (GB)",
                  lambda v: f"{v:.1f}")

    # Power law families = k * MAGs^gamma (Heaps' law), least squares in log space
    gamma, log_k = np.polyfit(np.log(mags), np.log(s["num_families"]), 1)
    xs = np.linspace(mags[0], mags[-1], 200)
    ax = axes[1][0]
    ax.plot(xs, np.exp(log_k) * xs ** gamma / 1e6, color=MUTED, linewidth=1, linestyle="--", zorder=2)
    ax.text(0.97, 0.06, f"fit: families ∝ MAGs$^{{{gamma:.2f}}}$", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color=INK_2)

    handles = [
        Line2D([], [], color=BLUE, linewidth=2, marker="o", markersize=5, markeredgecolor="white",
               label="MAG subsets (one run)"),
        Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE, markeredgewidth=1.4,
               linestyle="none", markersize=5, label="full dataset, separate build"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.04),
               handletextpad=0.4, columnspacing=1.5, fontsize=7)
    fig.tight_layout(h_pad=1.5, w_pad=1.5)
    print(f"families ~ MAGs^{gamma:.3f}")
    return fig


def main():
    for name, fig in [("scaling_latency", latency_figure(load_latency())),
                      ("scaling_build", build_figure(load_summary()))]:
        for ext in ("pdf", "png"):
            fig.savefig(HERE / f"{name}.{ext}")
        print("wrote", HERE / f"{name}.pdf", "and .png")


if __name__ == "__main__":
    main()
