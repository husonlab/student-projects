#!/usr/bin/env python3
"""
Figures for the buffer sweep (results/buffer_sweep.csv).

  buffer_sweep_io.pdf/.png            run time and peak memory over read/write buffer size
  buffer_sweep_rows_per_run.pdf/.png  run time and peak memory of stage 2b over rows_per_run

Points are medians over the repeats, vertical bars span min-max.

Usage: python3 plot_buffer_sweep.py [buffer_sweep.csv]
"""
import csv
import math
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "buffer_sweep.csv"

# Defaults at the time of the sweep, and the defaults chosen from it.
OLD_READ, OLD_WRITE, NEW_BUF = 64, 32, 8
MERGE_BUF = 4 * 1024 * 1024        # bytes per run during the merge
RECORD_BYTES_IN_MEMORY = 12        # (u32, u32, u16) tuple incl. padding
NNZ_SAMPLE = 49_660_727            # non-zero entries of the 50M-row sample
DEFAULT_ROWS_PER_RUN = 20_000_000

# Chart chrome (reference palette, light mode) on white paper.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"

# Fixed slot order of the categorical palette; marker shapes as secondary encoding.
STAGES = [
    ("1_main", "1 preprocessing", "#2a78d6", "o"),
    ("2a_csc", "2a CSC", "#eb6834", "s"),
    ("2b_csr", "2b CSR", "#1baf7a", "^"),
    ("3_members", "3 members", "#eda100", "D"),
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
    "axes.grid.axis": "y",
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "grid.linestyle": "-",
    "axes.axisbelow": True,
    "savefig.bbox": "tight",
    "savefig.dpi": 200,
    "pdf.fonttype": 42,
})


def load(path):
    rows = list(csv.DictReader(open(path)))
    by_stage = defaultdict(list)       # (experiment, config, stage) -> [(sec, mb)]
    totals = defaultdict(float)        # (config, repeat) -> total seconds of stages 1-4
    for r in rows:
        sec, mb = float(r["seconds"]), int(r["max_rss_kb"]) / 1024
        by_stage[(r["experiment"], r["config"], r["stage"])].append((sec, mb))
        if r["experiment"] == "io":
            totals[(r["config"], r["repeat"])] += sec
    total_by_config = defaultdict(list)
    for (config, _), t in totals.items():
        total_by_config[config].append(t)
    return by_stage, total_by_config


def summary(values):
    return st.median(values), min(values), max(values)


def config_for(kind, size):
    old = OLD_READ if kind == "read" else OLD_WRITE
    return "default" if size == old else f"{kind}_{size}MiB"


def style_log2_axis(ax, sizes, unit="MiB"):
    ax.set_xscale("log", base=2)
    ax.set_xticks(sizes)
    ax.set_xticklabels([str(s) for s in sizes])
    ax.minorticks_off()
    ax.set_xlabel(f"buffer size ({unit}, log scale)")


def mark_default(ax, x, label, strong=True, side="right", y=0.97):
    """Vertical reference line; the label sits left or right of it so two lines don't collide."""
    ax.axvline(x, color=INK_2 if strong else AXIS, linewidth=0.8, zorder=0)
    text = f" {label}" if side == "right" else f"{label} "
    ax.text(x, y, text, transform=ax.get_xaxis_transform(),
            color=INK_2 if strong else MUTED, fontsize=7, va="top",
            ha="left" if side == "right" else "right")


def plot_with_range(ax, xs, stats, color, marker, label=None):
    med = [s[0] for s in stats]
    ax.vlines(xs, [s[1] for s in stats], [s[2] for s in stats], color=color, linewidth=1.0, zorder=2)
    ax.plot(xs, med, color=color, linewidth=1.6, marker=marker, markersize=5,
            markeredgecolor="white", markeredgewidth=0.8, label=label, zorder=3)
    return med


def io_figure(by_stage, totals):
    sizes = {"read": [1, 2, 4, 8, 16, 32, 64, 128, 256], "write": [1, 2, 4, 8, 16, 32, 64, 128]}
    fig, axes = plt.subplots(2, 2, figsize=(6.3, 5.0), sharey="col")
    for row, kind in enumerate(["read", "write"]):
        xs = sizes[kind]
        old = OLD_READ if kind == "read" else OLD_WRITE

        # left: total run time of the pipeline (neutral ink: it is not one of the stages)
        ax = axes[row][0]
        stats = [summary(totals[config_for(kind, s)]) for s in xs]
        plot_with_range(ax, xs, stats, INK_2, "o")
        ax.set_ylim(0, 80)
        ax.set_ylabel("total run time (s)")
        ax.set_title(f"{kind.capitalize()} buffer: total run time", loc="left")
        style_log2_axis(ax, xs)
        mark_default(ax, NEW_BUF, "new default", side="left")
        mark_default(ax, old, "previous", strong=False)

        # right: peak memory per stage, direct labels at the right end
        ax = axes[row][1]
        ends = []
        for key, name, color, marker in STAGES:
            stats = [summary([v[1] for v in by_stage[("io", config_for(kind, s), key)]]) for s in xs]
            med = plot_with_range(ax, xs, stats, color, marker, label=name)
            ends.append((med[-1], name))
        ax.set_ylim(0, 600)
        ax.set_ylabel("peak memory (MB)")
        ax.set_title(f"{kind.capitalize()} buffer: peak memory per stage", loc="left")
        style_log2_axis(ax, xs)
        mark_default(ax, NEW_BUF, "new default", side="left")
        mark_default(ax, old, "previous", strong=False)
        # spread end labels so close lines don't collide
        placed = []
        for y, name in sorted(ends):
            y_lab = max(y, placed[-1] + 34) if placed else y
            placed.append(y_lab)
            ax.annotate(name, xy=(xs[-1], y), xytext=(xs[-1] * 1.18, y_lab),
                        color=INK_2, fontsize=7, va="center", annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", color=AXIS, linewidth=0.6))
        ax.set_xlim(xs[0] / 1.3, xs[-1] * 1.25)

    handles, labels = axes[0][1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(h_pad=2.2, w_pad=3.2)
    return fig


def rows_per_run_figure(by_stage):
    values = [1, 2, 4, 8, 16, 32, 64]                       # millions
    points = []
    for v in values:
        rec = by_stage[("rows_per_run", str(v * 1_000_000), "2b_csr")]
        runs = math.ceil(NNZ_SAMPLE / (v * 1_000_000))
        points.append((v, summary([r[0] for r in rec]), summary([r[1] for r in rec]), runs))
    # the default (20M) was measured in the io experiment with the same buffers
    rec = by_stage[("io", "default", "2b_csr")]
    points.append((20, summary([r[0] for r in rec]), summary([r[1] for r in rec]),
                   math.ceil(NNZ_SAMPLE / DEFAULT_ROWS_PER_RUN)))
    points.sort()
    xs = [p[0] for p in points]
    optimum = math.sqrt(NNZ_SAMPLE * MERGE_BUF / RECORD_BYTES_IN_MEMORY) / 1e6

    fig, axes = plt.subplots(1, 2, figsize=(6.3, 2.9))
    for ax, idx, ylabel, title, ymax in [
        (axes[0], 1, "run time (s)", "Stage 2b: run time", 14),
        (axes[1], 2, "peak memory (MB)", "Stage 2b: peak memory", 700),
    ]:
        plot_with_range(ax, xs, [p[idx] for p in points], BLUE, "o")
        ax.set_ylim(0, ymax)
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", pad=18)
        ax.set_xscale("log", base=2)
        ax.set_xticks(values)
        ax.set_xticklabels([f"{v}M" for v in values])
        ax.minorticks_off()
        ax.set_xlabel("entries per run (rows_per_run)")
        mark_default(ax, DEFAULT_ROWS_PER_RUN / 1e6, "default (20M)", strong=False, side="left", y=0.86)
        # number of runs on a secondary x axis (same positions, so it's not a second scale)
        top = ax.secondary_xaxis("top")
        top.set_xticks(values)
        top.set_xticks([], minor=True)
        top.set_xticklabels([str(math.ceil(NNZ_SAMPLE / (v * 1_000_000))) for v in values])
        top.tick_params(labelsize=7, colors=MUTED, labelcolor=MUTED, length=2)
        top.spines["top"].set_visible(False)
        ax.text(1.0, 1.02, "runs", transform=ax.transAxes, ha="left", va="bottom", fontsize=7, color=MUTED)
    mark_default(axes[1], optimum, f"predicted minimum ≈ {optimum:.1f}M")
    fig.tight_layout(w_pad=3.0)
    return fig


def main():
    by_stage, totals = load(CSV)
    for name, fig in [("buffer_sweep_io", io_figure(by_stage, totals)),
                      ("buffer_sweep_rows_per_run", rows_per_run_figure(by_stage))]:
        for ext in ("pdf", "png"):
            fig.savefig(HERE / f"{name}.{ext}")
        print("wrote", HERE / f"{name}.pdf", "and .png")


if __name__ == "__main__":
    main()
