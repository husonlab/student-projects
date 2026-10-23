#!/usr/bin/env python3
"""
Figures for the format comparison (bench_formats.py, results/formats_<N>.csv).

  formats_tradeoff.pdf/.png   size on disk vs. warm row-query latency, one point per variant
  formats_latency.pdf/.png    row and column latency per variant, first pass (from disk)
                              and second pass (from the page cache)

Both use the largest sample (200M input rows). Values are medians over the
three repeats of the per-pass medians.

Usage: python3 plot_formats.py [formats_200000000.csv]
"""
import csv
import statistics as st
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "formats_200000000.csv"

# Chart chrome (reference palette, light mode) on white paper, as in plot_buffer_sweep.py.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# One hue per format, fixed slot order of the categorical palette; marker shapes as
# secondary encoding (two of the hues are below 3:1 contrast, so every point is also labelled).
FORMATS = {
    "MAGrix": ("#2a78d6", "o"),
    "HDF5": ("#eb6834", "s"),
    "Zarr": ("#1baf7a", "^"),
    "SQLite": ("#eda100", "D"),
}

# variant in the CSV -> (format, label), in the order shown in formats_latency
VARIANTS = [
    ("magrix", "MAGrix", "MAGrix container"),
    ("hdf5_contiguous", "HDF5", "HDF5 contiguous"),
    ("hdf5_chunk64k", "HDF5", "HDF5 64k chunks"),
    ("hdf5_chunk64k_gzip", "HDF5", "HDF5 64k + gzip"),
    ("hdf5_chunk2M_gzip", "HDF5", "HDF5 2M chunks + gzip"),
    ("zarr_chunk64k", "Zarr", "Zarr 64k chunks"),
    ("zarr_chunk64k_zstd", "Zarr", "Zarr 64k + zstd"),
    ("zarr_chunk2M_zstd", "Zarr", "Zarr 2M chunks + zstd"),
    ("sqlite_rowid", "SQLite", "SQLite, one index per id"),
    ("sqlite_norowid", "SQLite", "SQLite WITHOUT ROWID"),
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


def load(path):
    """variant -> {column: median over repeats}"""
    rows = list(csv.DictReader(open(path)))
    out = {}
    for variant, _, _ in VARIANTS:
        reps = [r for r in rows if r["variant"] == variant]
        out[variant] = {k: st.median(float(r[k]) for r in reps)
                        for k in reps[0] if k not in ("variant",) and reps[0][k] != ""}
    return out, int(rows[0]["nnz"])


def time_label(us, _pos=None):
    if us >= 1000:
        return f"{us / 1000:g} ms"
    return f"{us:g} µs"


def log_time_axis(axis, lo, hi):
    ticks = [t for t in (1, 10, 100, 1_000, 10_000, 100_000) if lo <= t <= hi]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(time_label))
    axis.set_minor_locator(NullLocator())


def tradeoff_figure(data, nnz):
    fig, ax = plt.subplots(figsize=(6.3, 3.6))
    # label offsets in points, chosen so that labels of neighbouring points don't collide
    offsets = {
        "magrix": (7, 0, "left"),
        "hdf5_contiguous": (7, -5, "left"),
        "hdf5_chunk64k": (7, 5, "left"),
        "hdf5_chunk64k_gzip": (7, -6, "left"),
        "hdf5_chunk2M_gzip": (7, 0, "left"),
        "zarr_chunk64k": (7, 0, "left"),
        "zarr_chunk64k_zstd": (-7, 4, "right"),
        "zarr_chunk2M_zstd": (-7, 0, "right"),
        "sqlite_rowid": (0, -11, "center"),
        "sqlite_norowid": (0, 9, "center"),
    }
    for variant, fmt, label in VARIANTS:
        color, marker = FORMATS[fmt]
        x = data[variant]["bytes"] / 1e9
        y = data[variant]["row_second_median_us"]
        ax.plot(x, y, marker=marker, color=color, markersize=7, markeredgecolor="white",
                markeredgewidth=1.0, linestyle="none", zorder=3)
        dx, dy, ha = offsets[variant]
        ax.annotate(label, (x, y), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center", fontsize=7, color=INK_2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 13)  # room on the left for the labels of the compressed Zarr variants
    ax.set_ylim(2, 100_000)
    ax.xaxis.set_major_locator(FixedLocator([0.2, 0.5, 1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    log_time_axis(ax.yaxis, 1, 100_000)
    ax.set_xlabel("size on disk (GB, log scale)")
    ax.set_ylabel("row query, from page cache (median, log scale)")
    ax.set_title(f"Size vs. query latency, {nnz / 1e6:.0f} million non-zero entries "
                 "(smaller and lower is better)", loc="left")
    ax.legend(handles=format_handles(), loc="upper right", frameon=False, ncol=4,
              bbox_to_anchor=(1.0, 1.0), handletextpad=0.3, columnspacing=1.2)
    fig.tight_layout()
    return fig


def format_handles():
    return [Line2D([], [], marker=m, color=c, linestyle="none", markersize=6,
                   markeredgecolor="white", label=f) for f, (c, m) in FORMATS.items()]


def latency_figure(data, nnz):
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.8), sharey=True)
    ys = list(range(len(VARIANTS)))[::-1]  # first variant at the top
    for ax, kind, title in [(axes[0], "row", "Row query (families of a MAG)"),
                            (axes[1], "col", "Column query (MAGs of a family)")]:
        for y, (variant, fmt, _) in zip(ys, VARIANTS):
            color, marker = FORMATS[fmt]
            cold = data[variant][f"{kind}_first_median_us"]
            warm = data[variant][f"{kind}_second_median_us"]
            ax.plot([warm, cold], [y, y], color=AXIS, linewidth=1.2, zorder=1)
            ax.plot(cold, y, marker=marker, markersize=6.5, markerfacecolor="white",
                    markeredgecolor=color, markeredgewidth=1.4, linestyle="none", zorder=3)
            ax.plot(warm, y, marker=marker, markersize=6.5, color=color,
                    markeredgecolor="white", markeredgewidth=0.8, linestyle="none", zorder=4)
        ax.set_xscale("log")
        ax.set_xlim(1.5, 150_000)
        log_time_axis(ax.xaxis, 1, 100_000)
        ax.grid(axis="y", visible=False)
        ax.grid(axis="x", visible=True)
        ax.set_xlabel("latency (median, log scale)")
        ax.set_title(title, loc="left")
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([label for _, _, label in VARIANTS])
    axes[0].tick_params(axis="y", length=0)
    axes[1].tick_params(axis="y", length=0)
    pass_handles = [
        Line2D([], [], marker="o", color=INK_2, linestyle="none", markersize=6,
               markeredgecolor="white", label="from page cache (2nd pass)"),
        Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=INK_2,
               markeredgewidth=1.4, linestyle="none", markersize=6, label="from disk (1st pass)"),
    ]
    fig.legend(handles=format_handles() + pass_handles, loc="upper center", ncol=6, frameon=False,
               bbox_to_anchor=(0.55, 1.04), handletextpad=0.3, columnspacing=1.0, fontsize=7)
    fig.tight_layout(w_pad=1.5)
    return fig


def main():
    data, nnz = load(CSV)
    for name, fig in [("formats_tradeoff", tradeoff_figure(data, nnz)),
                      ("formats_latency", latency_figure(data, nnz))]:
        for ext in ("pdf", "png"):
            fig.savefig(HERE / f"{name}.{ext}")
        print("wrote", HERE / f"{name}.pdf", "and .png")


if __name__ == "__main__":
    main()
