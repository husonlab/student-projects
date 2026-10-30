#!/usr/bin/env python3
"""
Figures for the query latency (test_real_data.py) and for the format comparison
over the sample size (bench_formats.py).

  query_latency.pdf/.png    median latency per query type, first access (from disk)
                            and repeated (from the page cache), latency_50id80cov.csv
  formats_scaling.pdf/.png  size on disk and warm row/column latency over the number
                            of non-zero entries, one line per format, formats_<N>.csv

The format samples are the first N input rows, so they all contain (almost) all
MAGs and differ in the number of families. Values are medians over the three
repeats of the per-pass medians.

Usage: python3 plot_queries.py
"""
import csv
import statistics as st
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
LATENCY_CSV = HERE / "latency_50id80cov.csv"
FORMAT_CSVS = sorted(HERE.glob("formats_*.csv"), key=lambda p: int(p.stem.split("_")[1]))

# Chart chrome (reference palette, light mode) on white paper, as in plot_buffer_sweep.py.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
BLUE = "#2a78d6"

# query in the CSV -> label, in the order shown (top to bottom)
QUERIES = [
    ("row (families of a MAG)", "Row (families of a MAG)"),
    ("column (MAGs of a family)", "Column (MAGs of a family)"),
    ("get (one cell)", "Single cell"),
    ("members (genes of a cell)", "Member genes of a cell"),
    ("members worst case (largest family, last cell)", "Member genes, largest family"),
    ("submatrix (10 MAGs x 100 families)", "Submatrix (10 MAGs × 100 families)"),
]

# Same hue and marker per format as in plot_formats.py; line style tells the two
# HDF5 variants apart. One variant per format to keep the lines readable.
SCALING_VARIANTS = [
    ("magrix", "MAGrix", "#2a78d6", "o", "-"),
    ("hdf5_contiguous", "HDF5 contiguous", "#eb6834", "s", "-"),
    ("hdf5_chunk64k_gzip", "HDF5 64k + gzip", "#eb6834", "s", "--"),
    ("zarr_chunk64k_zstd", "Zarr 64k + zstd", "#1baf7a", "^", "-"),
    ("sqlite_norowid", "SQLite WITHOUT ROWID", "#eda100", "D", "-"),
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


def log_time_axis(axis, lo, hi):
    ticks = [t for t in (1, 10, 100, 1_000, 10_000, 100_000) if lo <= t <= hi]
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FuncFormatter(time_label))
    axis.set_minor_locator(NullLocator())


def latency_figure():
    rows = {(r["query"], r["pass"]): float(r["median_us"]) for r in csv.DictReader(open(LATENCY_CSV))}
    fig, ax = plt.subplots(figsize=(6.3, 2.6))
    ys = list(range(len(QUERIES)))[::-1]  # first query at the top
    for y, (query, _) in zip(ys, QUERIES):
        warm = rows[(query, "warm")]
        cold = rows.get((query, "cold"))  # the worst case was only measured warm
        if cold is not None:
            ax.plot([warm, cold], [y, y], color=AXIS, linewidth=1.2, zorder=1)
            ax.plot(cold, y, marker="o", markersize=6.5, markerfacecolor="white",
                    markeredgecolor=BLUE, markeredgewidth=1.4, linestyle="none", zorder=3)
            ax.annotate(time_label(round(cold, -1) if cold < 10_000 else round(cold, -2)),
                        (cold, y), xytext=(7, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=7, color=INK_2)
        ax.plot(warm, y, marker="o", markersize=6.5, color=BLUE,
                markeredgecolor="white", markeredgewidth=0.8, linestyle="none", zorder=4)
        ax.annotate(time_label(round(warm, 1) if warm < 100 else round(warm)),
                    (warm, y), xytext=(-7, 0), textcoords="offset points",
                    ha="right", va="center", fontsize=7, color=INK_2)
    ax.set_xscale("log")
    ax.set_xlim(0.8, 150_000)
    log_time_axis(ax.xaxis, 1, 100_000)
    ax.grid(axis="y", visible=False)
    ax.set_yticks(ys)
    ax.set_yticklabels([label for _, label in QUERIES])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("latency (median, log scale)")
    ax.set_title("Query latency on the full dataset", loc="left")
    handles = [
        Line2D([], [], marker="o", color=BLUE, linestyle="none", markersize=6,
               markeredgecolor="white", label="from page cache (repeated)"),
        Line2D([], [], marker="o", markerfacecolor="white", markeredgecolor=BLUE,
               markeredgewidth=1.4, linestyle="none", markersize=6, label="from disk (first access)"),
    ]
    # in the empty space right of the worst-case row (measured warm only)
    ax.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.27), frameon=False,
              handletextpad=0.3, fontsize=7)
    fig.tight_layout()
    return fig


def load_formats():
    """variant -> list of (nnz, {column: median over repeats}), sorted by nnz"""
    out = {v[0]: [] for v in SCALING_VARIANTS}
    for path in FORMAT_CSVS:
        rows = list(csv.DictReader(open(path)))
        for variant in out:
            reps = [r for r in rows if r["variant"] == variant]
            med = {k: st.median(float(r[k]) for r in reps)
                   for k in ("bytes", "row_second_median_us", "col_second_median_us")}
            out[variant].append((int(reps[0]["nnz"]), med))
    return out


def scaling_figure(data):
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.7))
    panels = [
        (axes[0], "bytes", 1e-9, "Size on disk (GB)"),
        (axes[1], "row_second_median_us", 1, "Row query, page cache"),
        (axes[2], "col_second_median_us", 1, "Column query, page cache"),
    ]
    for ax, key, scale, title in panels:
        # MAGrix last, so it is drawn on top of HDF5 contiguous (same size on disk)
        for variant, _, color, marker, style in SCALING_VARIANTS[1:] + SCALING_VARIANTS[:1]:
            xs = [nnz / 1e6 for nnz, _ in data[variant]]
            ys = [m[key] * scale for _, m in data[variant]]
            ax.plot(xs, ys, color=color, linestyle=style, linewidth=2, marker=marker,
                    markersize=6, markeredgecolor="white", markeredgewidth=1.0, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        nnz_ticks = [nnz / 1e6 for nnz, _ in data["magrix"]]
        ax.xaxis.set_major_locator(FixedLocator(nnz_ticks))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlim(nnz_ticks[0] / 1.5, nnz_ticks[-1] * 1.5)
        ax.set_xlabel("non-zero entries (millions)")
        ax.set_title(title, loc="left")
    axes[0].yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10]))
    axes[0].yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    axes[0].yaxis.set_minor_locator(NullLocator())
    axes[0].set_ylim(0.01, 10)
    axes[0].text(0.04, 0.97, "HDF5 contiguous has the\nsame size as MAGrix", transform=axes[0].transAxes,
                 ha="left", va="top", fontsize=7, color=INK_2)
    for ax in axes[1:]:
        ax.set_ylim(2, 10_000)
        log_time_axis(ax.yaxis, 1, 10_000)
    handles = [Line2D([], [], color=c, linestyle=s, linewidth=2, marker=m, markersize=5,
                      markeredgecolor="white", label=label)
               for _, label, c, m, s in SCALING_VARIANTS]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False,
               bbox_to_anchor=(0.5, 1.07), handletextpad=0.4, columnspacing=1.0, fontsize=7)
    fig.tight_layout(w_pad=1.2)
    return fig


def main():
    for name, fig in [("query_latency", latency_figure()),
                      ("formats_scaling", scaling_figure(load_formats()))]:
        for ext in ("pdf", "png"):
            fig.savefig(HERE / f"{name}.{ext}")
        print("wrote", HERE / f"{name}.pdf", "and .png")


if __name__ == "__main__":
    main()
