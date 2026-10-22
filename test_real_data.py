#!/usr/bin/env python3
"""
Sanity and performance check of the magrix query API on a real output directory.

Measures the time to open a dataset and the latency of each query type. Every
query type is run twice with the same random IDs: the first pass mostly reads
from disk ("cold"), the second finds the data in the page cache ("warm").
Also checks on a sample that member genes and copy numbers agree.

Usage:
    python3 test_real_data.py <data_dir> [--queries N] [--seed S] [--csv results.csv]
"""
import argparse
import csv
import random
import resource
import statistics as st
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from magrix import PanGenomeMatrices  # noqa: E402
from magrix.api import Container  # noqa: E402


def peak_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # Linux: KB


def latencies(fn, args_list):
    out = []
    for args in args_list:
        t0 = time.perf_counter()
        fn(*args)
        out.append(time.perf_counter() - t0)
    return out


def describe(times):
    us = sorted(t * 1e6 for t in times)
    p95 = us[min(len(us) - 1, int(round(0.95 * (len(us) - 1))))]
    return st.median(us), p95, us[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_dir")
    parser.add_argument("--queries", type=int, default=2000, help="queries per type (default 2000)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--csv", help="also write the latency results to this CSV file")
    args = parser.parse_args()
    data_dir, n = args.data_dir, args.queries

    # ── opening ──────────────────────────────────────────────────────────────
    # The container alone is cheap to open, so timing it separately costs nothing;
    # the rest of PanGenomeMatrices' open time is loading the dictionaries.
    t0 = time.perf_counter()
    Container(f"{data_dir}/paralog.container")
    t_container = time.perf_counter() - t0
    t0 = time.perf_counter()
    m = PanGenomeMatrices(data_dir)
    t_open = time.perf_counter() - t0
    print(f"Dataset: {data_dir}")
    print(f"  {m.dicts.num_mags:,} MAGs, {m.dicts.num_families:,} gene families, "
          f"{len(m.paralog.csc_indices):,} non-zero entries")
    print(f"Open: {t_open:.2f} s in total (container header {t_container*1000:.1f} ms, "
          f"the rest is loading the dictionaries); memory after opening {peak_rss_mb():.0f} MB\n")

    # ── query workload (same random IDs for both passes) ─────────────────────
    rng = random.Random(args.seed)
    mag_ids = [rng.randrange(m.dicts.num_mags) for _ in range(n)]
    family_ids = [rng.randrange(m.dicts.num_families) for _ in range(n)]

    def random_cell():
        while True:
            f = rng.randrange(m.dicts.num_families)
            s, e = int(m.paralog.csc_indptr[f]), int(m.paralog.csc_indptr[f + 1])
            if e > s:
                return int(m.paralog.csc_indices[rng.randrange(s, e)]), f

    cells = [random_cell() for _ in range(n)]
    n_sub = max(1, n // 20)
    submatrices = [([rng.randrange(m.dicts.num_mags) for _ in range(10)],
                    [rng.randrange(m.dicts.num_families) for _ in range(100)]) for _ in range(n_sub)]

    workloads = [
        ("row (families of a MAG)", m.paralog.row, [(i,) for i in mag_ids]),
        ("column (MAGs of a family)", m.paralog.column, [(j,) for j in family_ids]),
        ("get (one cell)", m.paralog.get, cells),
        ("members (genes of a cell)", m.members.gene_ids_for, cells),
        ("submatrix (10 MAGs x 100 families)", m.paralog.submatrix, submatrices),
    ]

    rows = []
    header = f"{'query':38} {'pass':5} {'n':>5} {'median us':>11} {'p95 us':>11} {'max us':>11}"
    print(header)
    print("-" * len(header))
    for name, fn, arglist in workloads:
        for pass_name in ("cold", "warm"):
            med, p95, mx = describe(latencies(fn, arglist))
            rows.append((name, pass_name, len(arglist), med, p95, mx))
            print(f"{name:38} {pass_name:5} {len(arglist):>5} {med:>11.1f} {p95:>11.1f} {mx:>11.1f}")

    # ── worst case for members: last cell of the largest family ─────────────
    prevalence = np.diff(m.paralog.csc_indptr)
    big = int(prevalence.argmax())
    last_mag = int(m.paralog.csc_indices[int(m.paralog.csc_indptr[big + 1]) - 1])
    worst = latencies(m.members.gene_ids_for, [(last_mag, big)] * 50)
    med, p95, mx = describe(worst[1:])  # first call reads from disk
    print(f"\nmembers worst case: family {big} carried by {int(prevalence[big]):,} MAGs, "
          f"last cell: first call {worst[0]*1e6:.0f} us, then median {med:.0f} us")
    rows.append(("members worst case (largest family, last cell)", "warm", len(worst) - 1, med, p95, mx))

    # ── consistency spot check ───────────────────────────────────────────────
    bad = 0
    for mag, fam in cells[:500]:
        ids = m.members.gene_ids_for(mag, fam)
        if len(ids) != m.paralog.get(mag, fam) or np.any((ids >> np.uint64(32)) != mag):
            bad += 1
    print(f"Consistency: {bad} mismatches in 500 cells (member genes vs. copy number and MAG)")

    mag, fam = cells[0]
    print(f"Example: {m.dicts.mag_name(mag)} in family {fam} "
          f"({m.dicts.family_representative_gene(fam)}): copy number {m.paralog.get(mag, fam)}, "
          f"genes {m.members.genes_for(mag, fam)[:3]}")
    print(f"Peak memory of this process: {peak_rss_mb():.0f} MB")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["query", "pass", "n", "median_us", "p95_us", "max_us"])
            w.writerows(rows)
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
