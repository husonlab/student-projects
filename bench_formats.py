#!/usr/bin/env python3
"""
Fair comparison of the MAGrix container with HDF5, Zarr and SQLite for the access
pattern MAGrix is built for: many independent single-row and single-column reads.

Every format stores exactly the same data, copied from an existing paralog.container:
the CSR and CSC arrays of the paralog matrix (indptr u64, indices u32, copy numbers u16).
The members arrays are left out, since the other formats would not store them either.

Variants:
  magrix               the container itself, read with numpy.memmap like the magrix API
  hdf5_contiguous      HDF5, no chunking, no compression
  hdf5_chunk64k        HDF5, chunks of 65,536 elements, no compression
  hdf5_chunk64k_gzip   HDF5, chunks of 65,536 elements, gzip level 4
  hdf5_chunk2M_gzip    HDF5, chunks of 2,000,000 elements, gzip level 4 (earlier benchmark)
  zarr_chunk64k        Zarr, chunks of 65,536 elements, no compression
  zarr_chunk64k_zstd   Zarr, chunks of 65,536 elements, Blosc zstd level 5
  zarr_chunk2M_zstd    Zarr, chunks of 2,000,000 elements, Blosc zstd level 5 (earlier benchmark)
  sqlite_rowid         table (mag_id, family_id, value) + one index per id (earlier benchmark)
  sqlite_norowid       WITHOUT ROWID table, primary key (mag_id, family_id), plus a covering
                       index (family_id, value) -- the best plain SQLite layout for both queries

For each variant the script measures the size on disk, the time and peak memory to build
it from the container (in its own process), and then --repeats times, each in a fresh
process: the time to open it, the latency of random row and column queries and the peak
memory. Before every repeat the variant's files are dropped from the page cache with
posix_fadvise(DONTNEED) (no root needed, best effort), so the first pass over the query IDs
reads from disk; the second pass repeats the same IDs and finds them in the page cache.

The build time is only the conversion from the container, not a full pipeline run.

Usage:
  python3 bench_formats.py <data_dir> <work_dir> [--queries 1000] [--repeats 3] [--seed 42]
                           [--variants magrix,hdf5_chunk64k,...] [--csv results.csv] [--keep]
"""
import argparse
import csv
import json
import os
import random
import shutil
import sqlite3
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from magrix.api import Container  # noqa: E402

ARRAYS = ["csr.indptr", "csr.indices", "csr.data", "csc.indptr", "csc.indices", "csc.data"]
COPY_STEP = 1 << 23  # elements per write while building a variant, keeps memory bounded
SQLITE_BATCH = 1 << 20  # rows per executemany

# variant -> (format, chunk size in elements, option)
VARIANTS = {
    "magrix": ("magrix", None, None),
    "hdf5_contiguous": ("hdf5", None, None),
    "hdf5_chunk64k": ("hdf5", 65_536, None),
    "hdf5_chunk64k_gzip": ("hdf5", 65_536, "gzip"),
    "hdf5_chunk2M_gzip": ("hdf5", 2_000_000, "gzip"),
    "zarr_chunk64k": ("zarr", 65_536, None),
    "zarr_chunk64k_zstd": ("zarr", 65_536, "zstd"),
    "zarr_chunk2M_zstd": ("zarr", 2_000_000, "zstd"),
    "sqlite_rowid": ("sqlite", None, "rowid"),
    "sqlite_norowid": ("sqlite", None, "norowid"),
}
EXTENSIONS = {"hdf5": ".h5", "zarr": ".zarr", "sqlite": ".sqlite"}


def variant_path(name, work_dir, container_path):
    fmt = VARIANTS[name][0]
    if fmt == "magrix":
        return container_path
    return os.path.join(work_dir, name + EXTENSIONS[fmt])


def files_of(path):
    if os.path.isdir(path):
        return [os.path.join(d, f) for d, _, fs in os.walk(path) for f in fs]
    return [path]


def remove(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def evict_from_page_cache(path):
    """Ask the kernel to drop the cached pages of these files (fsync first, dirty pages stay)."""
    for p in files_of(path):
        fd = os.open(p, os.O_RDONLY)
        try:
            os.fsync(fd)
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def copy_ranges(n, chunk):
    # write in whole chunks so that no chunk is compressed twice
    step = COPY_STEP if not chunk else chunk * max(1, COPY_STEP // chunk)
    for s in range(0, n, step):
        yield s, min(n, s + step)


# ── building ────────────────────────────────────────────────────────────────

def build_hdf5(c, path, chunk, option):
    import h5py
    with h5py.File(path, "w") as f:
        for name in ARRAYS:
            src = c[name]
            kw = {}
            if chunk:
                kw["chunks"] = (min(chunk, len(src)),)
            if option == "gzip":
                kw.update(compression="gzip", compression_opts=4)
            ds = f.create_dataset(name.replace(".", "/"), shape=src.shape, dtype=src.dtype, **kw)
            for s, e in copy_ranges(len(src), chunk):
                ds[s:e] = src[s:e]


def build_zarr(c, path, chunk, option):
    import zarr
    from zarr.codecs import BloscCodec
    compressors = [BloscCodec(cname="zstd", clevel=5)] if option == "zstd" else None
    root = zarr.open_group(store=path, mode="w")
    for name in ARRAYS:
        group, arr_name = name.split(".")
        src = c[name]
        arr = root.require_group(group).create_array(
            name=arr_name, shape=src.shape, dtype=src.dtype,
            chunks=(min(chunk, len(src)),), compressors=compressors)
        for s, e in copy_ranges(len(src), chunk):
            arr[s:e] = np.asarray(src[s:e])


def build_sqlite(c, path, option):
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    if option == "rowid":
        conn.execute("CREATE TABLE paralog (mag_id INTEGER, family_id INTEGER, value INTEGER)")
    else:
        conn.execute("CREATE TABLE paralog (mag_id INTEGER, family_id INTEGER, value INTEGER, "
                     "PRIMARY KEY (mag_id, family_id)) WITHOUT ROWID")
    # insert in CSR order, i.e. sorted by (mag_id, family_id): appends to the primary key
    indptr = np.asarray(c["csr.indptr"], dtype=np.int64)
    indices, data = c["csr.indices"], c["csr.data"]
    with conn:
        for s, e in copy_ranges(len(indices), SQLITE_BATCH):
            mags = np.searchsorted(indptr, np.arange(s, e), side="right") - 1
            conn.executemany("INSERT INTO paralog VALUES (?, ?, ?)",
                             zip(mags.tolist(), indices[s:e].tolist(), data[s:e].tolist()))
    if option == "rowid":
        conn.execute("CREATE INDEX idx_mag ON paralog(mag_id)")
        conn.execute("CREATE INDEX idx_family ON paralog(family_id)")
    else:
        # an index of a WITHOUT ROWID table also holds the primary key, so this one
        # answers column queries without touching the table
        conn.execute("CREATE INDEX idx_family ON paralog(family_id, value)")
    conn.commit()
    conn.close()


def child_build(args):
    fmt, chunk, option = VARIANTS[args.variant]
    c = Container(args.container)
    remove(args.path)
    t0 = time.perf_counter()
    if fmt == "hdf5":
        build_hdf5(c, args.path, chunk, option)
    elif fmt == "zarr":
        build_zarr(c, args.path, chunk, option)
    else:
        build_sqlite(c, args.path, option)
    print(json.dumps({"build_s": time.perf_counter() - t0}))


# ── querying ────────────────────────────────────────────────────────────────

class MagrixBackend:
    def __init__(self, path):
        c = Container(path)
        self.arrays = [c[name] for name in ARRAYS]

    def row(self, i):
        indptr, indices, data = self.arrays[0:3]
        s, e = int(indptr[i]), int(indptr[i + 1])
        return np.array(indices[s:e]), np.array(data[s:e])

    def column(self, j):
        indptr, indices, data = self.arrays[3:6]
        s, e = int(indptr[j]), int(indptr[j + 1])
        return np.array(indices[s:e]), np.array(data[s:e])


class ArrayStoreBackend:
    """HDF5 and Zarr: same lookups, only the library that opens the arrays differs."""

    def __init__(self, path, fmt):
        if fmt == "hdf5":
            import h5py
            self.f = h5py.File(path, "r")
            self.arrays = [self.f[name.replace(".", "/")] for name in ARRAYS]
        else:
            import zarr
            g = zarr.open_group(store=path, mode="r")
            self.arrays = [g[name.split(".")[0]][name.split(".")[1]] for name in ARRAYS]

    def _lookup(self, indptr, indices, data, k):
        bounds = indptr[k:k + 2]  # one read for both bounds
        s, e = int(bounds[0]), int(bounds[1])
        return indices[s:e], data[s:e]

    def row(self, i):
        return self._lookup(*self.arrays[0:3], i)

    def column(self, j):
        return self._lookup(*self.arrays[3:6], j)


class SQLiteBackend:
    def __init__(self, path):
        self.conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        self.conn.execute("SELECT count(*) FROM sqlite_master").fetchone()  # reads the schema

    # results stay a list of tuples (no conversion to numpy), which favours SQLite slightly
    def row(self, i):
        return self.conn.execute("SELECT family_id, value FROM paralog WHERE mag_id = ?", (i,)).fetchall()

    def column(self, j):
        return self.conn.execute("SELECT mag_id, value FROM paralog WHERE family_id = ?", (j,)).fetchall()


def open_backend(variant, path):
    fmt = VARIANTS[variant][0]
    if fmt == "magrix":
        return MagrixBackend(path)
    if fmt == "sqlite":
        return SQLiteBackend(path)
    return ArrayStoreBackend(path, fmt)


def as_sorted_pairs(result):
    if isinstance(result, list):
        pairs = np.array(result, dtype=np.int64).reshape(-1, 2)
    else:
        pairs = np.column_stack([np.asarray(result[0], dtype=np.int64), np.asarray(result[1], dtype=np.int64)])
    return pairs[np.argsort(pairs[:, 0], kind="stable")]


def latency_stats(times):
    us = sorted(t * 1e6 for t in times)
    p95 = us[min(len(us) - 1, int(round(0.95 * (len(us) - 1))))]
    return st.median(us), p95, us[-1]


def child_query(args):
    ids = json.loads(Path(args.ids).read_text())
    t0 = time.perf_counter()
    backend = open_backend(args.variant, args.path)
    out = {"open_ms": (time.perf_counter() - t0) * 1000}

    for kind, fn, id_list in (("row", backend.row, ids["rows"]), ("col", backend.column, ids["cols"])):
        for pass_name in ("first", "second"):
            times = []
            for k in id_list:
                t = time.perf_counter()
                fn(k)
                times.append(time.perf_counter() - t)
            med, p95, mx = latency_stats(times)
            out.update({f"{kind}_{pass_name}_median_us": med, f"{kind}_{pass_name}_p95_us": p95,
                        f"{kind}_{pass_name}_max_us": mx})

    # same answers as the container? (outside the timing)
    ref = MagrixBackend(args.container)
    bad = 0
    for i in ids["rows"][:50]:
        bad += not np.array_equal(as_sorted_pairs(backend.row(i)), as_sorted_pairs(ref.row(i)))
    for j in ids["cols"][:50]:
        bad += not np.array_equal(as_sorted_pairs(backend.column(j)), as_sorted_pairs(ref.column(j)))
    out["mismatches"] = bad
    print(json.dumps(out))


# ── driver ──────────────────────────────────────────────────────────────────

def run_child(argv):
    """Run this script in a fresh process; returns its JSON result and its peak memory in MB."""
    proc = subprocess.Popen([sys.executable, os.path.abspath(__file__)] + argv, stdout=subprocess.PIPE, text=True)
    stdout = proc.stdout.read()
    proc.stdout.close()
    _, status, usage = os.wait4(proc.pid, 0)
    proc.returncode = os.waitstatus_to_exitcode(status)
    if proc.returncode != 0:
        sys.exit(f"child process failed: {' '.join(argv)}")
    return json.loads(stdout.strip().splitlines()[-1]), usage.ru_maxrss / 1024  # Linux: KB


def library_versions():
    versions = [f"Python {sys.version.split()[0]}", f"numpy {np.__version__}", f"SQLite {sqlite3.sqlite_version}"]
    try:
        import h5py
        versions.append(f"h5py {h5py.__version__} (HDF5 {h5py.version.hdf5_version})")
    except ImportError:
        versions.append("h5py missing")
    try:
        import zarr
        versions.append(f"zarr {zarr.__version__}")
    except ImportError:
        versions.append("zarr missing")
    return ", ".join(versions)


def size_on_disk(variant, path, container):
    if VARIANTS[variant][0] == "magrix":
        # only the six matrix arrays, members are not part of the comparison
        return sum(len(container[name]) * container[name].itemsize for name in ARRAYS), 1
    files = files_of(path)
    return sum(os.path.getsize(p) for p in files), len(files)


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("_build", "_query"):
        p = argparse.ArgumentParser()
        p.add_argument("mode")
        p.add_argument("variant")
        p.add_argument("path")
        p.add_argument("container")
        p.add_argument("--ids")
        a = p.parse_args()
        (child_build if a.mode == "_build" else child_query)(a)
        return

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("data_dir", help="directory with paralog.container (output of magrix build)")
    parser.add_argument("work_dir", help="where the other formats are written (needs a few times the container size)")
    parser.add_argument("--queries", type=int, default=1000, help="row and column queries per pass (default 1000)")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--variants", default=",".join(VARIANTS), help="comma-separated, default: all")
    parser.add_argument("--csv", help="write one line per variant and repeat to this CSV file")
    parser.add_argument("--keep", action="store_true", help="keep the built files (default: delete after measuring)")
    args = parser.parse_args()

    variants = args.variants.split(",")
    unknown = [v for v in variants if v not in VARIANTS]
    if unknown:
        parser.error(f"unknown variants {unknown}, choose from {list(VARIANTS)}")
    container_path = os.path.abspath(os.path.join(args.data_dir, "paralog.container"))
    c = Container(container_path)
    num_mags, num_families = len(c["csr.indptr"]) - 1, len(c["csc.indptr"]) - 1
    nnz = len(c["csr.indices"])
    print(f"Dataset: {args.data_dir}: {num_mags:,} MAGs, {num_families:,} gene families, {nnz:,} non-zero entries")
    print(f"Container file {os.path.getsize(container_path) / 1e9:.3f} GB (incl. members); "
          f"{args.queries} row + {args.queries} column queries per pass, {args.repeats} repeats")
    print(f"Versions: {library_versions()}\n")

    os.makedirs(args.work_dir, exist_ok=True)
    rng = random.Random(args.seed)
    ids_path = os.path.join(args.work_dir, "query_ids.json")
    Path(ids_path).write_text(json.dumps({
        "rows": [rng.randrange(num_mags) for _ in range(args.queries)],
        "cols": [rng.randrange(num_families) for _ in range(args.queries)],
    }))

    results = []
    for variant in variants:
        path = variant_path(variant, args.work_dir, container_path)
        build = {"build_s": None, "build_peak_mb": None}
        if VARIANTS[variant][0] != "magrix":
            print(f"[{variant}] building ...", flush=True)
            out, peak = run_child(["_build", variant, path, container_path])
            build = {"build_s": out["build_s"], "build_peak_mb": peak}
        size, num_files = size_on_disk(variant, path, c)
        for r in range(1, args.repeats + 1):
            evict_from_page_cache(path)
            out, peak = run_child(["_query", variant, path, container_path, "--ids", ids_path])
            results.append({"variant": variant, "repeat": r, "num_mags": num_mags, "num_families": num_families,
                            "nnz": nnz, "bytes": size, "files": num_files, **build, **out, "query_peak_mb": peak})
        last = results[-1]
        print(f"[{variant}] {size / 1e9:.3f} GB in {num_files} files, row first/second "
              f"{last['row_first_median_us']:.0f}/{last['row_second_median_us']:.0f} us, column "
              f"{last['col_first_median_us']:.0f}/{last['col_second_median_us']:.0f} us, "
              f"{last['mismatches']} mismatches", flush=True)
        if not args.keep and VARIANTS[variant][0] != "magrix":
            remove(path)

    # ── summary: median over the repeats ────────────────────────────────────
    def med(variant, key):
        values = [r[key] for r in results if r["variant"] == variant and r[key] is not None]
        return st.median(values) if values else None

    print("\nMedian over repeats; latencies are medians over the queries of one pass, in microseconds")
    header = (f"{'variant':20} {'GB':>7} {'files':>6} {'build s':>8} {'build MB':>9} {'open ms':>8} "
              f"{'row 1st':>9} {'row 2nd':>8} {'col 1st':>9} {'col 2nd':>8} {'query MB':>9} {'bad':>4}")
    print(header)
    print("-" * len(header))
    for v in variants:
        b, bm = med(v, "build_s"), med(v, "build_peak_mb")
        print(f"{v:20} {med(v, 'bytes') / 1e9:>7.3f} {med(v, 'files'):>6.0f} "
              f"{'-' if b is None else f'{b:.1f}':>8} {'-' if bm is None else f'{bm:.0f}':>9} "
              f"{med(v, 'open_ms'):>8.1f} {med(v, 'row_first_median_us'):>9.0f} {med(v, 'row_second_median_us'):>8.0f} "
              f"{med(v, 'col_first_median_us'):>9.0f} {med(v, 'col_second_median_us'):>8.0f} "
              f"{med(v, 'query_peak_mb'):>9.0f} {sum(r['mismatches'] for r in results if r['variant'] == v):>4}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]))
            w.writeheader()
            w.writerows(results)
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
