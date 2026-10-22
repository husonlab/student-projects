#!/usr/bin/env python3
"""
Latency of a running magrix-serve over HTTP, e.g. from a laptop through an SSH tunnel.

Three endpoints: /info (no matrix data, so it measures the pure network and server
overhead), a row and a column query with random IDs. Each is measured twice:
  new      a new connection per request (like separate curl calls)
  reused   one connection kept open for all requests (HTTP keep-alive)

Uses only the Python standard library, so it runs on any machine.

Usage: python3 bench_http.py [--url http://localhost:8765] [--queries 50] [--seed 1] [--csv out.csv]
"""
import argparse
import csv
import http.client
import json
import random
import statistics as st
import time
from urllib.parse import urlparse


def timed_request(conn, path):
    t0 = time.perf_counter()
    conn.request("GET", path)
    resp = conn.getresponse()
    body = resp.read()
    elapsed = time.perf_counter() - t0
    if resp.status != 200:
        raise RuntimeError(f"GET {path} returned HTTP {resp.status}: {body[:200]!r}")
    return elapsed, len(body)


def measure(host, port, paths, reuse):
    times, sizes = [], []
    conn = http.client.HTTPConnection(host, port, timeout=60) if reuse else None
    for path in paths:
        if not reuse:
            conn = http.client.HTTPConnection(host, port, timeout=60)
        t, size = timed_request(conn, path)
        times.append(t)
        sizes.append(size)
        if not reuse:
            conn.close()
    if reuse:
        conn.close()
    return times, sizes


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--url", default="http://localhost:8765")
    parser.add_argument("--queries", type=int, default=50, help="requests per endpoint and mode (default 50)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--csv", help="also write the results to this CSV file")
    args = parser.parse_args()

    url = urlparse(args.url)
    host, port = url.hostname, url.port or 80

    conn = http.client.HTTPConnection(host, port, timeout=60)
    conn.request("GET", "/info")
    info = json.loads(conn.getresponse().read())
    conn.close()
    print(f"Server at {args.url}: {info['num_mags']:,} MAGs, {info['num_families']:,} gene families\n")

    rng = random.Random(args.seed)
    n = args.queries
    workloads = [
        ("/info (no data)", ["/info"] * n),
        ("row", [f"/paralog/row/{rng.randrange(info['num_mags'])}" for _ in range(n)]),
        ("column", [f"/paralog/column/{rng.randrange(info['num_families'])}" for _ in range(n)]),
    ]

    rows = []
    header = f"{'endpoint':16} {'connection':10} {'n':>4} {'median ms':>10} {'p95 ms':>9} {'max ms':>9} {'avg bytes':>10}"
    print(header)
    print("-" * len(header))
    for name, paths in workloads:
        for mode in ("new", "reused"):
            times, sizes = measure(host, port, paths, reuse=(mode == "reused"))
            ms = sorted(t * 1000 for t in times)
            p95 = ms[min(len(ms) - 1, int(round(0.95 * (len(ms) - 1))))]
            row = (name, mode, len(ms), st.median(ms), p95, ms[-1], st.mean(sizes))
            rows.append(row)
            print(f"{name:16} {mode:10} {row[2]:>4} {row[3]:>10.1f} {row[4]:>9.1f} {row[5]:>9.1f} {row[6]:>10.0f}")

    print("\nNote: the row/column IDs of the 'new' and 'reused' runs are the same, so the second run "
          "finds the data in the server's page cache; /info isolates the network overhead.")
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["endpoint", "connection", "n", "median_ms", "p95_ms", "max_ms", "avg_bytes"])
            w.writerows(rows)
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
