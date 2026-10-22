"""Builds the full pipeline once with `magrix build` against the small
committed clusters.tsv fixture, into a temp directory, and hands the result
to every test as a session-scoped fixture -- so the whole test session only
pays the pipeline build cost once.

It runs with --keep-intermediates, so mmseqs_long_paralog.tsv and
mmseqs_long_gene_list.tsv stay around and tests can check the API's answers
directly against that raw ground truth, instead of only checking internal
self-consistency.
"""
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from magrix import PanGenomeMatrices  # noqa: E402


def _run(cmd, **kwargs):
    result = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    assert result.returncode == 0, f"command failed: {cmd}\nstdout={result.stdout}\nstderr={result.stderr}"
    return result


@pytest.fixture(scope="session")
def built_pipeline(tmp_path_factory):
    # Always (re)build: a no-op when up to date, and it guarantees the tests
    # never run against stale binaries from an older version of the code.
    _run(["cargo", "build", "--release"], cwd=REPO_ROOT)
    bin_dir = REPO_ROOT / "target" / "release"

    out_dir = tmp_path_factory.mktemp("pipeline_out")
    _run([str(bin_dir / "magrix"), "build", str(REPO_ROOT / "clusters.tsv"), str(out_dir),
          "--keep-intermediates"])
    return out_dir


@pytest.fixture(scope="session")
def matrices(built_pipeline):
    return PanGenomeMatrices(str(built_pipeline))


@pytest.fixture(scope="session")
def multirun_csr(built_pipeline):
    """Rebuilds CSR for the same paralog.tsv but with a tiny rows_per_run, so
    Phase 1 is forced to produce several runs instead of the usual one --
    otherwise the external-merge (k-way merge across multiple runs) code
    path is never actually exercised, since the small fixture's ~11k rows
    all fit in a single default-sized run (20,000,000)."""
    bin_dir = REPO_ROOT / "target" / "release"
    num_mags = sum(1 for _ in open(built_pipeline / "mag_dictionary.tsv")) - 1
    prefix = built_pipeline / "paralog_multirun"
    result = _run([
        str(bin_dir / "build_csr"), str(built_pipeline / "mmseqs_long_paralog.tsv"),
        str(num_mags), str(prefix), "500",  # rows_per_run=500 -> forces multiple runs
    ])
    num_runs = int(result.stderr.split(" in ")[-1].split(" runs")[0])
    assert num_runs > 1, f"fixture didn't actually force multiple runs (got {num_runs}) -- test wouldn't cover the merge path"
    return prefix


@pytest.fixture(scope="session")
def gene_list_cells(built_pipeline):
    """Ground truth for members: {(family_id, mag_id): sorted packed gene ids},
    straight from Stage 1's gene_list.tsv. family_id is recovered through the
    family dictionary (representative gene -> family_id)."""
    rep_to_family = {}
    with open(built_pipeline / "family_dictionary.tsv") as f:
        next(f)
        for line in f:
            family_id, rep = line.rstrip("\n").split("\t")
            rep_to_family[int(rep)] = int(family_id)

    cells = {}
    with open(built_pipeline / "mmseqs_long_gene_list.tsv") as f:
        next(f)
        for line in f:
            rep, mag_id, member = line.rstrip("\n").split("\t")
            key = (rep_to_family[int(rep)], int(mag_id))
            cells.setdefault(key, []).append(int(member))
    return {k: sorted(v) for k, v in cells.items()}


@pytest.fixture(scope="session")
def paralog_tsv_rows(built_pipeline):
    """Ground-truth rows straight from Stage 1's paralog.tsv, kept around
    (not deleted) specifically so tests can check the API against them."""
    rows = []
    with open(built_pipeline / "mmseqs_long_paralog.tsv") as f:
        next(f)  # header
        for line in f:
            _rep, family_id, mag_id, value = line.rstrip("\n").split("\t")
            rows.append((int(family_id), int(mag_id), int(value)))
    return rows
