# MAGrix

**MAGrix** turns the gene clustering output of [MMseqs2](https://github.com/soedinglab/MMseqs2)
into a compact **MAG × gene family matrix** and lets you query it from Python without ever
loading the matrix into memory.

It is made for pan-genomes of hundreds of thousands of metagenome-assembled genomes (MAGs),
tens of millions of gene families and billions of non-zero entries: building the matrix needs
less than 1 GB of memory regardless of the input size, and a single row or column query takes
a few microseconds once the data is in the page cache.

- **Out-of-core construction.** No stage assumes the data fits into memory.
- **Lossless.** Stores copy numbers (not only presence/absence) *and* the individual member
  genes of every cell. Singleton clusters are kept.
- **Fast in both directions.** "Which families does this MAG have" and "which MAGs carry
  this family" both read only the requested part of the matrix.
- **Three output files.** One binary container plus two plain-text dictionaries.

---

## Contents

- [How the parts fit together](#how-the-parts-fit-together)
- [Quick start](#quick-start)
- [Installation](#installation)
- [Building a dataset](#building-a-dataset)
- [Querying from Python](#querying-from-python)
- [Querying over HTTP](#querying-over-http)
- [Reference](#reference): input format, output files, container format, configuration,
  tests, limits, repository layout

---

## How the parts fit together

```
clusters.tsv ──magrix build──►  paralog.container      ──►  Python API          ──►  HTTP server
 (MMseqs2)       (Rust)          mag_dictionary.tsv         from magrix import …      magrix-serve
                                 family_dictionary.tsv      (magrix/api.py)          (magrix/server.py)
```

| Part | What it does | Who needs it |
|---|---|---|
| `magrix build` (Rust) | builds the three output files from a clustering | whoever builds a dataset, once per dataset |
| Python API (`magrix/api.py`) | **the** query interface: `PanGenomeMatrices` | everyone who works with a built dataset |
| HTTP server (`magrix/server.py`) | a thin layer that calls the same Python API and returns JSON | only if you want to query a dataset from another machine |

There is only one query API. The HTTP server does no work of its own; every endpoint calls
the corresponding Python method. Querying needs no Rust at all.

## Quick start

From the repository folder:

```bash
cargo install --path .                       # build tool (Rust)
magrix build clusters.tsv output/my_dataset  # build a dataset
pip install .                                # query library (Python)
```

```python
from magrix import PanGenomeMatrices

m = PanGenomeMatrices("output/my_dataset")
families, copies = m.paralog.row(m.dicts.mag_id("Mag568575"))
```

## Installation

### Build tool (Rust)

Needs a Rust toolchain (`cargo`, edition 2021). No external crates are used.

```bash
cargo install --path .           # installs magrix and the stage programs into ~/.cargo/bin
```

or build inside the repository and call `./target/release/magrix`:

```bash
cargo build --release
```

### Query library (Python)

Needs Python ≥ 3.9; NumPy and SciPy are installed as dependencies.

```bash
pip install .                        # query library
pip install ".[server]"              # plus the HTTP server (FastAPI, uvicorn)
pip install ".[test]"                # plus pytest
pip install ".[bench]"               # plus h5py, Zarr and matplotlib for bench_formats.py and the plots (Python ≥ 3.11)
```

**Without installing**, e.g. when the repository sits in a shared folder you can only read:
`pip install` writes build files into the source folder, so instead tell Python where the
package is. Either in the shell

```bash
export PYTHONPATH=/path/to/repository    # the repository folder, which contains the package folder magrix/
```

or at the top of your script

```python
import sys
sys.path.insert(0, "/path/to/repository")
from magrix import PanGenomeMatrices
```

NumPy and SciPy then have to be in your environment already.

Tested with Rust 1.95, Python 3.13, NumPy 2.4 and SciPy 1.18; used with Python 3.11 on the
cluster.

## Building a dataset

```bash
magrix build clusters.tsv output/my_dataset
```

This runs the four stages in order, deletes intermediate files as soon as they are no longer
needed, and prints the run time of every stage and the size of the output files.

| Option | Meaning |
|---|---|
| `--rows-per-run N` | entries sorted in memory per run in stage 2b (default 20,000,000, ~240 MB) |
| `--keep-intermediates` | keep the intermediate tables and array files |
| `-h`, `--help` | help |

The input has to be an MMseqs2 cluster table with gene names of the form
`<MAG>|<MAG>_gene<N>` (see [Input format](#input-format)).

### On a SLURM cluster

Wrap the command in a small batch script:

```bash
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=4000m
magrix build datasets/clusters.tsv output/my_dataset
```

Memory stays below 1 GB for an input of about 1.4 billion rows, so a small `--mem` is enough;
`--time` depends on the input size and the file system (roughly one hour for 1.4 billion rows).

### Stage by stage

`magrix build` calls one program per stage, which can also be run on their own:

```bash
OUT=output/my_dataset
pangenome_matrices clusters.tsv $OUT                                 # 1: preprocessing, IDs, dictionaries
NUM_FAMILIES=$(($(wc -l < $OUT/family_dictionary.tsv) - 1))
NUM_MAGS=$(($(wc -l < $OUT/mag_dictionary.tsv) - 1))
build_csc $OUT/mmseqs_long_paralog.tsv $NUM_FAMILIES $OUT/paralog    # 2a: column-wise matrix (CSC)
build_csr $OUT/mmseqs_long_paralog.tsv $NUM_MAGS     $OUT/paralog    # 2b: row-wise matrix (CSR), external merge sort
build_members $OUT/mmseqs_long_gene_list.tsv $OUT/paralog            # 3: member genes of every cell
pack_container $OUT/paralog                                          # 4: pack everything into one file
rm $OUT/mmseqs_long_*.tsv $OUT/paralog.csc.* $OUT/paralog.csr.* $OUT/paralog.members.*
```

## Querying from Python

```python
from magrix import PanGenomeMatrices

m = PanGenomeMatrices("output/my_dataset")      # takes a few seconds: loads the dictionaries

mag = m.dicts.mag_id("Mag568575")               # MAG name -> id
families, copies = m.paralog.row(mag)           # gene families of this MAG + copy numbers
mags, copies = m.paralog.column(42)             # MAGs that carry family 42 + copy numbers
n = m.paralog.get(mag, 42)                      # one cell (0 if empty)
genes = m.members.genes_for(mag, 42)            # member genes of that cell, as names
sub = m.paralog.submatrix([0, 1, 2], [5, 6, 7]) # scipy.sparse COO matrix
```

Opening a dataset reads the 4 KB container header and the two dictionaries. The matrix itself
is memory-mapped: only the parts a query touches are read from disk, so memory use does not
grow with the size of the matrix.

| Call | Returns |
|---|---|
| `paralog.row(mag_id)` | `(family_ids, copy_numbers)` as NumPy arrays |
| `paralog.column(family_id)` | `(mag_ids, copy_numbers)` as NumPy arrays |
| `paralog.get(mag_id, family_id)` | copy number as `int`, 0 if the cell is empty |
| `paralog.submatrix(mag_ids, family_ids)` | `scipy.sparse.coo_matrix` |
| `paralog.family_prevalence(family_id)` | number of MAGs that carry the family |
| `paralog.mag_gene_count(mag_id)` | number of distinct gene families of the MAG |
| `paralog.families_by_prevalence(min_frac, max_frac)` | IDs of families found in that fraction of MAGs (scans all families) |
| `presence.*` | the same methods with values 0/1 |
| `members.gene_ids_for(mag_id, family_id)` | encoded gene IDs of one cell (`uint64`) |
| `members.genes_for(mag_id, family_id)` | gene names of one cell |
| `dicts.mag_id(name)` / `dicts.mag_name(mag_id)` | MAG name ↔ ID |
| `dicts.family_representative_gene(family_id)` | name of the family's representative gene |
| `dicts.num_mags`, `dicts.num_families` | size of the matrix |

Invalid IDs raise `ValueError`, unknown MAG names raise `KeyError`.

### Copy numbers: `uint16` on disk, `uint32` in results

Copy numbers are stored as 16-bit integers (0–65,535) to save space; the largest copy number
in the real data is 351. All methods above return them as `uint32`, so sums over many MAGs
are safe.

If you read the raw arrays directly (`m.paralog.csr_data`, `m.paralog.csc_data`), you get
`uint16`, and **NumPy wraps around silently** when the result of element-wise arithmetic
exceeds 65,535:

```python
a = m.paralog.csc_data[:2]          # say [300, 65300], dtype uint16
a + a                               # [600, 65064]   <- second value wrong, no warning
a.astype(np.uint32) + a             # [600, 130600]  <- correct
```

`.sum()` is safe (NumPy uses a wider type for it), but element-wise `+`, `*`, `np.add.at`
or accumulating into a `uint16` array are not. Convert with `.astype(np.uint32)` first, or use
`row()` / `column()`, which already do that.

## Querying over HTTP

With the `server` extra installed, `magrix-serve` offers the same queries over HTTP, so a
dataset can stay on a server while it is queried from elsewhere. It is a proof of concept:
no authentication, one process.

### Start, check, stop

On the machine that holds the data, inside `tmux` so it keeps running after you log out:

```bash
tmux new -s api                                  # (if it exists: tmux attach -t api)
magrix-serve output/my_dataset --port 8765       # ready when "Uvicorn running on ..." appears
# detach with Ctrl+B, then D
```

| Task | Command |
|---|---|
| Is a server running on this machine? | `pgrep -u $USER -af magrix-serve` (no output = no server) |
| Go back to it | `tmux attach -t api` |
| Stop it | in the tmux session: `Ctrl+C`, then `exit` |

`tmux` sessions and servers exist only on the machine where they were started; on a cluster
with several login nodes, note which one you are on (`hostname`).

### Connect from your own computer

The server only listens on `127.0.0.1` (`--host` changes that), so nothing is reachable from
outside. Open an SSH tunnel **to the host the server runs on** and keep that terminal open:

```bash
ssh -L 8765:localhost:8765 <user>@<host>
```

If the cluster has several login nodes behind one address that can't be reached individually,
jump over the shared address to the node the server runs on (`hostname` on the server shows it):

```bash
ssh -J <user>@<shared-login-address> -L 8765:localhost:8765 <user>@<login-node>
```

Then, in another terminal on your computer:

```bash
curl http://localhost:8765/info
curl http://localhost:8765/paralog/row/42
```

Interactive documentation of all endpoints is at `http://localhost:8765/docs`.

| Method | Path |
|---|---|
| GET | `/info` |
| GET | `/mag_name/{mag_id}`, `/mag_id/{name}`, `/family_representative_gene/{family_id}` |
| GET | `/{matrix}/row/{mag_id}`, `/{matrix}/column/{family_id}`, `/{matrix}/get/{mag_id}/{family_id}` |
| GET | `/{matrix}/family_prevalence/{family_id}`, `/{matrix}/mag_gene_count/{mag_id}` |
| GET | `/{matrix}/families_by_prevalence?min_frac=…&max_frac=…` |
| POST | `/{matrix}/submatrix` with body `{"mag_ids": [...], "family_ids": [...]}` |
| GET | `/members/{mag_id}/{family_id}` |

`{matrix}` is `paralog` (copy numbers) or `presence` (0/1). Responses are JSON; invalid IDs
return HTTP 400. The app can also be run directly with
`MAGRIX_DIR=output/my_dataset uvicorn magrix.server:app`.

Through an SSH tunnel a request takes 20–80 ms (median), almost all of it network: the query
itself takes microseconds. Keeping one connection open for many requests roughly halves that
(`bench_http.py` measures it).

---

## Reference

### Input format

A two-column, tab-separated table as written by MMseqs2 (`createtsv`): the representative gene
of a cluster, then one member gene. Surrounding double quotes are removed.

```
"Mag568047|Mag568047_gene2842"	"Mag568047|Mag568047_gene2842"
"Mag568047|Mag568047_gene2842"	"Mag568575|Mag568575_gene3424"
"Mag568047|Mag568047_gene2842"	"Mag565394|Mag565394_gene4306"
```

1. **Gene identifiers follow `<MAG>|<MAG>_gene<N>`.** The MAG and the gene number are taken
   from the identifier. A violation aborts the run with an error message.
2. **All rows of a cluster are contiguous.** MMseqs2 writes its output this way. This is *not*
   checked; interleaved clusters would be split into separate families.

Empty lines and lines without a tab are skipped.

### Output files

| File | Content |
|---|---|
| `paralog.container` | binary container with the CSC and CSR arrays and the member genes |
| `mag_dictionary.tsv` | `mag_id` → MAG name, e.g. `1  Mag568575` |
| `family_dictionary.tsv` | `family_id` → representative gene, as an encoded gene ID |

All IDs start at 0 and are assigned in order of first appearance in the input. A gene is
encoded as one 64-bit integer, `mag_id · 2³² + gene_number`; the MAG dictionary turns it back
into the original identifier:

```
4294970720 = 1 · 2³² + 3424   →   MAG 1 = Mag568575   →   Mag568575|Mag568575_gene3424
```

### Container format

A container starts with a 4,096-byte header, followed by the arrays. Each array starts at a
multiple of 4,096 bytes, so it can be memory-mapped directly at its offset. All values are
little-endian.

| Bytes | Content |
|---|---|
| 0–3 | magic number `PANG` |
| 4–7 | number of arrays (`u32`) |
| 8– | one 64-byte record per array: name (32 bytes, zero-padded), type code (1 byte: 0 = `u32`, 1 = `u64`, 2 = `u16`), offset (`u64`), number of elements (`u64`), 15 reserved bytes |

| Array | Type | Length | Content |
|---|---|---|---|
| `csc.indptr` | u64 | families + 1 | start of each family's entries |
| `csc.indices` | u32 | nnz | MAG IDs, ascending within each family |
| `csc.data` | u16 | nnz | copy numbers |
| `csr.indptr` | u64 | MAGs + 1 | start of each MAG's entries |
| `csr.indices` | u32 | nnz | family IDs, ascending within each MAG |
| `csr.data` | u16 | nnz | copy numbers |
| `members.family_ptr` | u64 | families + 1 | start of each family's genes in `gene_numbers` |
| `members.gene_numbers` | u32 | number of genes | local gene numbers, ordered by family, then MAG |

The genes of a cell are found through the CSC column of its family: they start at
`family_ptr[family]` plus the copy numbers of the cells before it in that column. Their MAG is
the MAG of the cell, so only the local gene number is stored. Looking up a cell near the end of
a very large family therefore sums more copy numbers; on the largest family (83,781 MAGs) this
takes 64 µs instead of the usual ~17 µs.

### Configuration

Buffer sizes are fixed defaults that can be overridden with environment variables (values in
MiB); `magrix build` passes them on to the stages:

| Variable | Default | Used by |
|---|---|---|
| `MAGRIX_READ_BUF_MB` | 8 | stages 1–3, input buffer |
| `MAGRIX_WRITE_BUF_MB` | 8 | stages 1–3, per output file |
| `MAGRIX_MERGE_BUF_MB` | 4 | stage 2b, per run during the merge |
| `MAGRIX_COPY_BUF_MB` | 16 | stage 4 |

```bash
MAGRIX_READ_BUF_MB=32 magrix build clusters.tsv output/my_dataset --rows-per-run 5000000
```

Measurements of these settings are in `results/buffer_sweep.csv`. Read and write buffers between 1 and 256 MiB made no measurable difference in run
time while memory grew with the buffer size, hence the 8 MiB defaults. For `--rows-per-run`,
memory is lowest when the sort buffer of phase 1 and the merge buffers of phase 2 are about
equal, at roughly √(nnz · merge buffer / 12 bytes); for 1.35 billion entries this is about
22 million, close to the default.

### Tests

From the repository folder:

```bash
cargo test --release            # unit tests of the Rust programs
python3 -m pytest tests/        # builds the pipeline on clusters.tsv and checks the API
```

The Python tests rebuild the Rust programs (a no-op if they are up to date), run
`magrix build` on the small `clusters.tsv` fixture and check the query results against the
intermediate tables, including the member genes of every single cell.

Scripts for real datasets:

| Script | Measures |
|---|---|
| `test_real_data.py <data_dir>` | open time, latency of every query type (cold and warm), worst-case member query, consistency check |
| `bench_http.py --url http://localhost:8765` | latency of a running `magrix-serve` over HTTP |
| `bench_formats.py <data_dir> <work_dir>` | the same data in HDF5, Zarr and SQLite: size, build time, query latency (needs `pip install ".[bench]"`) |

### Limits and assumptions

- At most 2³² MAGs, gene families and genes per MAG.
- Copy numbers must not exceed 65,535 (stored as 16-bit integers; larger values abort the run).
- The input must be grouped by cluster (see [Input format](#input-format)).
- Opening a dataset loads both dictionaries into memory: for 96 million families this takes
  about 7 seconds and 870 MB.
- All stages are single-threaded.

### Repository layout

```
Cargo.toml, pyproject.toml
src/main.rs              stage 1
src/bin/                 magrix (whole pipeline) and stages 2–4
magrix/                  Python package: api.py (queries), server.py (HTTP)
tests/                   Python tests
clusters.tsv             small test fixture
test_real_data.py, bench_http.py, bench_formats.py   measurement scripts (see Tests)
results/                 measured data as CSV, and the scripts that plot it
```
