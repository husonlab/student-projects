"""
Thin HTTP layer over magrix.PanGenomeMatrices, so the data can stay on the
server while a client (e.g. a laptop, over an SSH tunnel) queries it over the
network -- the library itself only works locally (numpy.memmap needs the
container file on the same machine).

Run on the machine that holds the data:
    magrix-serve output/my_dataset --port 8765
(equivalent: MAGRIX_DIR=output/my_dataset python3 -m uvicorn magrix.server:app --port 8765)

From your laptop, tunnel to it over SSH (data never leaves the server, only
the query results do):
    ssh -L 8765:localhost:8765 <user>@<host>
    curl http://localhost:8765/paralog/row/42

Whether this sits behind a reverse proxy or is reached directly through the
tunnel is independent of this file -- it just serves HTTP on one port.

Every endpoint calls one method of PanGenomeMatrices and returns its result
as JSON; there is no query logic in this file.
"""
import argparse
import os
from contextlib import asynccontextmanager
from typing import List

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from .api import PanGenomeMatrices

m: PanGenomeMatrices = None  # loaded once at startup, see lifespan below


# Runs once when the server starts (before the first request): opening the dataset
# takes a few seconds because of the dictionaries, so it is done once, not per request.
@asynccontextmanager
async def lifespan(_app):
    global m
    m = PanGenomeMatrices(os.environ.get("MAGRIX_DIR", "output"))
    yield


app = FastAPI(title="MAGrix API", lifespan=lifespan)


def _matrix(name: str):
    """The {matrix} part of the URL -> the matching matrix object."""
    if name == "paralog":
        return m.paralog
    if name == "presence":
        return m.presence
    raise HTTPException(status_code=400, detail=f"unknown matrix '{name}', use 'paralog' or 'presence'")


def _guard(fn, *args):
    """Runs an API call, translating the library's ValueError/KeyError (bad
    ids) into proper HTTP 400s instead of a raw 500 crash."""
    try:
        return fn(*args)
    except (ValueError, KeyError) as e:
        raise HTTPException(status_code=400, detail=str(e))


# ---- dictionaries -----------------------------------------------------------

@app.get("/info")
def info():
    return {"num_mags": m.dicts.num_mags, "num_families": m.dicts.num_families}


@app.get("/mag_name/{mag_id}")
def mag_name(mag_id: int):
    return {"mag_id": mag_id, "name": _guard(m.dicts.mag_name, mag_id)}


@app.get("/mag_id/{name}")
def mag_id(name: str):
    return {"name": name, "mag_id": _guard(m.dicts.mag_id, name)}


@app.get("/family_representative_gene/{family_id}")
def family_representative_gene(family_id: int):
    return {"family_id": family_id, "representative_gene": _guard(m.dicts.family_representative_gene, family_id)}


# ---- row / column / get, for either matrix ----------------------------------
# NumPy arrays are not JSON-serialisable, hence the .tolist() calls.

@app.get("/{matrix}/row/{mag_id}")
def row(matrix: str, mag_id: int):
    mtx = _matrix(matrix)
    family_ids, values = _guard(mtx.row, mag_id)
    return {"mag_id": mag_id, "family_ids": family_ids.tolist(), "values": values.tolist()}


@app.get("/{matrix}/column/{family_id}")
def column(matrix: str, family_id: int):
    mtx = _matrix(matrix)
    mag_ids, values = _guard(mtx.column, family_id)
    return {"family_id": family_id, "mag_ids": mag_ids.tolist(), "values": values.tolist()}


@app.get("/{matrix}/get/{mag_id}/{family_id}")
def get_cell(matrix: str, mag_id: int, family_id: int):
    mtx = _matrix(matrix)
    return {"mag_id": mag_id, "family_id": family_id, "value": _guard(mtx.get, mag_id, family_id)}


@app.get("/{matrix}/family_prevalence/{family_id}")
def family_prevalence(matrix: str, family_id: int):
    mtx = _matrix(matrix)
    return {"family_id": family_id, "prevalence": _guard(mtx.family_prevalence, family_id)}


@app.get("/{matrix}/mag_gene_count/{mag_id}")
def mag_gene_count(matrix: str, mag_id: int):
    mtx = _matrix(matrix)
    return {"mag_id": mag_id, "gene_count": _guard(mtx.mag_gene_count, mag_id)}


@app.get("/{matrix}/families_by_prevalence")
def families_by_prevalence(matrix: str, min_frac: float = 0.0, max_frac: float = 1.0):
    mtx = _matrix(matrix)
    ids = _guard(mtx.families_by_prevalence, min_frac, max_frac)
    return {"family_ids": ids.tolist()}


# Body of POST /{matrix}/submatrix; FastAPI checks that both are lists of integers
class SubmatrixRequest(BaseModel):
    mag_ids: List[int]
    family_ids: List[int]


@app.post("/{matrix}/submatrix")
def submatrix(matrix: str, req: SubmatrixRequest):
    mtx = _matrix(matrix)
    sub = _guard(mtx.submatrix, req.mag_ids, req.family_ids)
    coo = sub.tocoo()
    # Sparse result as coordinate lists: entry k is at (rows[k], cols[k]) with value data[k],
    # where row/col index into the requested mag_ids/family_ids
    return {
        "shape": list(sub.shape),
        "mag_ids": req.mag_ids,
        "family_ids": req.family_ids,
        "rows": coo.row.tolist(),
        "cols": coo.col.tolist(),
        "data": coo.data.tolist(),
    }


# ---- members: exact gene ids/names in one cell -------------------------------

@app.get("/members/{mag_id}/{family_id}")
def members(mag_id: int, family_id: int):
    genes = _guard(m.members.genes_for, mag_id, family_id)
    return {"mag_id": mag_id, "family_id": family_id, "genes": genes}


def main():
    """Entry point of the `magrix-serve` command."""
    parser = argparse.ArgumentParser(
        prog="magrix-serve",
        description="Serve a MAGrix output directory over HTTP.",
    )
    parser.add_argument("data_dir", help="directory with paralog.container and the two dictionaries")
    parser.add_argument("--host", default="127.0.0.1", help="interface to listen on (default: 127.0.0.1, local only)")
    parser.add_argument("--port", type=int, default=8000, help="port to listen on (default: 8000)")
    args = parser.parse_args()

    if not os.path.isfile(os.path.join(args.data_dir, "paralog.container")):
        parser.error(f"no paralog.container in {args.data_dir}")
    # lifespan() reads the data directory from this variable
    os.environ["MAGRIX_DIR"] = args.data_dir

    # imported here so that the query library works without the server extra installed
    import uvicorn

    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
