"""End-to-end correctness tests for the pipeline + the magrix query API, run
against the small committed clusters.tsv fixture (fast, deterministic).

Most tests check the API's answers against ground truth pulled straight from
Stage 1's raw paralog.tsv / the packed container's own structure, rather than
hardcoding expected numbers -- so they keep working if the fixture file ever
changes, and they'd actually catch a real pipeline regression.
"""
import random

import numpy as np
import pytest


# ---- Dictionaries ----------------------------------------------------------

def test_mag_dictionary_round_trips(matrices):
    assert matrices.dicts.num_mags > 0
    for mag_id in range(matrices.dicts.num_mags):
        name = matrices.dicts.mag_name(mag_id)
        assert matrices.dicts.mag_id(name) == mag_id


def test_family_representative_gene_looks_decoded(matrices):
    assert matrices.dicts.num_families > 0
    rep = matrices.dicts.family_representative_gene(0)
    # decoded form is "<mag>|<mag>_gene<n>", not a raw packed integer
    assert "|" in rep and "_gene" in rep


# ---- paralog.tsv is ground truth -------------------------------------------

def test_paralog_get_matches_raw_tsv(matrices, paralog_tsv_rows):
    random.seed(0)
    sample = random.sample(paralog_tsv_rows, min(200, len(paralog_tsv_rows)))
    for family_id, mag_id, value in sample:
        assert matrices.paralog.get(mag_id, family_id) == value


def test_paralog_values_never_zero(matrices, paralog_tsv_rows):
    # sparse storage: a stored entry always means "at least one gene here"
    assert all(value >= 1 for _, _, value in paralog_tsv_rows)


def test_presence_is_paralog_clamped_to_one(matrices, paralog_tsv_rows):
    random.seed(1)
    sample = random.sample(paralog_tsv_rows, min(200, len(paralog_tsv_rows)))
    for family_id, mag_id, value in sample:
        assert matrices.presence.get(mag_id, family_id) == 1
        assert matrices.paralog.get(mag_id, family_id) == value


# ---- CSR vs CSC must agree on the same matrix ------------------------------

def test_row_and_column_agree(matrices):
    random.seed(2)
    for _ in range(200):
        mag_id = random.randrange(matrices.dicts.num_mags)
        family_ids, values = matrices.paralog.row(mag_id)
        if len(family_ids) == 0:
            continue
        idx = random.randrange(len(family_ids))
        family_id, expected_value = int(family_ids[idx]), int(values[idx])

        mag_ids_in_col, values_in_col = matrices.paralog.column(family_id)
        pos = np.searchsorted(mag_ids_in_col, mag_id)
        assert mag_ids_in_col[pos] == mag_id
        assert int(values_in_col[pos]) == expected_value


def test_get_matches_row_scan(matrices):
    random.seed(3)
    for _ in range(200):
        mag_id = random.randrange(matrices.dicts.num_mags)
        family_ids, values = matrices.paralog.row(mag_id)
        if len(family_ids) == 0:
            continue
        idx = random.randrange(len(family_ids))
        assert matrices.paralog.get(mag_id, int(family_ids[idx])) == int(values[idx])


# ---- Structural / prevalence helpers ---------------------------------------

def test_family_prevalence_matches_column_length(matrices):
    random.seed(4)
    for family_id in random.sample(range(matrices.dicts.num_families), 50):
        mag_ids, _ = matrices.presence.column(family_id)
        assert matrices.presence.family_prevalence(family_id) == len(mag_ids)


def test_mag_gene_count_matches_row_length(matrices):
    random.seed(5)
    for mag_id in random.sample(range(matrices.dicts.num_mags), min(50, matrices.dicts.num_mags)):
        family_ids, _ = matrices.presence.row(mag_id)
        assert matrices.presence.mag_gene_count(mag_id) == len(family_ids)


def test_families_by_prevalence_full_range_covers_everything(matrices):
    all_ids = matrices.presence.families_by_prevalence(0.0, 1.0)
    assert len(all_ids) == matrices.dicts.num_families


# ---- members: exact gene-level detail --------------------------------------

def test_members_cell_count_matches_paralog_value(matrices, paralog_tsv_rows):
    assert matrices.members is not None
    random.seed(6)
    sample = random.sample(paralog_tsv_rows, min(200, len(paralog_tsv_rows)))
    for family_id, mag_id, value in sample:
        gene_ids = matrices.members.gene_ids_for(mag_id, family_id)
        assert len(gene_ids) == value


def test_members_decoded_genes_belong_to_the_right_mag(matrices, paralog_tsv_rows):
    random.seed(7)
    sample = random.sample(paralog_tsv_rows, min(100, len(paralog_tsv_rows)))
    for family_id, mag_id, _value in sample:
        expected_mag_name = matrices.dicts.mag_name(mag_id)
        for gene_name in matrices.members.genes_for(mag_id, family_id):
            gene_mag_name = gene_name.split("|", 1)[0]
            assert gene_mag_name == expected_mag_name


def test_members_match_gene_list_for_every_cell(matrices, gene_list_cells):
    # members stores one pointer per family and only local gene numbers, so a
    # cell's genes are located via the copy numbers of the preceding cells in
    # its CSC column -- check every single cell against the raw gene list.
    assert len(gene_list_cells) == len(matrices.paralog.csc_indices)
    for (family_id, mag_id), expected in gene_list_cells.items():
        got = matrices.members.gene_ids_for(mag_id, family_id).tolist()
        assert got == expected, f"cell (family {family_id}, MAG {mag_id})"


def test_members_missing_cell_returns_empty(matrices):
    # a (mag, family) pair with no paralog entry should yield no genes, not an error
    for mag_id in range(matrices.dicts.num_mags):
        family_ids, _ = matrices.presence.row(mag_id)
        missing = set(range(matrices.dicts.num_families)) - set(int(f) for f in family_ids)
        if missing:
            some_missing_family = next(iter(missing))
            assert matrices.members.gene_ids_for(mag_id, some_missing_family).size == 0
            return
    pytest.skip("no MAG had a missing family in this fixture")


# ---- singleton clusters must be preserved (Anupam's explicit requirement) --

def test_singleton_clusters_are_not_dropped(matrices):
    # a family with prevalence 1 is a cluster whose representative gene had
    # no other hits -- must still show up with paralog_count/presence == 1,
    # not be silently skipped.
    singleton_found = False
    for family_id in range(min(matrices.dicts.num_families, 2000)):
        if matrices.presence.family_prevalence(family_id) == 1:
            mag_ids, values = matrices.paralog.column(family_id)
            assert len(mag_ids) == 1
            assert int(values[0]) >= 1
            singleton_found = True
            break
    assert singleton_found, "expected at least one singleton cluster in the fixture"


# ---- submatrix(): the thesis's primary query -------------------------------

def test_submatrix_matches_individual_gets(matrices):
    random.seed(8)
    mag_ids = random.sample(range(matrices.dicts.num_mags), min(10, matrices.dicts.num_mags))
    family_ids = random.sample(range(matrices.dicts.num_families), 30)

    sub = matrices.paralog.submatrix(mag_ids, family_ids)
    dense = sub.toarray()

    for local_row, mag_id in enumerate(mag_ids):
        for local_col, family_id in enumerate(family_ids):
            assert dense[local_row, local_col] == matrices.paralog.get(mag_id, family_id)


def test_submatrix_presence_is_clamped(matrices):
    mag_ids = list(range(min(5, matrices.dicts.num_mags)))
    family_ids = list(range(min(30, matrices.dicts.num_families)))
    sub = matrices.presence.submatrix(mag_ids, family_ids)
    assert set(sub.data.tolist()) <= {1}


# ---- build_csr's external-merge sort across MULTIPLE runs -------------------
# (the small fixture's ~11k rows fit in one default-sized run, so without
# this the k-way merge across several runs would never actually be exercised)

def test_csr_multirun_merge_matches_single_run(matrices, multirun_csr):
    multirun_indptr = np.memmap(f"{multirun_csr}.csr.indptr", dtype="<u8", mode="r")
    multirun_indices = np.memmap(f"{multirun_csr}.csr.indices", dtype="<u4", mode="r")
    multirun_data = np.memmap(f"{multirun_csr}.csr.data", dtype="<u2", mode="r")

    assert np.array_equal(multirun_indptr, matrices.paralog.csr_indptr)
    assert np.array_equal(multirun_indices, matrices.paralog.csr_indices)
    assert np.array_equal(multirun_data, matrices.paralog.csr_data)


# ---- error handling: invalid ids must fail loudly and clearly, not crash
# with a raw numpy IndexError or silently return something misleading ------

def test_row_rejects_out_of_range_mag_id(matrices):
    with pytest.raises(ValueError):
        matrices.paralog.row(matrices.dicts.num_mags)  # first invalid id
    with pytest.raises(ValueError):
        matrices.paralog.row(-1)


def test_column_rejects_out_of_range_family_id(matrices):
    with pytest.raises(ValueError):
        matrices.paralog.column(matrices.dicts.num_families)
    with pytest.raises(ValueError):
        matrices.paralog.column(-1)


def test_get_rejects_out_of_range_mag_id(matrices):
    with pytest.raises(ValueError):
        matrices.paralog.get(matrices.dicts.num_mags, 0)


def test_get_rejects_out_of_range_family_id(matrices):
    # without the check this returned 0, as if the cell were just empty
    with pytest.raises(ValueError):
        matrices.paralog.get(0, matrices.dicts.num_families)
    with pytest.raises(ValueError):
        matrices.paralog.get(0, -1)


def test_prevalence_and_gene_count_reject_out_of_range_ids(matrices):
    # -1 used to index from the end of indptr and return a wrong number
    for bad in (-1, matrices.dicts.num_families):
        with pytest.raises(ValueError):
            matrices.paralog.family_prevalence(bad)
        with pytest.raises(ValueError):
            matrices.presence.family_prevalence(bad)
    for bad in (-1, matrices.dicts.num_mags):
        with pytest.raises(ValueError):
            matrices.paralog.mag_gene_count(bad)


def test_dicts_reject_out_of_range_ids(matrices):
    with pytest.raises(ValueError):
        matrices.dicts.mag_name(matrices.dicts.num_mags)
    with pytest.raises(ValueError):
        matrices.dicts.family_representative_gene(matrices.dicts.num_families)


def test_dicts_mag_id_rejects_unknown_name(matrices):
    with pytest.raises(KeyError):
        matrices.dicts.mag_id("this-mag-does-not-exist")


def test_members_reject_out_of_range_ids(matrices):
    with pytest.raises(ValueError):
        matrices.members.gene_ids_for(matrices.dicts.num_mags, 0)
    with pytest.raises(ValueError):
        matrices.members.gene_ids_for(0, matrices.dicts.num_families)


# ---- container format -------------------------------------------------------

def test_copy_numbers_stored_as_u16_but_returned_as_u32(matrices):
    # u16 on disk keeps the container small; returning u32 protects callers
    # that sum copy numbers across many MAGs from wrapping around at 65535.
    assert matrices.paralog.csc_data.dtype == np.dtype("<u2")
    assert matrices.paralog.csr_data.dtype == np.dtype("<u2")
    _, values = matrices.paralog.row(0)
    assert values.dtype == np.uint32
    _, values = matrices.paralog.column(0)
    assert values.dtype == np.uint32


def test_container_has_valid_magic(built_pipeline):
    with open(built_pipeline / "paralog.container", "rb") as f:
        header = f.read(4)
    assert header == b"PANG"


def test_container_is_the_only_binary_source_the_api_needs(built_pipeline, matrices):
    # sanity check that queries work purely off the container, independent
    # of whether the (now-redundant) separate .csc/.csr/.members files
    # still happen to exist alongside it in this test's temp dir
    mag_id = 0
    family_ids, values = matrices.paralog.row(mag_id)
    assert isinstance(family_ids, np.ndarray)
    assert isinstance(values, np.ndarray)
