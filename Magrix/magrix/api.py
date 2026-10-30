"""
Query API over the CSR/CSC/members arrays built by `magrix build` and packed
into one paralog.container file.

Uses numpy.memmap (at an offset into the container) for indptr/indices/data,
so opening a matrix costs almost nothing — only the specific rows/columns
actually queried get paged in from disk by the OS. This is the read-side
counterpart to the out-of-core Rust builders: memory usage stays bounded by
the size of a query, not by the size of the whole matrix. (What does take time
and memory when opening a dataset are the two dictionaries, which are loaded
completely: about 7 s and 870 MB for 96 million families.)

Structure of this file:
    Container           reads the header of paralog.container, one memmap per array
    Dictionaries        MAG names and family representatives (the two TSV files)
    SparseMatrix        row/column/cell queries on the CSR + CSC arrays (copy numbers)
    PresenceView        the same queries with 0/1 values, computed from SparseMatrix
    MembersView         the individual genes behind one cell
    PanGenomeMatrices   opens everything; the only class a user needs

Usage:
    from magrix import PanGenomeMatrices
    m = PanGenomeMatrices("output/my_dataset")
    m.paralog.row(mag_id)                        # family ids + copy numbers in one MAG
    m.presence.row(mag_id)                       # family ids in one MAG
    m.presence.column(family_id)                 # MAGs carrying one family
    m.presence.get(mag_id, family_id)             # single cell
    m.presence.submatrix(mag_ids, family_ids)     # scipy COO submatrix
    m.presence.family_prevalence(family_id)       # how many MAGs carry it
    m.presence.families_by_prevalence(0.99, 1.0)  # e.g. "core genome" ids
    m.members.genes_for(mag_id, family_id)        # exact member genes in one cell
    m.dicts.mag_name(mag_id) / m.dicts.mag_id(name)
    m.dicts.family_representative_gene(family_id)
"""
import numpy as np
import scipy.sparse as sp


def _decode_gene_id(gene_id, mag_names):
    """Packed gene id (mag_id << 32 | gene number) -> original name "<MAG>|<MAG>_gene<N>"."""
    mag_id = gene_id >> 32
    local = gene_id & 0xFFFFFFFF
    name = mag_names[mag_id]
    return f"{name}|{name}_gene{local}"


class Container:
    """Reads a single packed .container file (built by pack_container) and
    exposes each of its arrays as a numpy.memmap view directly into the
    file -- no data is read until an array is actually accessed.
    See src/bin/pack_container.rs for the exact header/layout format."""

    # dtype code in the header -> little-endian numpy type (u32, u64, u16)
    _DTYPES = {0: "<u4", 1: "<u8", 2: "<u2"}

    def __init__(self, path):
        # The whole header is one 4096-byte block at the start of the file
        with open(path, "rb") as f:
            header = f.read(4096)
        assert header[0:4] == b"PANG", f"not a pangenome container file: {path}"
        num_arrays = int.from_bytes(header[4:8], "little")
        self._arrays = {}
        # One 64-byte record per array: name, dtype code, byte offset, number of elements
        for i in range(num_arrays):
            base = 8 + i * 64
            name = header[base:base + 32].rstrip(b"\x00").decode("utf-8")
            dtype_code = header[base + 32]
            offset = int.from_bytes(header[base + 33:base + 41], "little")
            count = int.from_bytes(header[base + 41:base + 49], "little")
            dtype = self._DTYPES[dtype_code]
            # Maps the array's part of the file; nothing is read until it is indexed
            self._arrays[name] = np.memmap(path, dtype=dtype, mode="r", offset=offset, shape=(count,))

    def __getitem__(self, name):
        return self._arrays[name]

    def __contains__(self, name):
        return name in self._arrays


class Dictionaries:
    """MAG name <-> id, and family_id -> representative_gene lookups.
    Both dictionaries are loaded completely into memory when the dataset is opened."""

    def __init__(self, out_dir):
        # mag_names[i] is the name of MAG i; the dict is the reverse direction
        self.mag_names = []
        self._mag_name_to_id = {}
        with open(f"{out_dir}/mag_dictionary.tsv") as f:
            next(f)  # header
            for line in f:
                mid, name = line.rstrip("\n").split("\t")
                mid = int(mid)
                assert mid == len(self.mag_names), "mag_dictionary.tsv must be id-ordered"
                self.mag_names.append(name)
                self._mag_name_to_id[name] = mid

        # One-time load of the family_id -> representative_gene column, as packed u64 ids
        # (8 bytes per family, ~770 MB for 96 million families; this is most of the open time).
        # Names are only decoded on request, see family_representative_gene.
        self.family_rep_ids = np.loadtxt(
            f"{out_dir}/family_dictionary.tsv",
            delimiter="\t", skiprows=1, usecols=1, dtype=np.uint64,
        )

    @property
    def num_mags(self):
        return len(self.mag_names)

    @property
    def num_families(self):
        return len(self.family_rep_ids)

    def mag_name(self, mag_id):
        if not (0 <= mag_id < self.num_mags):
            raise ValueError(f"mag_id {mag_id} out of range (0..{self.num_mags})")
        return self.mag_names[mag_id]

    def mag_id(self, name):
        if name not in self._mag_name_to_id:
            raise KeyError(f"unknown MAG name: {name!r}")
        return self._mag_name_to_id[name]

    def family_representative_gene(self, family_id):
        if not (0 <= family_id < self.num_families):
            raise ValueError(f"family_id {family_id} out of range (0..{self.num_families})")
        return _decode_gene_id(int(self.family_rep_ids[family_id]), self.mag_names)


class SparseMatrix:
    """One CSR+CSC matrix pair (the paralog matrix), backed by memmapped
    views into a single packed container file (see Container above).

    Copy numbers are stored as u16 but returned as u32, so that callers who
    add up values across many MAGs can't silently overflow at 65535."""

    def __init__(self, container, dicts):
        self.dicts = dicts
        self.csr_indptr = container["csr.indptr"]
        self.csr_indices = container["csr.indices"]
        self.csr_data = container["csr.data"]
        self.csc_indptr = container["csc.indptr"]
        self.csc_indices = container["csc.indices"]
        self.csc_data = container["csc.data"]

    # ---- Tier 1: core access ----

    def row(self, mag_id):
        """Family ids (+ values) present in one MAG. Only reads that row's byte range."""
        if not (0 <= mag_id < self.dicts.num_mags):
            raise ValueError(f"mag_id {mag_id} out of range (0..{self.dicts.num_mags})")
        # indptr[i] .. indptr[i+1] is the range of row i in indices/data
        start, end = int(self.csr_indptr[mag_id]), int(self.csr_indptr[mag_id + 1])
        # np.array copies the slice out of the memmap, so the result doesn't depend on the open file
        return np.array(self.csr_indices[start:end]), np.array(self.csr_data[start:end], dtype=np.uint32)

    def column(self, family_id):
        """MAG ids (+ values) that carry one gene family. Only reads that column's byte range."""
        if not (0 <= family_id < self.dicts.num_families):
            raise ValueError(f"family_id {family_id} out of range (0..{self.dicts.num_families})")
        start, end = int(self.csc_indptr[family_id]), int(self.csc_indptr[family_id + 1])
        return np.array(self.csc_indices[start:end]), np.array(self.csc_data[start:end], dtype=np.uint32)

    def get(self, mag_id, family_id):
        """Single-cell lookup: binary search within the (sorted) row."""
        # row() checks mag_id; without this check an invalid family_id would just return 0
        if not (0 <= family_id < self.dicts.num_families):
            raise ValueError(f"family_id {family_id} out of range (0..{self.dicts.num_families})")
        # Reads the whole row of the MAG; family ids within a CSR row are ascending
        family_ids, values = self.row(mag_id)
        pos = np.searchsorted(family_ids, family_id)
        if pos < len(family_ids) and family_ids[pos] == family_id:
            return int(values[pos])
        return 0

    def submatrix(self, mag_ids, family_ids):
        """MAG subset x gene-family subset -> matching entries as a small scipy
        COO matrix. Cost scales with the size of the selected rows, not with
        the full matrix. Row i / column j of the result correspond to
        mag_ids[i] / family_ids[j]."""
        mag_ids = list(mag_ids)
        family_ids = list(family_ids)
        wanted = np.sort(np.asarray(family_ids))
        # global family id -> column index in the result
        fam_to_local = {int(f): i for i, f in enumerate(family_ids)}

        rows, cols, data = [], [], []
        # Read each selected row once and keep only the requested families
        for local_row, mag_id in enumerate(mag_ids):
            fam_in_row, vals_in_row = self.row(int(mag_id))
            mask = np.isin(fam_in_row, wanted)
            for f, v in zip(fam_in_row[mask], vals_in_row[mask]):
                rows.append(local_row)
                cols.append(fam_to_local[int(f)])
                data.append(int(v))

        return sp.coo_matrix((data, (rows, cols)), shape=(len(mag_ids), len(family_ids)))

    # ---- Tier 2: prevalence — just indptr differences, no indices/data are read ----

    def family_prevalence(self, family_id):
        """How many MAGs carry this gene family."""
        # Checked explicitly: a negative id would otherwise index from the end and return a wrong value
        if not (0 <= family_id < self.dicts.num_families):
            raise ValueError(f"family_id {family_id} out of range (0..{self.dicts.num_families})")
        # length of the family's column = number of MAGs with an entry
        return int(self.csc_indptr[family_id + 1] - self.csc_indptr[family_id])

    def mag_gene_count(self, mag_id):
        """How many distinct gene families this MAG has (not genes: for those,
        sum the copy numbers of row(mag_id))."""
        if not (0 <= mag_id < self.dicts.num_mags):
            raise ValueError(f"mag_id {mag_id} out of range (0..{self.dicts.num_mags})")
        return int(self.csr_indptr[mag_id + 1] - self.csr_indptr[mag_id])

    def families_by_prevalence(self, min_frac=0.0, max_frac=1.0):
        """Family ids whose fraction of carrying MAGs falls in [min_frac, max_frac].
        e.g. families_by_prevalence(0.99, 1.0) ~= core genome,
             families_by_prevalence(0.0, 1/num_mags) ~= unique/cloud genes."""
        # Prevalence of every family at once: reads the whole csc.indptr
        # (8 bytes per family, ~770 MB for 96 million families), but never indices/data
        counts = np.diff(self.csc_indptr)
        frac = counts / self.dicts.num_mags
        return np.nonzero((frac >= min_frac) & (frac <= max_frac))[0]


class PresenceView:
    """Presence matrix derived from paralog's structure, with no storage of its
    own: a cell is stored in the paralog matrix exactly when the MAG has at
    least one gene of the family, so "presence" is just "does paralog have an
    entry here" — any stored cell implies presence=1, since only nonzero
    copy numbers are stored."""

    def __init__(self, paralog):
        self._paralog = paralog
        self.dicts = paralog.dicts

    def row(self, mag_id):
        family_ids, values = self._paralog.row(mag_id)
        return family_ids, np.ones_like(values, dtype=np.uint8)

    def column(self, family_id):
        mag_ids, values = self._paralog.column(family_id)
        return mag_ids, np.ones_like(values, dtype=np.uint8)

    def get(self, mag_id, family_id):
        return 1 if self._paralog.get(mag_id, family_id) > 0 else 0

    def submatrix(self, mag_ids, family_ids):
        sub = self._paralog.submatrix(mag_ids, family_ids)
        sub.data = np.ones_like(sub.data)
        return sub

    # Structural queries are identical to paralog's — same indices/indptr.
    def family_prevalence(self, family_id):
        return self._paralog.family_prevalence(family_id)

    def mag_gene_count(self, mag_id):
        return self._paralog.mag_gene_count(mag_id)

    def families_by_prevalence(self, min_frac=0.0, max_frac=1.0):
        return self._paralog.families_by_prevalence(min_frac, max_frac)


class MembersView:
    """Individual member genes for each (MAG, family) cell — the detail the
    matrix can't hold, since CSR/CSC only store one scalar value per cell.
    Built by build_members from gene_list.tsv.

    Genes are stored per family in the same cell order as the CSC
    representation (MAG-ascending within each family), with only one pointer
    per family. The genes of the cell at position i in a CSC column start at
    the family pointer plus the copy numbers of cells 0..i-1 of that column.
    Only the local gene number is stored; the MAG part of the gene id is the
    MAG of the cell and is added back when the ids are returned."""

    def __init__(self, container, paralog, dicts):
        self.paralog = paralog
        self.dicts = dicts
        self.family_ptr = container["members.family_ptr"]
        self.gene_numbers = container["members.gene_numbers"]

    def _gene_range(self, mag_id, family_id):
        """Start/end position of one cell's genes in gene_numbers, or None if the cell is empty."""
        if not (0 <= mag_id < self.dicts.num_mags):
            raise ValueError(f"mag_id {mag_id} out of range (0..{self.dicts.num_mags})")
        if not (0 <= family_id < self.dicts.num_families):
            raise ValueError(f"family_id {family_id} out of range (0..{self.dicts.num_families})")
        # 1. Find the cell in the family's CSC column: the MAG ids there are ascending, so binary search
        start, end = int(self.paralog.csc_indptr[family_id]), int(self.paralog.csc_indptr[family_id + 1])
        mags_in_col = self.paralog.csc_indices[start:end]
        pos = int(np.searchsorted(mags_in_col, mag_id))
        if pos >= len(mags_in_col) or mags_in_col[pos] != mag_id:
            return None  # the MAG has no gene of this family
        # 2. Copy numbers of the cells up to and including this one (= their number of genes)
        counts = self.paralog.csc_data[start:start + pos + 1]
        # 3. The cell's genes start after the genes of all earlier cells of the family.
        #    Summed as u64: the stored u16 values could overflow otherwise.
        first = int(self.family_ptr[family_id]) + int(counts[:pos].sum(dtype=np.uint64))
        return first, first + int(counts[pos])

    def gene_ids_for(self, mag_id, family_id):
        """Packed gene ids (u64, mag_id << 32 | gene number) for one (MAG, family) cell."""
        cell = self._gene_range(mag_id, family_id)
        if cell is None:
            return np.array([], dtype=np.uint64)
        s, e = cell
        numbers = np.array(self.gene_numbers[s:e], dtype=np.uint64)
        # Put the MAG part back in front of the stored gene numbers -> the same packed ids as in stage 1
        return (np.uint64(mag_id) << np.uint64(32)) | numbers

    def genes_for(self, mag_id, family_id):
        """Same, decoded to readable gene name strings."""
        return [_decode_gene_id(int(g), self.dicts.mag_names) for g in self.gene_ids_for(mag_id, family_id)]


class PanGenomeMatrices:
    """Top-level entry point: loads dictionaries once, exposes presence + paralog matrices
    and the member genes. paralog, presence and members all share the same memmapped arrays."""

    def __init__(self, out_dir):
        self.dicts = Dictionaries(out_dir)  # the slow part of opening (see module docstring)
        container = Container(f"{out_dir}/paralog.container")  # only reads the 4 KB header
        self.paralog = SparseMatrix(container, self.dicts)
        self.presence = PresenceView(self.paralog)
        self.members = MembersView(container, self.paralog, self.dicts)
