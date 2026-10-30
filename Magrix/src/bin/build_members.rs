// Build the "members" structure from mmseqs_long_gene_list.tsv: for each
// (family, MAG) cell, the exact list of individual member genes — the detail
// the matrix can't hold, since CSR/CSC only store one scalar value per cell.
//
// gene_list.tsv columns: representative_gene, MAG, member_gene (packed u64
// gene ids; the MAG column is unused here, mag_id is derived from member_gene).
// It's grouped by representative_gene (family), so — like build_csc — a single
// streaming pass with a one-family buffer is enough; no external sort needed.
//
// Within a family, genes are sorted by (mag_id, gene number). This is the same
// cell order as the CSC representation (family, then MAG ascending), so the
// genes of the cell at position i in a CSC column start right after the genes
// of positions 0..i-1 of that column. Two consequences keep the output small:
//   - only ONE pointer per family is stored; the start of a cell inside the
//     family is the family pointer plus the copy numbers of the preceding
//     cells in that CSC column (computed at query time)
//   - only the local gene number (lower 32 bits) is stored; the MAG part is
//     the MAG of the cell, which the query already knows
//
// Output (raw little-endian binary):
//   <out_prefix>.members.family_ptr    — u64[num_families + 1]
//   <out_prefix>.members.gene_numbers  — u32[total gene_list rows]
//
// Memory: only the genes of the current family are buffered (8 bytes each),
// so usage scales with the largest family, not with the input size.

// Imports
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

// Buffer sizes in MiB, overridable via MAGRIX_READ_BUF_MB / MAGRIX_WRITE_BUF_MB.
// 8 MiB each: no measurable speed difference between 1 and 256 MiB in the buffer sweep (results/).
const DEFAULT_READ_BUF_MB: usize = 8;
const DEFAULT_WRITE_BUF_MB: usize = 8;

// Buffer size in bytes: the environment variable `var` (in MiB) if set, else the default.
fn buffer_size(var: &str, default_mb: usize) -> usize {
    let mb = match env::var(var) {
        Ok(v) => v.parse::<usize>().ok().filter(|&n| n > 0)
            .unwrap_or_else(|| panic!("{} must be a positive integer (MiB), got {:?}", var, v)),
        Err(_) => default_mb,
    };
    mb * 1024 * 1024
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: build_members <gene_list.tsv> <out_prefix>");
        std::process::exit(1);
    }
    let tsv_path = &args[1];
    let out_prefix = &args[2];

    let read_buf = buffer_size("MAGRIX_READ_BUF_MB", DEFAULT_READ_BUF_MB);
    let write_buf = buffer_size("MAGRIX_WRITE_BUF_MB", DEFAULT_WRITE_BUF_MB);
    eprintln!("Buffers: read={} MiB write={} MiB", read_buf >> 20, write_buf >> 20);

    let infile = File::open(tsv_path)?;
    let reader = BufReader::with_capacity(read_buf, infile);

    let mut family_ptr_w = BufWriter::with_capacity(write_buf, File::create(format!("{}.members.family_ptr", out_prefix))?);
    let mut gene_numbers_w = BufWriter::with_capacity(write_buf, File::create(format!("{}.members.gene_numbers", out_prefix))?);

    // Packed gene ids (mag_id << 32 | gene number) of the family currently being read.
    let mut buffer: Vec<u64> = Vec::new();
    let mut current_rep: Option<u64> = None;
    let mut genes_written: u64 = 0; // running value of family_ptr
    let mut families: u64 = 0;
    let mut cells: u64 = 0; // (family, MAG) pairs, for cross-checking against the CSC nnz
    let mut max_gene_number: u32 = 0; // diagnostic: would gene numbers fit into u16?
    let mut rows_read: u64 = 0;

    // Flush one family: write its start pointer, then its gene numbers in CSC cell order.
    let mut flush_family = |buffer: &mut Vec<u64>,
                            family_ptr_w: &mut BufWriter<File>,
                            gene_numbers_w: &mut BufWriter<File>|
     -> io::Result<()> {
        if buffer.is_empty() {
            return Ok(());
        }
        // Sorting the packed ids sorts by mag_id first, then by gene number,
        // because mag_id sits in the upper 32 bits.
        buffer.sort_unstable();

        // The family starts where the previous family's genes ended
        family_ptr_w.write_all(&genes_written.to_le_bytes())?;
        families += 1;

        let mut prev_mag: Option<u32> = None;
        for &gene_id in buffer.iter() {
            // Split the packed id: upper 32 bits = mag_id, lower 32 bits = gene number
            let mag_id = (gene_id >> 32) as u32;
            let gene_number = (gene_id & 0xFFFF_FFFF) as u32;
            // A new mag_id within the family = a new (family, MAG) cell
            if prev_mag != Some(mag_id) {
                cells += 1;
                prev_mag = Some(mag_id);
            }
            max_gene_number = max_gene_number.max(gene_number);
            // Only the gene number is written; the mag_id is known from the cell at query time
            gene_numbers_w.write_all(&gene_number.to_le_bytes())?;
        }
        genes_written += buffer.len() as u64;
        buffer.clear();
        Ok(())
    };

    for line in reader.lines() {
        let line = line?;
        // Skip empty lines and the header line
        if line.is_empty() || line.starts_with("representative_gene") {
            continue;
        }
        let mut parts = line.split('\t');
        let rep_id: u64 = parts.next().unwrap().parse().expect("bad representative_gene");
        let _mag_col = parts.next().unwrap(); // unused: mag_id is part of member_gene
        let member_id: u64 = parts.next().unwrap().parse().expect("bad member_gene");
        rows_read += 1;

        // Cluster boundary: same invariant main.rs relies on (rows for the same
        // representative_gene are contiguous), NOT that rep_id itself increases —
        // the packed gene id has no reason to be monotonic across clusters.
        if current_rep != Some(rep_id) {
            flush_family(&mut buffer, &mut family_ptr_w, &mut gene_numbers_w)?;
            current_rep = Some(rep_id);
        }
        buffer.push(member_id);
    }

    // Flush the last family (the loop only flushes a family when the next one starts)
    flush_family(&mut buffer, &mut family_ptr_w, &mut gene_numbers_w)?;
    // Final pointer: total number of genes (closes the last family).
    family_ptr_w.write_all(&genes_written.to_le_bytes())?;

    family_ptr_w.flush()?;
    gene_numbers_w.flush()?;

    eprintln!(
        "members done. rows_read={} genes={} families={} cells={} max_gene_number={}",
        rows_read, genes_written, families, cells, max_gene_number
    );
    Ok(())
}
