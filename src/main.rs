// Stage 1 of the pipeline (preprocessing and ID encoding).
// The Program takes a tab-delimited file of MMseqs2 cluster output (representative_gene \t member_gene) and produces a set of files suitable for constructing a large sparse matrix of gene-family copy numbers across MAGs
// The input has to be grouped by representative_gene, so that all members of a cluster are contiguous in the input (MMseqs2 writes it this way; it does not have to be sorted)
// Every cluster becomes one gene family; MAGs and families get sequential u32 ids, genes are packed into one u64 (see encode_gene_id)
// The input is read once, line by line; only the current cluster and the two dictionaries are kept in memory
// generates (in out_dir):
//   mmseqs_long_gene_list.tsv : one row per input line (representative_gene, mag_id, member_gene), input for build_members
//   mmseqs_long_paralog.tsv   : one row per (family, MAG) cell with the copy number, input for build_csc and build_csr
//   mag_dictionary.tsv        : mag_id -> MAG name
//   family_dictionary.tsv     : family_id -> packed id of the representative gene
// (no presence table: presence is derived from paralog at query time, paralog_count > 0)

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

// Buffer sizes for reading and writing operations, in MiB.
// Can be overridden via the environment variables MAGRIX_READ_BUF_MB / MAGRIX_WRITE_BUF_MB.
// 8 MiB each: the buffer sweep (results/) showed no measurable speed difference
// between 1 and 256 MiB, while memory grows with the buffer size.
const DEFAULT_READ_BUF_MB: usize = 8;   // read buffer for the input file
const DEFAULT_WRITE_BUF_MB: usize = 8;  // write buffer per output file

// Buffer size in bytes: the environment variable `var` (in MiB) if set, else the default.
fn buffer_size(var: &str, default_mb: usize) -> usize {
    let mb = match env::var(var) {
        Ok(v) => v.parse::<usize>().ok().filter(|&n| n > 0)
            .unwrap_or_else(|| panic!("{} must be a positive integer (MiB), got {:?}", var, v)),
        Err(_) => default_mb,
    };
    mb * 1024 * 1024
}

// Remove surrounding quotes from a string, if present
#[inline(always)]
fn unquote(s: &str) -> &str {
    let s = s.trim();
    if s.starts_with('"') && s.ends_with('"') && s.len() >= 2 {
        &s[1..s.len() - 1]
    } 
    else {
        s
    }
}

// Extract the MAG identifier from a gene identifier (everything before the first "|" )
#[inline(always)]
fn extract_mag(gene: &str) -> &str {
    match gene.find('|') {
        Some(pos) => &gene[..pos],
        None => gene,
    }
}

// Assign a u32 id for a MAG name, in order of the first appearance. If the MAG has been seen before, return its existing id
// Hashmap stores the mapping name -> id, Vector stores the names in order of their assigned ids.
#[inline(always)]
fn mag_id_for(mag: &str, dict: &mut HashMap<String, u32>, names: &mut Vec<String>) -> u32 {
    // Check if the MAG name already has an assigned id in the dictionary
    if let Some(&id) = dict.get(mag) {
        id
    } 
    // If new MAG name, assign a new id based on the current length of the names vector.
    else {
        let id = names.len() as u32;
        dict.insert(mag.to_string(), id);
        names.push(mag.to_string());
        id
    }
}

// Pack a gene id "<MAG>|<MAG>_gene<N>" into a single u64 as mag_id << 32 | N
// A Dictionary is not needed for individual genes: the combination of mag_id and local number uniquely identifies the gene
// The original string can be reconstructed using the MAG dictionary alone.
// Program panics if the name does not follow the assumed naming convention
#[inline(always)]
fn encode_gene_id(gene: &str, mag: &str, mag_id: u32) -> u64 {
    let pipe = gene
        .find('|')
        .unwrap_or_else(|| panic!("gene id missing '|': {}", gene));
    let after = &gene[pipe + 1..];
    let rest = after
        // Cuts off the MAG prefix from the gene suffix, leaving "_gene<N>"
        .strip_prefix(mag)
        .unwrap_or_else(|| panic!("gene suffix does not start with MAG '{}': {}", mag, gene));
    let local_str = rest
        // Cuts off everything except the local number after "_gene"
        .strip_prefix("_gene")
        .unwrap_or_else(|| panic!("gene id missing '_gene' suffix: {}", gene));
    let local: u32 = local_str
        .parse()
        .unwrap_or_else(|_| panic!("gene local number not an integer: {}", gene));
    // Pack the MAG id and local number into a single u64
    // The MAG id is shifted left by 32 bits and combined with the local number
    ((mag_id as u64) << 32) | (local as u64)
}


fn main() -> io::Result<()> {
    // Collect command-line arguments for input and output paths
    // 0: program name, 1: input file path, 2: output directory
    let args: Vec<String> = env::args().collect();
    let input_path = if args.len() > 1 {
        args[1].clone()
    } 
    else {
        // Default input file name if not provided
        "clusters.tsv".to_string()
    };
    // Output directory is optional
    let out_dir = if args.len() > 2 {
        args[2].clone()
    } 
    else {
        ".".to_string()
    };

    // Print input and output paths to the standard error for logging purposes
    eprintln!("Input : {}", input_path);
    eprintln!("Output: {}/", out_dir);
    let read_buf = buffer_size("MAGRIX_READ_BUF_MB", DEFAULT_READ_BUF_MB);
    let write_buf = buffer_size("MAGRIX_WRITE_BUF_MB", DEFAULT_WRITE_BUF_MB);
    eprintln!("Buffers: read={} MiB write={} MiB", read_buf >> 20, write_buf >> 20);

    // Open the input file
    let infile = File::open(&input_path)
        // If the file cannot be opened, print an error message and exit the program
        .unwrap_or_else(|e| { eprintln!("Cannot open '{}': {}", input_path, e); std::process::exit(1); });
    // Create a buffered reader with a specified capacity for efficient reading
    let reader = BufReader::with_capacity(read_buf, infile);

    let mut gl = BufWriter::with_capacity(
        write_buf,
        File::create(Path::new(&out_dir).join("mmseqs_long_gene_list.tsv"))?,
    );
    let mut pa = BufWriter::with_capacity(
        write_buf,
        File::create(Path::new(&out_dir).join("mmseqs_long_paralog.tsv"))?,
    );

    // Header lines for the output files
    // representative_gene / member_gene are packed u64 ids: (mag_id << 32 | local_number), MAG is the mag_id
    // MAGs and Families are assigned sequential u32 ids in order of first appearance, and stored in mag_dictionary.tsv and family_dictionary.tsv respectively
    writeln!(gl, "representative_gene\tMAG\tmember_gene")?;
    writeln!(pa, "representative_gene\tfamily_id\tMAG\tparalog_count")?;

    // MAG dictionary: mag_name -> mag_id, and a vector of names in order of their assigned ids.
    let mut mag_dict: HashMap<String, u32> = HashMap::new();
    let mut mag_names: Vec<String> = Vec::new();

    // Gene-family dictionary: family_id -> representative_gene packed id.
    // The family_id is assigned sequentially in order of first appearance of a representative_gene.
    // 8 bytes per family: this vector is the largest memory user of stage 1 (~770 MB for 96 million families)
    let mut family_reps: Vec<u64> = Vec::new();

    // Per MAG within the current cluster: member count so far (paralog count).
    // Only holds the MAGs of one cluster, it is cleared at every cluster boundary.
    let mut mag_counts: HashMap<u32, u32> = HashMap::new();
    let mut current_rep: Option<u64> = None;
    let mut current_family_id: u32 = 0;

    let mut rows_read: u64 = 0;
    let mut rows_written: u64 = 0;
    // Timer
    let start = Instant::now();

    // Flush the current cluster to the paralog file, and reset the mag_counts for the next cluster
    // Gets the information about the current cluster's representative gene and its family id, along with the counts of member genes per MAG
    // Writes one row per MAG (one matrix cell) to the paralog file, and clears the mag_counts for the next cluster.
    // The MAGs come out in HashMap order, i.e. not sorted; build_csc sorts them per family.
    let flush = |rep_id: u64,
                 family_id: u32,
                 mag_counts: &mut HashMap<u32, u32>,
                 pa: &mut BufWriter<File>,
                 rows_written: &mut u64| -> io::Result<()> {
        // For each MAG in the current cluster, write a row to the paralog file
        for (mag_id, cnt) in mag_counts.iter() {
            writeln!(pa, "{}\t{}\t{}\t{}", rep_id, family_id, mag_id, cnt)?;
            *rows_written += 1;
        }
        mag_counts.clear();
        Ok(())
    };

    // Read the input file line by line, processing each representative and member gene pair
    // all rows for the same representative_gene must be contiguous in the input
    // if they are not, a cluster is split into several families (no panic will occur)
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() { continue; }

        // Split the line into representative and member gene identifiers, unquote them if necessary
        // Lines without a tab are skipped
        let mut parts = line.splitn(2, '\t');
        let rep    = match parts.next() { Some(s) => unquote(s), None => continue };
        let member = match parts.next() { Some(s) => unquote(s), None => continue };

        rows_read += 1;

        // Extract the MAG identifier from the representative gene, get its assigned id, and encode the representative gene into a packed u64 id
        let rep_mag_str = extract_mag(rep);
        let rep_mag_id = mag_id_for(rep_mag_str, &mut mag_dict, &mut mag_names);
        let rep_id = encode_gene_id(rep, rep_mag_str, rep_mag_id);

        // If the current representative gene is different from the previous one, a new cluster starts:
        // flush the previous cluster's data to the output files and prepare for the new cluster
        if current_rep != Some(rep_id) {
            if let Some(prev_rep_id) = current_rep {
                // Write the paralog counts for the previous cluster to the output file
                flush(prev_rep_id, current_family_id, &mut mag_counts, &mut pa, &mut rows_written)?;
            }
            // Update the current representative gene and assign a new family_id for the new cluster
            current_rep = Some(rep_id);
            current_family_id = family_reps.len() as u32;
            family_reps.push(rep_id);
        }

        // Extract the MAG identifier from the member gene, get its assigned id, and encode the member gene into a packed u64 id
        let member_mag_str = extract_mag(member);
        let member_mag_id = mag_id_for(member_mag_str, &mut mag_dict, &mut mag_names);
        let member_id = encode_gene_id(member, member_mag_str, member_mag_id);

        // Write gene_list row immediately
        writeln!(gl, "{}\t{}\t{}", rep_id, member_mag_id, member_id)?;
        rows_written += 1;

        // Increment the paralog counter for this MAG
        *mag_counts.entry(member_mag_id).or_insert(0) += 1;

        // Progress every 10M rows
        if rows_read % 10_000_000 == 0 {
            let elapsed = start.elapsed().as_secs_f64();
            eprintln!(
                "  {:>12} rows read | {:>12} rows written | {:.1}s elapsed | {:.1}M rows/s",
                rows_read,
                rows_written,
                elapsed,
                rows_read as f64 / elapsed / 1_000_000.0
            );
        }
    }

    // Flush the last cluster
    if let Some(rep_id) = current_rep {
        flush(rep_id, current_family_id, &mut mag_counts, &mut pa, &mut rows_written)?;
    }

    // Force all buffers to disk
    gl.flush()?;
    pa.flush()?;

    // MAG dictionary: write the mapping of mag_id to mag_name to a TSV file
    let mut dict_file = BufWriter::new(
        File::create(Path::new(&out_dir).join("mag_dictionary.tsv"))?,
    );
    writeln!(dict_file, "mag_id\tmag_name")?;
    for (id, name) in mag_names.iter().enumerate() {
        writeln!(dict_file, "{}\t{}", id, name)?;
    }
    dict_file.flush()?;

    // Family dictionary: write the mapping of family_id to representative_gene packed id to a TSV file
    let mut family_dict_file = BufWriter::new(
        File::create(Path::new(&out_dir).join("family_dictionary.tsv"))?,
    );
    writeln!(family_dict_file, "family_id\trepresentative_gene")?;
    for (id, rep_id) in family_reps.iter().enumerate() {
        writeln!(family_dict_file, "{}\t{}", id, rep_id)?;
    }
    family_dict_file.flush()?;

    // Print summary statistics to standard error
    let elapsed = start.elapsed().as_secs_f64();
    eprintln!("\nDone.");
    eprintln!("  Rows read   : {}", rows_read);
    eprintln!("  Rows written: {}", rows_written);
    eprintln!("  MAGs seen   : {}", mag_names.len());
    eprintln!("  Families    : {}", family_reps.len());
    eprintln!("  Time        : {:.1}s", elapsed);
    eprintln!("  Throughput  : {:.1}M rows/s", rows_read as f64 / elapsed / 1_000_000.0);

    Ok(())
}

// Unit tests for the utility functions
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unquote_strips_matching_quotes() {
        assert_eq!(unquote("\"Mag1|Mag1_gene1\""), "Mag1|Mag1_gene1");
    }

    #[test]
    fn unquote_leaves_unquoted_text_alone() {
        assert_eq!(unquote("Mag1|Mag1_gene1"), "Mag1|Mag1_gene1");
    }

    #[test]
    fn unquote_leaves_a_lone_quote_char_alone() {
        // len < 2, so starts_with('"') && ends_with('"') can't both hold
        // meaningfully -- must not strip anything down to an empty string.
        assert_eq!(unquote("\""), "\"");
    }

    #[test]
    fn extract_mag_takes_everything_before_first_pipe() {
        assert_eq!(extract_mag("Mag75|Mag75_gene51"), "Mag75");
    }

    #[test]
    fn extract_mag_returns_whole_string_if_no_pipe() {
        assert_eq!(extract_mag("no_pipe_here"), "no_pipe_here");
    }

    #[test]
    fn mag_id_for_assigns_sequential_ids_in_first_seen_order() {
        let mut dict = HashMap::new();
        let mut names = Vec::new();
        assert_eq!(mag_id_for("Mag1", &mut dict, &mut names), 0);
        assert_eq!(mag_id_for("Mag75", &mut dict, &mut names), 1);
        // seeing "Mag1" again must return its existing id, not a new one
        assert_eq!(mag_id_for("Mag1", &mut dict, &mut names), 0);
        assert_eq!(names, vec!["Mag1".to_string(), "Mag75".to_string()]);
    }

    #[test]
    fn encode_gene_id_matches_the_known_real_example() {
        // Mag75|Mag75_gene51, mag_id=1 -> 1*4294967296 + 51 = 4294967347
        // (the exact value this project has cross-checked by hand repeatedly)
        assert_eq!(encode_gene_id("Mag75|Mag75_gene51", "Mag75", 1), 4_294_967_347);
    }

    #[test]
    fn encode_gene_id_round_trips_through_the_packing_scheme() {
        let mag_id = 42u32;
        let local = 12345u32;
        let packed = encode_gene_id("Magx|Magx_gene12345", "Magx", mag_id);
        assert_eq!((packed >> 32) as u32, mag_id);
        assert_eq!((packed & 0xFFFF_FFFF) as u32, local);
    }

    #[test]
    #[should_panic(expected = "gene id missing '|'")]
    fn encode_gene_id_panics_without_a_pipe() {
        encode_gene_id("Mag1_gene1", "Mag1", 0);
    }

    #[test]
    #[should_panic(expected = "gene suffix does not start with MAG")]
    fn encode_gene_id_panics_on_mag_mismatch() {
        encode_gene_id("Mag1|Mag2_gene1", "Mag1", 0);
    }

    #[test]
    #[should_panic(expected = "gene id missing '_gene' suffix")]
    fn encode_gene_id_panics_without_gene_suffix() {
        encode_gene_id("Mag1|Mag1_foo1", "Mag1", 0);
    }

    #[test]
    #[should_panic(expected = "gene local number not an integer")]
    fn encode_gene_id_panics_on_non_numeric_local() {
        encode_gene_id("Mag1|Mag1_geneABC", "Mag1", 0);
    }
}