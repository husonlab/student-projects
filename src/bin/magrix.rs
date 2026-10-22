// Command-line entry point: runs the whole pipeline with one command.
//
//   magrix build <clusters.tsv> <out_dir> [--rows-per-run N] [--keep-intermediates]
//
// It calls the stage programs (pangenome_matrices, build_csc, build_csr,
// build_members, pack_container), which are expected next to this executable,
// in order, derives the number of gene families and MAGs from the dictionaries
// written by stage 1, deletes intermediate files once they are no longer needed
// and prints the run time of every stage. Buffer sizes can be set through the
// MAGRIX_*_BUF_MB environment variables, which the stages inherit.
//
// Each stage is a separate program, so its memory is freed when it exits and
// a failing stage stops the build with its intermediate files left in place.

// Imports
use std::env;
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{self, Command};
use std::time::Instant;

const USAGE: &str = "\
magrix - build and query large MAG x gene family matrices

Usage:
  magrix build <clusters.tsv> <out_dir> [options]

Builds paralog.container, mag_dictionary.tsv and family_dictionary.tsv in
<out_dir> from an MMseqs2 clustering table (representative<TAB>member).

Options:
  --rows-per-run N       entries sorted in memory per run when building the
                         row-wise matrix (default 20000000, about 240 MB)
  --keep-intermediates   keep the intermediate tables and array files
  -h, --help             show this help
  -V, --version          show the version

Environment variables (sizes in MiB):
  MAGRIX_READ_BUF_MB (8)   MAGRIX_WRITE_BUF_MB (8)
  MAGRIX_MERGE_BUF_MB (4)  MAGRIX_COPY_BUF_MB (16)";

// Stage programs in the order they run
const STAGES: [&str; 5] = ["pangenome_matrices", "build_csc", "build_csr", "build_members", "pack_container"];

// Parsed arguments of `magrix build`
struct BuildOptions {
    input: PathBuf,
    out_dir: PathBuf,
    rows_per_run: Option<u64>,
    keep_intermediates: bool,
}

// Print an error message and exit with code 1
fn fail(msg: &str) -> ! {
    eprintln!("magrix: {}", msg);
    process::exit(1);
}

fn main() {
    // Skip the program name; the first argument is the command ("build", --help, --version)
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("build") => match parse_build_args(&args[1..]) {
            Ok(Some(opts)) => build(&opts),
            Ok(None) => println!("{}", USAGE),
            Err(msg) => {
                eprintln!("magrix: {}\n\n{}", msg, USAGE);
                process::exit(2);
            }
        },
        Some("-h") | Some("--help") => println!("{}", USAGE),
        Some("-V") | Some("--version") => println!("magrix {}", env!("CARGO_PKG_VERSION")),
        Some(other) => {
            eprintln!("magrix: unknown command '{}'\n\n{}", other, USAGE);
            process::exit(2);
        }
        None => {
            eprintln!("{}", USAGE);
            process::exit(2);
        }
    }
}

// Ok(None) means help was requested.
fn parse_build_args(args: &[String]) -> Result<Option<BuildOptions>, String> {
    let mut positional: Vec<&String> = Vec::new();
    let mut rows_per_run = None;
    let mut keep_intermediates = false;
    let mut i = 0;
    while i < args.len() {
        let arg = args[i].as_str();
        match arg {
            "-h" | "--help" => return Ok(None),
            "--keep-intermediates" => keep_intermediates = true,
            "--rows-per-run" => {
                i += 1;
                let value = args.get(i).ok_or("--rows-per-run needs a value")?;
                rows_per_run = Some(parse_rows_per_run(value)?);
            }
            _ if arg.starts_with("--rows-per-run=") => {
                rows_per_run = Some(parse_rows_per_run(&arg["--rows-per-run=".len()..])?);
            }
            _ if arg.starts_with('-') => return Err(format!("unknown option '{}'", arg)),
            _ => positional.push(&args[i]),
        }
        i += 1;
    }
    if positional.len() != 2 {
        return Err("build needs exactly two arguments: <clusters.tsv> <out_dir>".to_string());
    }
    Ok(Some(BuildOptions {
        input: PathBuf::from(positional[0]),
        out_dir: PathBuf::from(positional[1]),
        rows_per_run,
        keep_intermediates,
    }))
}

fn parse_rows_per_run(value: &str) -> Result<u64, String> {
    value
        .parse::<u64>()
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| format!("--rows-per-run must be a positive integer, got '{}'", value))
}

// The stage programs are installed next to this executable.
fn stage_path(name: &str) -> PathBuf {
    let exe = env::current_exe().unwrap_or_else(|e| fail(&format!("cannot locate own executable: {}", e)));
    exe.parent().expect("executable has a parent directory").join(name)
}

// Run one stage program, wait for it and return its run time in seconds; exit if it fails.
// The stage's own progress output goes straight to the terminal.
fn run_stage(label: &str, program: &Path, args: &[String]) -> f64 {
    eprintln!("\n== {} ==", label);
    let start = Instant::now();
    let status = Command::new(program)
        .args(args)
        .status()
        .unwrap_or_else(|e| fail(&format!("cannot start {}: {}", program.display(), e)));
    if !status.success() {
        fail(&format!("stage '{}' failed ({}); intermediate files were left in place", label, status));
    }
    start.elapsed().as_secs_f64()
}

// Number of data rows of a TSV file written by stage 1 (all lines minus the header).
fn count_entries(path: &Path) -> u64 {
    let mut file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {}: {}", path.display(), e)));
    let mut buf = vec![0u8; 1 << 20];
    let mut lines: u64 = 0;
    loop {
        let n = file.read(&mut buf).unwrap_or_else(|e| fail(&format!("cannot read {}: {}", path.display(), e)));
        if n == 0 {
            break;
        }
        lines += buf[..n].iter().filter(|&&b| b == b'\n').count() as u64;
    }
    lines.saturating_sub(1)
}

// Delete intermediate files; a file that can't be deleted only gives a warning
fn remove(paths: &[PathBuf]) {
    for p in paths {
        if let Err(e) = fs::remove_file(p) {
            eprintln!("magrix: warning: could not delete {}: {}", p.display(), e);
        }
    }
}

// Format a byte count as e.g. "23.2 GiB" for the summary
fn human_size(bytes: u64) -> String {
    let units = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit < units.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    format!("{:.1} {}", value, units[unit])
}

// Run the whole pipeline. Each intermediate file is deleted right after the last stage that reads it,
// so the peak disk usage stays well below the sum of all intermediates.
fn build(opts: &BuildOptions) {
    // Check everything that can be checked up front, so a typo doesn't fail after a long stage 1
    if !opts.input.is_file() {
        fail(&format!("input file not found: {}", opts.input.display()));
    }
    let programs: Vec<PathBuf> = STAGES.iter().map(|s| stage_path(s)).collect();
    for p in &programs {
        if !p.is_file() {
            fail(&format!("stage program not found: {} (build all programs with `cargo build --release`)", p.display()));
        }
    }
    fs::create_dir_all(&opts.out_dir)
        .unwrap_or_else(|e| fail(&format!("cannot create {}: {}", opts.out_dir.display(), e)));

    let out = &opts.out_dir;
    // Paths are passed to the stages as strings
    let s = |p: &Path| p.to_string_lossy().into_owned();
    // All array files are named <out_dir>/paralog.<array>
    let prefix = s(&out.join("paralog"));
    let paralog_tsv = out.join("mmseqs_long_paralog.tsv");
    let gene_list_tsv = out.join("mmseqs_long_gene_list.tsv");
    let mut times: Vec<(&str, f64)> = Vec::new();

    times.push(("1 preprocessing", run_stage("1 preprocessing", &programs[0], &[s(&opts.input), s(out)])));

    // build_csc and build_csr need the matrix dimensions: one dictionary line per family / MAG
    let num_families = count_entries(&out.join("family_dictionary.tsv"));
    let num_mags = count_entries(&out.join("mag_dictionary.tsv"));
    eprintln!("\n{} MAGs, {} gene families", num_mags, num_families);

    // Both 2a and 2b read paralog.tsv, so it is only deleted after 2b
    times.push(("2a CSC", run_stage("2a CSC", &programs[1], &[s(&paralog_tsv), num_families.to_string(), prefix.clone()])));

    let mut csr_args = vec![s(&paralog_tsv), num_mags.to_string(), prefix.clone()];
    if let Some(n) = opts.rows_per_run {
        csr_args.push(n.to_string());
    }
    times.push(("2b CSR", run_stage("2b CSR", &programs[2], &csr_args)));
    if !opts.keep_intermediates {
        remove(&[paralog_tsv]);
    }

    times.push(("3 members", run_stage("3 members", &programs[3], &[s(&gene_list_tsv), prefix.clone()])));
    if !opts.keep_intermediates {
        remove(&[gene_list_tsv]);
    }

    // Once packed, the separate array files are no longer needed
    times.push(("4 container", run_stage("4 container", &programs[4], &[prefix.clone()])));
    if !opts.keep_intermediates {
        let arrays = [
            "csc.indptr", "csc.indices", "csc.data",
            "csr.indptr", "csr.indices", "csr.data",
            "members.family_ptr", "members.gene_numbers",
        ];
        let paths: Vec<PathBuf> = arrays.iter().map(|a| PathBuf::from(format!("{}.{}", prefix, a))).collect();
        remove(&paths);
    }

    let total: f64 = times.iter().map(|(_, t)| t).sum();
    eprintln!("\n== done ==");
    for (label, t) in &times {
        eprintln!("  {:<16} {:>9.1} s", label, t);
    }
    eprintln!("  {:<16} {:>9.1} s", "total", total);
    eprintln!("\nOutput in {}:", out.display());
    for name in ["paralog.container", "mag_dictionary.tsv", "family_dictionary.tsv"] {
        let size = fs::metadata(out.join(name)).map(|m| m.len()).unwrap_or(0);
        eprintln!("  {:<24} {:>12}", name, human_size(size));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_positional_arguments_and_options() {
        let opts = parse_build_args(&args(&["in.tsv", "out", "--rows-per-run", "500", "--keep-intermediates"]))
            .unwrap()
            .unwrap();
        assert_eq!(opts.input, PathBuf::from("in.tsv"));
        assert_eq!(opts.out_dir, PathBuf::from("out"));
        assert_eq!(opts.rows_per_run, Some(500));
        assert!(opts.keep_intermediates);
    }

    #[test]
    fn accepts_option_with_equals_sign() {
        let opts = parse_build_args(&args(&["in.tsv", "out", "--rows-per-run=7"])).unwrap().unwrap();
        assert_eq!(opts.rows_per_run, Some(7));
    }

    #[test]
    fn help_returns_none() {
        assert!(parse_build_args(&args(&["--help"])).unwrap().is_none());
    }

    #[test]
    fn rejects_wrong_argument_count_and_bad_values() {
        assert!(parse_build_args(&args(&["only_one"])).is_err());
        assert!(parse_build_args(&args(&["in", "out", "--rows-per-run", "0"])).is_err());
        assert!(parse_build_args(&args(&["in", "out", "--unknown"])).is_err());
    }
}
