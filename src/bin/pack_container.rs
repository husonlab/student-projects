// Pack the 8 separate binary arrays (csc.indptr/indices/data,
// csr.indptr/indices/data, members.family_ptr/gene_numbers) for one prefix
// into a single container file. The two text dictionaries stay separate --
// so 10 files total becomes 3: one container + 2 TSVs.
//
// Container layout:
//   [0, HEADER_SIZE)   fixed-size header block (one page, 4096 bytes)
//     bytes [0,4)   magic "PANG"
//     bytes [4,8)   u32 num_arrays
//     then num_arrays fixed 64-byte records:
//       [0,32)  array name, UTF-8, zero-padded (e.g. "csc.indptr")
//       [32,33) dtype code: 0 = u32, 1 = u64, 2 = u16
//       [33,41) offset (u64 LE) -- byte offset into this file
//       [41,49) count (u64 LE)  -- number of elements, not bytes
//       [49,64) reserved / zero padding
//   [HEADER_SIZE, ...) array data, back to back. Each array is individually
//                      padded at the end to the next 4096-byte boundary, so
//                      every array starts on its own page: the mmap system
//                      call only accepts page-aligned offsets, so any reader
//                      can map an array directly at its offset (numpy.memmap
//                      would also cope with unaligned offsets by rounding
//                      down internally, but other readers would not).
//
// The arrays are copied in chunks of MAGRIX_COPY_BUF_MB, so memory stays at a
// few tens of MB regardless of the container size.

// Imports
use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};

const HEADER_SIZE: u64 = 4096;
const RECORD_SIZE: usize = 64;
const PAGE: u64 = 4096;
// Copy buffer in MiB, overridable via MAGRIX_COPY_BUF_MB.
const DEFAULT_COPY_BUF_MB: usize = 16;

// Buffer size in bytes: the environment variable `var` (in MiB) if set, else the default.
fn buffer_size(var: &str, default_mb: usize) -> usize {
    let mb = match env::var(var) {
        Ok(v) => v.parse::<usize>().ok().filter(|&n| n > 0)
            .unwrap_or_else(|| panic!("{} must be a positive integer (MiB), got {:?}", var, v)),
        Err(_) => default_mb,
    };
    mb * 1024 * 1024
}

// One array to pack: its file name suffix (= its name in the header) and its element type
struct ArraySpec {
    suffix: &'static str,
    dtype_code: u8, // 0 = u32, 1 = u64, 2 = u16
}

// All arrays in the order they are written into the container
// (the reader looks them up by name, so the order only affects the layout)
const ARRAYS: [ArraySpec; 8] = [
    ArraySpec { suffix: "csc.indptr", dtype_code: 1 },
    ArraySpec { suffix: "csc.indices", dtype_code: 0 },
    ArraySpec { suffix: "csc.data", dtype_code: 2 },
    ArraySpec { suffix: "csr.indptr", dtype_code: 1 },
    ArraySpec { suffix: "csr.indices", dtype_code: 0 },
    ArraySpec { suffix: "csr.data", dtype_code: 2 },
    ArraySpec { suffix: "members.family_ptr", dtype_code: 1 },
    ArraySpec { suffix: "members.gene_numbers", dtype_code: 0 },
];

// Size of one element in bytes for a dtype code
fn elem_size(dtype_code: u8) -> u64 {
    match dtype_code {
        0 => 4,
        1 => 8,
        2 => 2,
        _ => panic!("unknown dtype code {}", dtype_code),
    }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: pack_container <prefix>   (e.g. output/paralog)");
        std::process::exit(1);
    }
    let prefix = &args[1];
    let out_path = format!("{}.container", prefix);

    // Pass 1: read each source file's size (no data read yet) and lay out
    // offsets ahead of time, so the header can be written before the data.
    let mut offset = HEADER_SIZE;
    let mut records: Vec<(String, u8, u64, u64)> = Vec::new(); // (suffix, dtype_code, offset, count)
    for spec in ARRAYS.iter() {
        let path = format!("{}.{}", prefix, spec.suffix);
        let len = std::fs::metadata(&path)
            .unwrap_or_else(|e| panic!("cannot stat {}: {}", path, e))
            .len();
        let esz = elem_size(spec.dtype_code);
        assert!(len % esz == 0, "{} size {} not a multiple of element size {}", path, len, esz);
        let count = len / esz;
        records.push((spec.suffix.to_string(), spec.dtype_code, offset, count));
        // Round the length up to the next multiple of 4096, so the next array starts on a page boundary
        let padded_len = ((len + PAGE - 1) / PAGE) * PAGE;
        offset += padded_len;
    }

    let copy_buf = buffer_size("MAGRIX_COPY_BUF_MB", DEFAULT_COPY_BUF_MB);
    eprintln!("Buffers: copy={} MiB", copy_buf >> 20);

    let mut out = BufWriter::with_capacity(copy_buf, File::create(&out_path)?);

    // Write header block: starts as 4096 zero bytes, so unused space and reserved bytes stay zero.
    let mut header = vec![0u8; HEADER_SIZE as usize];
    header[0..4].copy_from_slice(b"PANG");
    header[4..8].copy_from_slice(&(records.len() as u32).to_le_bytes());
    for (i, (suffix, dtype_code, off, count)) in records.iter().enumerate() {
        let base = 8 + i * RECORD_SIZE;
        let name_bytes = suffix.as_bytes();
        assert!(name_bytes.len() <= 32, "array name too long: {}", suffix);
        header[base..base + name_bytes.len()].copy_from_slice(name_bytes);
        header[base + 32] = *dtype_code;
        header[base + 33..base + 41].copy_from_slice(&off.to_le_bytes());
        header[base + 41..base + 49].copy_from_slice(&count.to_le_bytes());
    }
    out.write_all(&header)?;

    // Pass 2: stream-copy each array's bytes, then pad to the next page boundary.
    let mut written: u64 = HEADER_SIZE;
    for (spec, (_, _, target_offset, _)) in ARRAYS.iter().zip(records.iter()) {
        // Safety check: the array must land exactly at the offset already written into the header
        assert_eq!(written, *target_offset, "internal offset mismatch for {}", spec.suffix);
        let path = format!("{}.{}", prefix, spec.suffix);
        // Read straight into one large buffer; a BufReader would add a second, unused one.
        let mut infile = File::open(&path)?;
        let mut buf = vec![0u8; copy_buf];
        let mut this_len: u64 = 0;
        loop {
            let n = infile.read(&mut buf)?;
            if n == 0 {
                break;
            }
            out.write_all(&buf[..n])?;
            this_len += n as u64;
        }
        written += this_len;

        let padded_len = ((this_len + PAGE - 1) / PAGE) * PAGE;
        let pad = padded_len - this_len;
        if pad > 0 {
            out.write_all(&vec![0u8; pad as usize])?;
            written += pad;
        }
        eprintln!("packed {} ({} bytes)", path, this_len);
    }

    out.flush()?;
    eprintln!("Container written: {} ({} bytes total)", out_path, written);
    Ok(())
}
