use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

/// Lookup table: index = ASCII byte → complement byte.
/// Handles upper and lower case IUPAC bases; everything else maps to itself.
static COMP: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    t[b'A' as usize] = b'T';
    t[b'a' as usize] = b't';
    t[b'T' as usize] = b'A';
    t[b't' as usize] = b'a';
    t[b'C' as usize] = b'G';
    t[b'c' as usize] = b'g';
    t[b'G' as usize] = b'C';
    t[b'g' as usize] = b'c';
    t[b'N' as usize] = b'N';
    t[b'n' as usize] = b'n';
    t[b'R' as usize] = b'Y';
    t[b'r' as usize] = b'y';
    t[b'Y' as usize] = b'R';
    t[b'y' as usize] = b'r';
    t[b'S' as usize] = b'S';
    t[b's' as usize] = b's';
    t[b'W' as usize] = b'W';
    t[b'w' as usize] = b'w';
    t[b'K' as usize] = b'M';
    t[b'k' as usize] = b'm';
    t[b'M' as usize] = b'K';
    t[b'm' as usize] = b'k';
    t[b'B' as usize] = b'V';
    t[b'b' as usize] = b'v';
    t[b'V' as usize] = b'B';
    t[b'v' as usize] = b'b';
    t[b'D' as usize] = b'H';
    t[b'd' as usize] = b'h';
    t[b'H' as usize] = b'D';
    t[b'h' as usize] = b'd';
    t
};

pub fn revcomp(input: &Path, output: &mut dyn Write) -> Result<u64> {
    if std::fs::metadata(input).is_ok_and(|m| m.len() == 0) {
        return Ok(0);
    }
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    // Single reusable scratch buffer: no per-record allocation in the hot path.
    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let seq = rec.seq();

        buf.clear();
        buf.reserve(seq.len());
        // Reverse-complement in one pass using the lookup table.
        for &b in seq.iter().rev() {
            buf.push(COMP[b as usize]);
        }

        out.write_all(b">").map_err(RsomicsError::Io)?;
        out.write_all(rec.id()).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        out.write_all(&buf).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
        count += 1;
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}
