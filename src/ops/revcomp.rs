use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

use crate::alphabet::{ComplementKind, complement};

pub fn revcomp(input: &Path, output: &mut dyn Write) -> Result<u64> {
    if std::fs::metadata(input).is_ok_and(|m| m.len() == 0) {
        return Ok(0);
    }
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let mut out = BufWriter::with_capacity(256 * 1024, output);
    // Single reusable scratch buffer: no per-record allocation in the hot path.
    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    // seqkit fixes the complement table from the first record's alphabet, so an
    // RNA file pairs A with U and a protein file is reversed without complementing.
    let mut kind: Option<ComplementKind> = None;
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let rec = record.map_err(|e| RsomicsError::InvalidInput(format!("parsing: {e}")))?;
        let seq = rec.seq();
        let k = *kind.get_or_insert_with(|| ComplementKind::guess(&seq));

        buf.clear();
        buf.reserve(seq.len());
        for &b in seq.iter().rev() {
            buf.push(complement(b, k));
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
