//! Alphabet-aware complement, matching how `seqkit seq -rp` picks one complement
//! table for a whole file from the first record's alphabet.

/// Which complement to apply, decided once from the first record. RNA pairs A
/// with U; a non-nucleotide alphabet is reversed without complementing.
#[derive(Clone, Copy)]
pub(crate) enum ComplementKind {
    Dna,
    Rna,
    None,
}

// seqkit's DNAredundant / RNAredundant letter sets: bases, IUPAC codes, gaps.
const DNA_ALPHABET: &[u8] = b"acgtryswkmbdhvACGTRYSWKMBDHVnN -.";
const RNA_ALPHABET: &[u8] = b"acguryswkmbdhvACGURYSWKMBDHVnN -.";

impl ComplementKind {
    pub(crate) fn guess(first_seq: &[u8]) -> Self {
        // seqkit guesses from the first 10 kb of the first record.
        let head = &first_seq[..first_seq.len().min(10_000)];
        if head.is_empty() {
            ComplementKind::None
        } else if head.iter().all(|b| DNA_ALPHABET.contains(b)) {
            ComplementKind::Dna
        } else if head.iter().all(|b| RNA_ALPHABET.contains(b)) {
            ComplementKind::Rna
        } else {
            ComplementKind::None
        }
    }
}

pub(crate) fn complement(b: u8, kind: ComplementKind) -> u8 {
    match kind {
        ComplementKind::None => b,
        ComplementKind::Dna => match b {
            b'A' => b'T',
            b'a' => b't',
            b'T' => b'A',
            b't' => b'a',
            other => complement_shared(other),
        },
        ComplementKind::Rna => match b {
            b'A' => b'U',
            b'a' => b'u',
            b'U' => b'A',
            b'u' => b'a',
            other => complement_shared(other),
        },
    }
}

/// Bases whose complement is identical in DNA and RNA — C/G/N and the IUPAC
/// ambiguity codes. An unmapped byte passes through unchanged, as seqkit leaves
/// a letter outside the alphabet as-is.
fn complement_shared(b: u8) -> u8 {
    match b {
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'N' => b'N',
        b'n' => b'n',
        b'R' => b'Y',
        b'r' => b'y',
        b'Y' => b'R',
        b'y' => b'r',
        b'S' => b'S',
        b's' => b's',
        b'W' => b'W',
        b'w' => b'w',
        b'K' => b'M',
        b'k' => b'm',
        b'M' => b'K',
        b'm' => b'k',
        b'B' => b'V',
        b'b' => b'v',
        b'V' => b'B',
        b'v' => b'b',
        b'D' => b'H',
        b'd' => b'h',
        b'H' => b'D',
        b'h' => b'd',
        other => other,
    }
}
