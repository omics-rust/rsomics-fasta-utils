# rsomics-fasta-utils

A toolkit of small FASTA operations under one binary — counting, filtering,
converting, and reshaping FASTA records. Covers the everyday `seqkit` / `seqtk`
verbs.

## Install

```
cargo install rsomics-fasta-utils
```

## Usage

```
rsomics-fasta-utils count genome.fa
rsomics-fasta-utils revcomp genome.fa -o rc.fa
rsomics-fasta-utils filter -m 500 -M 2000 genome.fa -o filtered.fa
```

Subcommands: `chroms`, `composition`, `count`, `dedup`, `extract`, `filter`,
`gc-content`, `grep`, `head`, `kmers`, `len`, `rename`, `revcomp`, `sample`,
`shuffle`, `sort`, `split`, `tab`, `to-bed`, `to-fastq`, `unique`, `upper`,
`window`, `wrap`. Run `rsomics-fasta-utils <cmd> --help` for each command's
flags; most take a positional FASTA and write to `-o` (`-` = stdout).

## Origin

Independent Rust reimplementation of common FASTA operations offered by `seqkit`
and `seqtk`. `count`, `head`, and `revcomp` are checked against `seqkit stats`,
`seqkit head`, and `seqkit seq -rp` respectively; the remaining subcommands are
self-tested for correctness rather than byte-compared to a specific upstream, so
they are compatible in operation rather than guaranteed byte-identical (in
particular `sample` / `shuffle` use their own RNG and will not reproduce
seqkit's exact selection for a given seed).

License: MIT OR Apache-2.0.
Upstream credit: [seqkit](https://github.com/shenwei356/seqkit) (MIT),
[seqtk](https://github.com/lh3/seqtk) (MIT).
