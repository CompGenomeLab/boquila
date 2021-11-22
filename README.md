# boquila

Tool for generating synthetic next-generation sequencing reads with same nucleotide distribution as model reads.

## Usage

```
boquila 0.5.1
Generate NGS reads with same nucleotide distribution as input file
Generated reads will be written to stdout
By default input and output format is FASTQ

USAGE:
    boquila [FLAGS] [OPTIONS] <src>

ARGS:
    <src>    Model file

FLAGS:
        --fasta         Change input and output format to FASTA
        --inseqFasta    Change the input sequencing format to FASTA
    -h, --help          Prints help information
    -V, --version       Prints version information

OPTIONS:
        --inseq <FILE>      Input sequencing reads to be used instead of reference genome
        --kmer <INT>        Kmer size to be used while calculating frequency [default: 1]
        --bed <FILE>        File name in which the simulated reads will be saved in BED format
        --ref <FILE>        Reference FASTA
        --regions <FILE>    RON formatted file containing genomic regions that generated reads will
                            be selected from
        --seed <INT>        Random number seed. If not provided system's default source of entropy
                            will be used instead.
```

Generated reads will be written to stdout in FASTA or FASTQ format.

If `--bed` option is provided, generated reads also will be written to given file in `BED6` format.

Sample `regions` file for Homo sapiens (human) genome assembly GRCh38 (hg38) is provided as `GRCh38.ron`

If Input Sequencing reads will be used for simulation, they should be provided with `--inseq` argument, instead of using `--ref` and `--regions`.

### Examples

> More detailed example can be found in the [examples](./examples) directory

Simple usage
```
boquila input_reads.fq --ref ref_genome.fa --regions GRCh38.ron > out.fq
```

Using seed for RNG
```
boquila input_reads.fq --ref ref_genome.fa --regions GRCh38.ron --seed 7 > out.fq
```

Using reads that are in FASTA format
```
boquila input_reads.fa --fasta --ref ref_genome.fa --regions GRCh38.ron > out.fa
```

Saving output in BED format
```
boquila input_reads.fq --ref ref_genome.fa --regions GRCh38.ron --bed out.bed > out.fq
```

Using Input Sequencing instead of reference genome
```
boquila input_reads.fq --inseq inputseq_reads.fq > out.fq
```

Using Input Sequencing reads which are in FASTA format
```
boquila input_reads.fq --inseqFasta --inseq inputseq_reads.fa > out.fq
```

## Installation

boquila is written in Rust, so you'll need to grab a [Rust installation](https://www.rust-lang.org/) in order to install or compile it.

The current minimum Rust version is `1.55.0`

Use `cargo` to install

```
$ cargo install --branch main --git https://github.com/CompGenomeLab/boquila.git
```

Or build from source

```
$ git clone https://github.com/CompGenomeLab/boquila.git
$ cd boquila
$ cargo build --release
$ ./target/release/boquila --version
0.5.1
```