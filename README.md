# boquila

Tool for generating synthetic next-generation sequencing reads with same nucleotide distribution as model reads.

## Usage

```
boquila 0.3.1
Generate NGS reads with same nucleotide distribution as input file
Generated reads will be written to stdout
By default input and output format is FASTQ

USAGE:
    boquila [FLAGS] [OPTIONS] <src> --ref <FILE> --regions <FILE>

ARGS:
    <src>    Model file

FLAGS:
        --fasta      Change input and output format to FASTA
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --bed <FILE>        File name in which the simulated reads will be saved in BED format
        --ref <FILE>        Reference FASTA
        --regions <FILE>    RON formatted file containing genomic regions that generated reads will
                            be selected from
        --seed <INT>        Random number seed. If not provided system's default source of entropy
                            will be used instead.
```

Generated reads will be written to stdout

Sample `regions` file for Homo sapiens (human) genome assembly GRCh38 (hg38) is provided as `GRCh38.ron`

## Installation

boquila is written in Rust, so you'll need to grab a [Rust installation](https://www.rust-lang.org/) in order to install or compile it.

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
0.3.1
```