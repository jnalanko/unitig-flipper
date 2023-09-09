# Unitig Flipper

This program takes in a set of unitigs as a FASTA file, and re-orients the unitigs heuristically in an attempt to minimize the number of dummy nodes in the SBWT of the k-mers.

## Compiling

First, install the [Rust](https://www.rust-lang.org/tools/install). Then:

```
git submodule update --init
cargo build --release
```

This produces the binary to `target/release/unitig_flipper`. If you want to install the program to `$PATH`, run `cargo install --path .`

## Usage

```
Usage: unitig_flipper --input <input> --output <output> -k <k>

Options:
  -i, --input <input>    Input FASTA or FASTQ file, possibly gzipped
  -o, --output <output>  Output FASTA or FASTQ file, possibly gzipped
  -k <k>                 k-mer length
  -h, --help             Print help
```
