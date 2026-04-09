# Genomic Feature Analyzer (GFA) - R Refactor

This project contains the **gfaR** R package, a refactor of the original Genomic Feature Analyzer (GFA) C-program. GFA is designed to identify non-B DNA motifs in genomic sequences, which are associated with various biological processes and genomic instabilities. The original GFA program was authored by 
Cer, R. Z., Donohue, D. E., Mudunuri, U. S., Temiz, N. A., Loss, M. A., Starner, N. J., ... & Stephens, R. M. (2013). **Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools.** *Nucleic Acids Research*, 41(D1), D94-D100. doi: 10.1093/nar/gks955.

The gfaR pacakge and this repo was generated through the use of Gemini CLI

## Project Structure

- `gfaR/`: The complete R package source code, including C++ implementations of the core algorithms via Rcpp.
- `gfa_c_program/`: The original C source code used as the foundation for this refactor.
- `test_data/`: Contains gold-standard FASTA files and expected GFF/TSV outputs for validation.
- `references/`: Documentation and original publications related to the GFA program.
- `compare_gfa_versions.R`: A specialized script to run both versions side-by-side and verify functional parity.

## The gfaR R Package

The `gfaR` package provides an idiomatic R interface to the optimized GFA algorithms. It supports the detection of seven major non-B DNA motif types:

1.  **A-Phased Repeats (APR)**: Associated with bent DNA.
2.  **Direct Repeats (DR)**: Associated with slipped DNA.
3.  **G-quadruplexes (GQ)**: Four-stranded DNA structures.
4.  **Inverted Repeats (IR)**: Associated with cruciforms.
5.  **Mirror Repeats (MR)**: Associated with triplexes.
6.  **Short Tandem Repeats (STR)**: Microsatellites.
7.  **Z-DNA Motifs**: Left-handed helical DNA.

### Installation

To install the package from this repository, use `devtools`:

```r
# Install devtools if you haven't already
# install.packages("devtools")

# Install gfaR
devtools::install("gfaR")
```

### Quick Start

```r
library(gfaR)

# Analyze a sequence for all motif types
sequence <- "AGTGCAACCCAGAGGGCAGGATTTCCTGCTGGACTTTGAAATCCAACCCGGTCACCTACCCGCGCGACTG"
results <- gfa_analyze(sequence)

# Results are returned as a named list of data.frames
print(results$gq)  # View G-quadruplex results
```

## Validation & Parity

A core requirement of this project was maintaining 100% functional parity with the original C implementation. You can verify this by running the included comparison script:

```bash
Rscript compare_gfa_versions.R
```

This script automatically compiles the original C program, runs both versions on the same test sequence, and performs a side-by-side coordinate comparison.

## Citation

If you use this software in your research, please cite the original GFA publication:

> Cer, R. Z., Donohue, D. E., Mudunuri, U. S., Temiz, N. A., Loss, M. A., Starner, N. J., ... & Stephens, R. M. (2013). **Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools.** *Nucleic Acids Research*, 41(D1), D94-D100. doi: 10.1093/nar/gks955.

You can also view the formal citation in R using:
```r
citation("gfaR")
```
