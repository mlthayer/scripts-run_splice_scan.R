This repository contains a reproducible R script to:
1) sanitize WT and MUT cDNA sequences to A/C/G/T,
2) find the first nucleotide difference,
3) scan a ±N nt window around the difference for donor (GT) and acceptor (AG) motifs,
4) compute simple heuristic splice-site scores, and
5) output a WT vs MUT comparison table (CSV).

## Requirements
- R (>= 4.0 recommended)
- No external packages required (base R only)

## Run
From the repository root:

```bash
Rscript scripts/run_splice_scan.R --wt "ACGT..." --mut "ACGT..." --radius 200 --out outputs/results.csv

## Code availability

The analysis pipeline is archived on Zenodo:
mlthayer. (2026). mlthayer/scripts-run_splice_scan.R: WT vs MUT splice motif scan (v0.1.0) (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.18295243
