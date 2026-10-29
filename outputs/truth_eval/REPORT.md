# Phase 1 — Truth-labelled correctness at locus scale

Generated 2026-10-29 on the laptop (WSL2 Ubuntu 24.04). Code: `scripts/simulate_truth_reads.py`, `scripts/truth_eval.py` (uncommitted when run; see each `results/*/run_manifest.json` for git state, tool versions and input checksums).

## What was done

Method follows the paper's own scripts (`PanAligner/data/simulated_reads.sh`, `evaluate.sh`):
PBSIM2 (commit eeb5a19, P6C4 HMM model), fixed 10 kb reads at 95% accuracy, true intervals embedded with `paftools.js pbsim2fq`, mapping with `-cx lr`, scoring with `paftools.js mapeval` (paftools 2.30-r1287). Aligners: PanAligner v1.0 and minigraph 0.21-r606. Graphs: `graphs/train/<gene>.combined.gfa` and the linear reference FASTA.

| | 1A `reference` | 1B `heldout` |
|---|---|---|
| Reads simulated from | gene reference backbone (as the paper used CHM13) | the 46 held-out test haplotypes (not in the train graphs) |
| Coverage | 30x of the reference | 10x of each haplotype |
| Seeds | 1, 2, 3 | 1, 2, 3 |
| Reads per seed (APP / PSEN1 / PSEN2) | 872 / 262 / 173 | 2910 / 2816 / 232 |

Coverage differs from the paper (0.5x over a whole genome), because 0.5x of a 57–290 kb locus is 2–14 reads.

## Results (3 seeds pooled; % of simulated reads)

| Reads | Gene | Graph | Aligned (both aligners) | Strict mapeval wrong, PanAligner / minigraph | Wrong after projection |
|---|---|---|---|---|---|
| 1A | APP | combined | 100% | 1.68% / 1.72% | 0 / 0 |
| 1A | PSEN1 | combined | 100% | 15.39% / 15.52% | 0 / 0 |
| 1A | PSEN2 | combined | 100% | 2.70% / 2.89% | 0 / 0 |
| 1A | all | linear | 100% | 0 / 0 | 0 / 0 |
| 1B | APP | combined | 99.99% (1 read unaligned) | 1.65% / 1.70% | 0 / 0 |
| 1B | PSEN1 | combined | 100% / 99.99% | 16.19% / 16.68% | 0 / 1 read |
| 1B | PSEN2 | combined | 100% | 3.02% / 3.16% | 0 / 0 |
| 1B | all | linear | 99.99–100% | 0 / 0 | 0 / 0 |

Full numbers: `results/<source>/by_condition.csv` (mean and SD over seeds), `per_run.csv`, `mapq_curve.csv`, and `<gene>.mapq_curve.png` (Fig. 6-style, strict mapeval).

## How to read this

- **Every strict mapeval error is a scoring artefact.** All of them are alignments whose path passes through a non-reference graph node. mapeval compares the target name with the reference name, so any such path counts as wrong. Projecting those paths back onto the reference (`projected_*`) leaves **no misplaced read** in 1A, and one in 1B (minigraph, PSEN1: a 415 bp partial alignment at MAPQ 28 inside the true locus, below the 10% overlap rule). The paper's numbers use the same strict rule.
- Unaligned reads: one APP read (from `APP_U_16`) that neither aligner places on either graph, and one PSEN1 read that minigraph does not place.
- PanAligner and minigraph are effectively tied here. **This benchmark is too easy to separate them.** Each held-out haplotype is the reference with one SNV, two SNVs (`APP_U_16`) or one small indel (three PSEN1 sequences), and the regions are 57–290 kb. Differences of this kind are expected to appear on whole-genome, structurally diverse graphs (the paper's 10H–95H), which this project does not have.
- `wall_seconds` in `per_run.csv` is one laptop/WSL call and is not a runtime benchmark.

## Scope

This is a locus-scale correctness check using the paper's method. It is **not** a reproduction of the paper's whole-genome benchmark tables or Fig. 6, and GraphAligner is not included (Phase 2, HPC).

## Files

`reads/` (211 MB) and the per-seed GAFs under `results/<source>/<gene>/` (273 MB) are regenerable from the scripts and seeds; keep them out of Git. The CSV, PNG, JSON and this report (~0.4 MB) are the results to keep.
