# Proposed-methodology final evaluation

This is a proposed-only, alignment-free evaluation bundle. It is not a biological benchmark and does not overwrite or merge with the minigraph/PanAligner baseline.

## Datasets and results

| Gene | Assemblies / paths | Nodes | Edges | Overall structural quality | Backend |
|---|---:|---:|---:|---:|---|
| APP | 36 | 72646 | 79053 | 0.983994 | cpu |
| PSEN1 | 124 | 21397 | 23118 | 0.979987 | cpu |
| PSEN2 | 14 | 15635 | 16089 | 0.986480 | cpu |

All assemblies currently listed in `data/metadata/train_manifest.json` were used, without the older eight-per-bucket cap.

## Pipeline

FASTA inputs → deterministic minimizers → path-preserving proposed graph → Step 2 structural/path coherence → Step 3 anchors → Step 4 PG-SCUnK controller.

## Interpretation and limitations

Graph/path consistency is structural coherence derived from observed minimizer paths, not independent biological correctness. Minimizer backend requested: `cpu`; every gene used the CPU, so this evaluation has no GPU acceleration. Baseline alignment metrics and proposed structural metrics are not directly interchangeable.

Runtime: 194.49 seconds. Existing baseline files were read only and are described separately in `baseline_proposed_comparison.md`.
