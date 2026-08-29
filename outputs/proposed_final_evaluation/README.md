# Proposed-methodology final evaluation

This is a proposed-only, alignment-free evaluation bundle. It is not a biological benchmark and does not overwrite or merge with the minigraph/PanAligner baseline.

## Datasets and results

| Gene | Assemblies / paths | Nodes | Edges | Overall structural quality | Backend |
|---|---:|---:|---:|---:|---|
| APP | 36 | 72646 | 79053 | 0.886395 | cpu |
| PSEN1 | 124 | 21397 | 23118 | 0.848715 | cpu |
| PSEN2 | 14 | 15635 | 16089 | 0.988508 | cpu |

All assemblies currently listed in `data/metadata/train_manifest.json` were used, without the older eight-per-bucket cap.

## Pipeline

FASTA inputs → deterministic CPU minimizers → path-preserving proposed graph → Step 2 structural/path coherence → Step 3 anchors → Step 4 PG-SCUnK controller.

## Interpretation and limitations

Graph/path consistency is structural coherence derived from observed minimizer paths, not independent biological correctness. The proposed method has no GPU acceleration in this evaluation; the recorded backend is CPU. Baseline alignment metrics and proposed structural metrics are not directly interchangeable.

Runtime: 995.17 seconds. Existing baseline files were read only and are described separately in `baseline_proposed_comparison.md`.
