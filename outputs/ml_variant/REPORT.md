# Model B feasibility study (research label recovery, not diagnosis)

**Label source, version and criteria are unknown; results describe recovery of these labels, not disease association.**

Data: 220 sequences, each the GRCh38 locus reference with one inserted variant ({'SNV': 205, 'multi_SNV': 3, 'insertion': 4, 'deletion': 8}); 110 labelled unhealthy.

Evaluation: StratifiedGroupKFold(5) over 1000 bp genomic windows, seeds [0, 1, 2, 3, 4]; mean ± SD over seeds. `uncertain` margin chosen by inner grouped CV for target accuracy 0.8.

| Model | Balanced accuracy | MCC | ROC-AUC | Confident predictions (mean) | Accuracy when confident |
|---|---|---|---|---|---|
| majority | 0.481 ± 0.002 | -0.039 ± 0.004 | – | – | – |
| nearest_neighbour | 0.628 ± 0.009 | 0.256 ± 0.019 | – | – | – |
| position_only | 0.500 ± 0.026 | -0.000 ± 0.053 | 0.490 ± 0.044 | 0.0 of 220 | – |
| context_only | 0.550 ± 0.015 | 0.103 ± 0.032 | 0.556 ± 0.021 | 0.2 of 220 | 1.000 ± 0.000 |
| context_and_position | 0.555 ± 0.009 | 0.109 ± 0.017 | 0.582 ± 0.015 | 0.4 of 220 | 1.000 ± 0.000 |

Original seed-42 per-sequence split, nearest training variant: 33/46 correct (balanced accuracy 0.717); the random split leaks through position.

Read `context_only` against `position_only` and `nearest_neighbour`: anything the position models achieve reflects where the labelled variants were chosen, not a property learned from sequence context.
