# Model B feasibility study (research label recovery, not diagnosis)

**Label source, version and criteria are unknown; results describe recovery of these labels, not disease association.** No scientifically valid classifier is claimed, and nothing here is a clinical finding.

Data: 220 sequences, each the GRCh38 locus reference with one inserted variant.

Evaluation: StratifiedGroupKFold(5) over 1000 bp genomic windows, seeds [0, 1, 2, 3, 4]; mean ± SD over seeds. The 95% interval resamples whole windows using the first seed's predictions. `uncertain` margin chosen by inner grouped CV for target accuracy 0.8.

## Class balance

| Gene | Healthy label | Unhealthy label |
|---|---|---|
| APP | 23 | 23 |
| PSEN1 | 78 | 78 |
| PSEN2 | 9 | 9 |

## Overall results

| Model | Balanced accuracy | 95% interval | MCC | ROC-AUC | Confident predictions (mean) | Accuracy when confident |
|---|---|---|---|---|---|---|
| majority | 0.481 ± 0.002 | 0.359–0.604 | -0.039 ± 0.004 | – | – | – |
| nearest_neighbour | 0.628 ± 0.009 | 0.505–0.752 | 0.256 ± 0.019 | – | – | – |
| position_knn | 0.646 ± 0.014 | 0.551–0.761 | 0.338 ± 0.035 | – | – | – |
| position_only | 0.500 ± 0.026 | 0.406–0.645 | -0.000 ± 0.053 | 0.490 ± 0.044 | 0.0 of 220 | – |
| context_only | 0.550 ± 0.015 | 0.502–0.649 | 0.103 ± 0.032 | 0.556 ± 0.021 | 0.2 of 220 | 1.000 ± 0.000 |
| context_and_position | 0.555 ± 0.009 | 0.483–0.648 | 0.109 ± 0.017 | 0.582 ± 0.015 | 0.4 of 220 | 1.000 ± 0.000 |

## Per-gene results

Balanced accuracy, with the MCC after the slash and the 95% interval in brackets.

| Model | APP | PSEN1 | PSEN2 |
|---|---|---|---|
| majority | 0.500 / 0.003 [0.026–0.714] | 0.477 / -0.046 [0.350–0.601] | 0.467 / -0.075 [0.247–0.875] |
| nearest_neighbour | 0.783 / 0.628 [0.610–0.969] | 0.577 / 0.156 [0.454–0.696] | 0.678 / 0.385 [0.556–1.000] |
| position_knn | 0.696 / 0.493 [0.556–0.906] | 0.632 / 0.306 [0.509–0.771] | 0.644 / 0.313 [0.490–1.000] |
| position_only | 0.343 / -0.425 [0.250–0.500] | 0.579 / 0.163 [0.468–0.733] | 0.211 / -0.587 [0.000–0.350] |
| context_only | 0.470 / -0.069 [0.404–0.905] | 0.587 / 0.178 [0.538–0.691] | 0.433 / -0.137 [0.125–0.677] |
| context_and_position | 0.552 / 0.112 [0.440–0.928] | 0.577 / 0.154 [0.507–0.688] | 0.367 / -0.278 [0.077–0.583] |

## Window-size sensitivity

Balanced accuracy when whole windows of this size are kept out of training.

| Window | majority | nearest_neighbour | position_knn | position_only | context_only |
|---|---|---|---|---|---|
| 1000 bp | 0.481 ± 0.002 | 0.628 ± 0.009 | 0.646 ± 0.014 | 0.500 ± 0.026 | 0.550 ± 0.015 |
| 5000 bp | 0.471 ± 0.012 | 0.570 ± 0.045 | 0.602 ± 0.008 | 0.505 ± 0.034 | 0.528 ± 0.024 |
| 20000 bp | 0.383 ± 0.005 | 0.450 ± 0.000 | 0.591 ± 0.000 | 0.398 ± 0.002 | 0.453 ± 0.010 |

## Position separation

AUC of position against the label in each gene (0.5 = none, 1.0 = complete, either direction).

| Gene | Position AUC | Healthy-label positions (min / median / max) | Unhealthy-label positions (min / median / max) |
|---|---|---|---|
| APP | 0.669 | 153 / 31260 / 101774 | 11228 / 11243 / 17076 |
| PSEN1 | 0.622 | 6 / 61728 / 87019 | 34528 / 56324 / 82786 |
| PSEN2 | 0.951 | 415 / 5344 / 11877 | 11704 / 15362 / 19951 |

## Variant type by label

| Type | Healthy label | Unhealthy label |
|---|---|---|
| SNV | 106 | 99 |
| deletion | 2 | 6 |
| insertion | 2 | 2 |
| multi_SNV | 0 | 3 |

## Leakage in the original split

On the original seed-42 per-sequence split, the nearest training variant gets 33/46 test labels right (balanced accuracy 0.717): neighbouring variants fall on both sides.

## Limitations

- The labels have no recorded source, so they cannot be treated as ground truth about biology.
- Labels are confounded with position, and variant type is unevenly distributed between them. Any signal may reflect how the variants were selected.
- Few independent regions exist: variants cluster in a small number of windows, and PSEN2 has only 18 variants, so per-gene intervals are wide.
- Sequence context here is minimal; no exon, consequence, conservation or allele-frequency annotation is available.
- Each sequence is the reference plus one inserted variant, not a real haplotype, so nothing generalises to real genomes or individuals.
