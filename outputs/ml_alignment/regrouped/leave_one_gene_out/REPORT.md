# Model A: alignment reliability

Labels come from simulation truth (PBSIM2 reads from the held-out haplotypes), never from the aligners being judged.

Rows: 113295 correct, 14073 incorrect, 259920 unaligned (not used), 0 unprojectable (not used).

13977 of 14073 incorrect alignments are off-target: a read aligned to another gene's graph, where its true locus was not available to the aligner. The rest are within-locus misplacements.

Evaluation: each gene's reads tested in turn, with that gene's reads and graph removed from training; isotonic calibration and thresholds on a group-disjoint part of each training fold. Targets: error among accepted <= 0.01, rejection precision >= 0.9.

| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|---|---|---|---|
| mapq_only | 0.549 | 0.789 | 0.0623 | 0.0478 | 0.437 | 0.0679 | 0.000 |
| logistic | 0.956 | 0.982 | 0.0111 | 0.0097 | 1.000 | 0.0109 | 0.911 |
| gradient_boosting | 0.929 | 0.957 | 0.0116 | 0.0102 | 1.000 | 0.0132 | 0.893 |

## Within-locus and off-target errors separately

Both columns share the same correct rows, but within-locus incorrect rows are far rarer, so compare each column with its own prevalence.

| Model | Within-locus PR-AUC | Within-locus ROC-AUC | Off-target PR-AUC | Off-target ROC-AUC |
|---|---|---|---|---|
| mapq_only | 0.038 | 0.942 | 0.546 | 0.788 |
| logistic | 0.261 | 0.998 | 0.956 | 0.982 |
| gradient_boosting | 0.502 | 1.000 | 0.928 | 0.957 |

## Plain MAPQ thresholds (no training)

Accept when MAPQ >= Q. `gradient_boosting at same coverage` is that model's error among the same share of alignments, chosen by lowest p(incorrect).

| Q | Coverage | Error among accepted | Incorrect caught | gradient_boosting at same coverage |
|---|---|---|---|---|
| 0 | 1.000 | 0.1105 | 0.000 | 0.1105 |
| 1 | 0.997 | 0.1087 | 0.019 | 0.1086 |
| 5 | 0.989 | 0.1026 | 0.082 | 0.1011 |
| 10 | 0.979 | 0.0953 | 0.156 | 0.0918 |
| 20 | 0.964 | 0.0859 | 0.250 | 0.0782 |
| 30 | 0.949 | 0.0753 | 0.353 | 0.0636 |
| 40 | 0.934 | 0.0644 | 0.455 | 0.0485 |
| 50 | 0.922 | 0.0556 | 0.537 | 0.0354 |
| 60 | 0.896 | 0.0424 | 0.657 | 0.0104 |

`minigraph` only:

| Q | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|
| 0 | 1.000 | 0.1100 | 0.000 |
| 1 | 0.995 | 0.1064 | 0.038 |
| 5 | 0.979 | 0.0942 | 0.162 |
| 10 | 0.958 | 0.0795 | 0.307 |
| 20 | 0.930 | 0.0599 | 0.493 |
| 30 | 0.905 | 0.0423 | 0.651 |
| 40 | 0.888 | 0.0312 | 0.748 |
| 50 | 0.872 | 0.0219 | 0.826 |
| 60 | 0.861 | 0.0160 | 0.875 |

`panaligner` only:

| Q | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|
| 0 | 1.000 | 0.1110 | 0.000 |
| 1 | 1.000 | 0.1110 | 0.000 |
| 5 | 1.000 | 0.1110 | 0.000 |
| 10 | 1.000 | 0.1109 | 0.001 |
| 20 | 1.000 | 0.1109 | 0.001 |
| 30 | 0.995 | 0.1063 | 0.048 |
| 40 | 0.982 | 0.0953 | 0.156 |
| 50 | 0.972 | 0.0867 | 0.241 |
| 60 | 0.932 | 0.0675 | 0.434 |

## Notes

- `mapq_only` on within-locus alignments: 96 incorrect of 113391; caught 0.000, 11 accepted as correct, 0 flagged incorrect in total.
- `logistic`: abstention band collapsed in 3 of 3 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `logistic` on within-locus alignments: 96 incorrect of 113391; caught 0.896, 10 accepted as correct, 544 flagged incorrect in total.
- `gradient_boosting`: abstention band collapsed in 3 of 3 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `gradient_boosting` on within-locus alignments: 96 incorrect of 113391; caught 1.000, 0 accepted as correct, 513 flagged incorrect in total.

## Limits

- The incorrect class is dominated by off-target reads on single-locus graphs, so the models mainly learn that low identity and coverage mean wrong. This is evidence about detecting off-target or contaminating reads, not about whole-genome alignment.
- Simulated reads (PBSIM2, locus scale); real sequencing data and whole-genome graphs are not covered.
