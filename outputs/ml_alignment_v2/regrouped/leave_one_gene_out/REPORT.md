# Model A: alignment reliability

Labels come from simulation truth (PBSIM2 reads from the held-out haplotypes), never from the aligners being judged.

Rows: 189184 correct, 1588 incorrect, 664204 unaligned (not used), 0 unprojectable (not used).

1052 of 1588 incorrect alignments are off-target: a read aligned to another gene's graph, where its true locus was not available to the aligner. The rest are within-locus misplacements.

Evaluation: each gene's reads tested in turn, with that gene's reads and graph removed from training; isotonic calibration and thresholds on a group-disjoint part of each training fold. Targets: error among accepted <= 0.01, rejection precision >= 0.9.

| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|---|---|---|---|
| mapq_only | 0.105 | 0.908 | 0.0079 | 0.0056 | 1.000 | 0.0083 | 0.000 |
| logistic | 0.701 | 0.927 | 0.0035 | 0.0025 | 1.000 | 0.0053 | 0.370 |
| gradient_boosting | 0.758 | 0.919 | 0.0027 | 0.0025 | 1.000 | 0.0028 | 0.666 |

## Within-locus and off-target errors separately

Both columns share the same correct rows, but within-locus incorrect rows are far rarer, so compare each column with its own prevalence.

| Model | Within-locus PR-AUC | Within-locus ROC-AUC | Off-target PR-AUC | Off-target ROC-AUC |
|---|---|---|---|---|
| mapq_only | 0.038 | 0.915 | 0.074 | 0.905 |
| logistic | 0.525 | 0.996 | 0.631 | 0.892 |
| gradient_boosting | 0.946 | 1.000 | 0.632 | 0.878 |

## Plain MAPQ thresholds (no training)

Accept when MAPQ >= Q. `gradient_boosting at same coverage` is that model's error among the same share of alignments, chosen by lowest p(incorrect).

| Q | Coverage | Error among accepted | Incorrect caught | gradient_boosting at same coverage |
|---|---|---|---|---|
| 0 | 1.000 | 0.0083 | 0.000 | 0.0083 |
| 1 | 0.998 | 0.0079 | 0.051 | 0.0061 |
| 5 | 0.990 | 0.0068 | 0.195 | 0.0021 |
| 10 | 0.980 | 0.0056 | 0.340 | 0.0017 |
| 20 | 0.962 | 0.0040 | 0.533 | 0.0018 |
| 30 | 0.946 | 0.0030 | 0.664 | 0.0018 |
| 40 | 0.930 | 0.0024 | 0.730 | 0.0015 |
| 50 | 0.916 | 0.0019 | 0.790 | 0.0014 |
| 60 | 0.867 | 0.0013 | 0.866 | 0.0013 |

`minigraph` only:

| Q | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|
| 0 | 1.000 | 0.0135 | 0.000 |
| 1 | 0.996 | 0.0128 | 0.059 |
| 5 | 0.981 | 0.0107 | 0.222 |
| 10 | 0.963 | 0.0086 | 0.390 |
| 20 | 0.928 | 0.0056 | 0.612 |
| 30 | 0.898 | 0.0038 | 0.746 |
| 40 | 0.870 | 0.0029 | 0.814 |
| 50 | 0.846 | 0.0022 | 0.860 |
| 60 | 0.822 | 0.0017 | 0.895 |

`panaligner` only:

| Q | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|
| 0 | 1.000 | 0.0024 | 0.000 |
| 1 | 1.000 | 0.0024 | 0.000 |
| 5 | 1.000 | 0.0024 | 0.019 |
| 10 | 1.000 | 0.0024 | 0.023 |
| 20 | 1.000 | 0.0023 | 0.032 |
| 30 | 1.000 | 0.0021 | 0.148 |
| 40 | 0.999 | 0.0019 | 0.194 |
| 50 | 0.996 | 0.0016 | 0.343 |
| 60 | 0.918 | 0.0008 | 0.685 |

## Notes

- `mapq_only` on within-locus alignments: 536 incorrect of 189720; caught 0.000, 536 accepted as correct, 0 flagged incorrect in total.
- `logistic`: abstention band collapsed in 3 of 3 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `logistic` on within-locus alignments: 536 incorrect of 189720; caught 0.289, 381 accepted as correct, 288 flagged incorrect in total.
- `gradient_boosting`: abstention band collapsed in 3 of 3 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `gradient_boosting` on within-locus alignments: 536 incorrect of 189720; caught 0.890, 59 accepted as correct, 493 flagged incorrect in total.

## Limits

- 66% of incorrect alignments are off-target reads on single-locus graphs, where the true locus was not available to the aligner. This says nothing direct about whole-genome alignment.
- Of 536 within-locus errors, 96% from minigraph, so those results mostly describe that aligner. Errors from the same regions repeat, so the effective sample is smaller than the count.
- Simulated reads (PBSIM2, locus scale); real sequencing data and whole-genome graphs are not covered.
