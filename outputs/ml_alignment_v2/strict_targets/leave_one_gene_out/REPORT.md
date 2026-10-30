# Model A strict targets: leave_one_gene_out

Grouping: each gene's reads tested in turn, with that gene's reads and graph removed from training. Thresholds come only from a validation part of each training fold, disjoint by group from the rows that fitted the model, the rows that fitted the calibrator, and the test fold.

Held-out calibration of gradient boosting: ECE 0.0026, Brier 0.0028, mean predicted 0.67% against observed 0.83%.

## Targets and abstention

A band is meaningful only if it opens in every fold, covers at least 0.1% of alignments, and its error rate lies between accepted and rejected.

| Target | Model | Accepted | Abstentions | Rejected | Error among accepted | Error among abstained | Precision of rejections | Folds with band | Meaningful |
|---|---|---|---|---|---|---|---|---|---|
| error_0.5pct | gradient_boosting | 190746 | 0 | 26 | 0.819% | – | 96.2% | 2/3 | no |
| error_0.5pct | mapq_only | 190772 | 0 | 0 | 0.832% | – | – | 3/3 | no |
| error_0.2pct | gradient_boosting | 189651 | 1095 | 26 | 0.275% | 95.16% | 96.2% | 3/3 | yes |
| error_0.2pct | mapq_only | 186888 | 3884 | 0 | 0.571% | 13.41% | – | 3/3 | yes |
| error_0.1pct | gradient_boosting | 189262 | 1484 | 26 | 0.205% | 79.18% | 96.2% | 3/3 | yes |
| error_0.1pct | mapq_only | 178932 | 11840 | 0 | 0.265% | 9.40% | – | 3/3 | yes |
| error_0.05pct | gradient_boosting | 101630 | 89116 | 26 | 0.286% | 1.43% | 96.2% | 2/3 | no |
| error_0.05pct | mapq_only | 166834 | 23938 | 0 | 0.161% | 5.51% | – | 3/3 | yes |
| coverage_95pct | gradient_boosting | 189161 | 1585 | 26 | 0.187% | 76.34% | 96.2% | 3/3 | yes |
| coverage_95pct | mapq_only | 183226 | 7546 | 0 | 0.403% | 11.25% | – | 3/3 | yes |
| coverage_90pct | gradient_boosting | 189161 | 1585 | 26 | 0.187% | 76.34% | 96.2% | 3/3 | yes |
| coverage_90pct | mapq_only | 175541 | 15231 | 0 | 0.203% | 8.09% | – | 3/3 | yes |

## Within-locus and off-target separately (gradient boosting)

| Target | Within-locus errors accepted / abstained / rejected | Within-locus correct abstained / rejected | Off-target errors accepted / abstained / rejected |
|---|---|---|---|
| error_0.5pct | 523 / 0 / 13 | 0 / 1 | 1040 / 0 / 12 |
| error_0.2pct | 67 / 456 / 13 | 53 / 1 | 454 / 586 / 12 |
| error_0.1pct | 5 / 518 / 13 | 309 / 1 | 383 / 657 / 12 |
| error_0.05pct | 5 / 518 / 13 | 87844 / 1 | 286 / 754 / 12 |
| coverage_95pct | 1 / 522 / 13 | 375 / 1 | 352 / 688 / 12 |
| coverage_90pct | 1 / 522 / 13 | 375 / 1 | 352 / 688 / 12 |

## Gradient boosting against plain MAPQ thresholds at equal coverage

Difference is model minus MAPQ error, with a 95% interval resampling whole 20 kb windows. Accepted errors are split into within-locus and off-target.

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) | MAPQ accepted errors: within / off-target | Model accepted errors: within / off-target |
|---|---|---|---|---|---|---|
| 10 | 98.0% | 0.561% | 0.181% | -0.00380 (-0.00538 to -0.00246) | 355 / 693 | 1 / 337 |
| 20 | 96.2% | 0.404% | 0.171% | -0.00233 (-0.00374 to -0.00128) | 231 / 510 | 1 / 313 |
| 30 | 94.6% | 0.295% | 0.167% | -0.00129 (-0.00245 to -0.00053) | 147 / 386 | 1 / 299 |
| 40 | 93.0% | 0.242% | 0.167% | -0.00075 (-0.00153 to -0.00020) | 104 / 325 | 1 / 295 |
| 50 | 91.6% | 0.191% | 0.167% | -0.00025 (-0.00100 to +0.00031) | 70 / 264 | 1 / 290 |
| 60 | 86.7% | 0.128% | 0.167% | +0.00038 (-0.00030 to +0.00089) | 35 / 177 | 1 / 275 |

**Gradient boosting does not have reliably lower error than the plain MAPQ threshold, not at every compared coverage, so these results do not support using it instead of MAPQ.**

`minigraph` only:

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) |
|---|---|---|---|---|
| 10 | 96.3% | 0.857% | 0.199% | -0.00893 to -0.00448 |
| 20 | 92.8% | 0.565% | 0.187% | -0.00590 to -0.00223 |
| 30 | 89.8% | 0.383% | 0.187% | -0.00339 to -0.00099 |
| 40 | 87.0% | 0.289% | 0.187% | -0.00197 to -0.00036 |
| 50 | 84.6% | 0.224% | 0.187% | -0.00110 to +0.00011 |
| 60 | 82.2% | 0.173% | 0.187% | -0.00052 to +0.00069 |

`panaligner` only:

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) |
|---|---|---|---|---|
| 10 | 100.0% | 0.236% | 0.236% | -0.00003 to +0.00002 |
| 20 | 100.0% | 0.234% | 0.233% | -0.00007 to +0.00002 |
| 30 | 100.0% | 0.206% | 0.203% | -0.00015 to +0.00012 |
| 40 | 99.9% | 0.195% | 0.170% | -0.00054 to +0.00020 |
| 50 | 99.6% | 0.160% | 0.153% | -0.00122 to +0.00078 |
| 60 | 91.8% | 0.083% | 0.143% | -0.00011 to +0.00125 |

## Limits

- Simulated PBSIM2 reads at 75% accuracy on three loci; real data and whole-genome graphs are not covered.
- Errors from the same regions repeat, so the effective sample is smaller than the row counts, and there are few independent windows.
- Rankings pool probabilities from different folds' calibrators.
- Targets are compared across this one dataset; no separate untouched test set was held back for choosing among them.
