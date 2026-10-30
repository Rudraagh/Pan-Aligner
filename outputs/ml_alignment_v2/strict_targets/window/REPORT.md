# Model A strict targets: window

Grouping: GroupKFold over 20 kb windows of the read's true position, with adjacent windows removed from training. Thresholds come only from a validation part of each training fold, disjoint by group from the rows that fitted the model, the rows that fitted the calibrator, and the test fold.

Held-out calibration of gradient boosting: ECE 0.0015, Brier 0.0026, mean predicted 0.98% against observed 0.83%.

## Targets and abstention

A band is meaningful only if it opens in every fold, covers at least 0.1% of alignments, and its error rate lies between accepted and rejected.

| Target | Model | Accepted | Abstentions | Rejected | Error among accepted | Error among abstained | Precision of rejections | Folds with band | Meaningful |
|---|---|---|---|---|---|---|---|---|---|
| error_0.5pct | gradient_boosting | 189748 | 0 | 1024 | 0.301% | – | 99.2% | 0/5 | no |
| error_0.5pct | mapq_only | 188412 | 2360 | 0 | 0.733% | 8.73% | – | 5/5 | yes |
| error_0.2pct | gradient_boosting | 189678 | 234 | 860 | 0.292% | 77.78% | 99.1% | 4/5 | no |
| error_0.2pct | mapq_only | 183581 | 7191 | 0 | 0.442% | 10.81% | – | 5/5 | yes |
| error_0.1pct | gradient_boosting | 189596 | 316 | 860 | 0.273% | 68.99% | 99.1% | 4/5 | no |
| error_0.1pct | mapq_only | 177795 | 12977 | 0 | 0.290% | 8.26% | – | 5/5 | yes |
| error_0.05pct | gradient_boosting | 152448 | 37464 | 860 | 0.272% | 0.86% | 99.1% | 5/5 | yes |
| error_0.05pct | mapq_only | 108213 | 82559 | 0 | 0.146% | 1.73% | – | 3/5 | no |
| coverage_95pct | gradient_boosting | 185211 | 4704 | 857 | 0.147% | 9.89% | 99.2% | 5/5 | yes |
| coverage_95pct | mapq_only | 182589 | 8183 | 0 | 0.376% | 11.01% | – | 5/5 | yes |
| coverage_90pct | gradient_boosting | 183291 | 6624 | 857 | 0.147% | 7.07% | 99.2% | 5/5 | yes |
| coverage_90pct | mapq_only | 175881 | 14891 | 0 | 0.221% | 8.05% | – | 5/5 | yes |

## Within-locus and off-target separately (gradient boosting)

| Target | Within-locus errors accepted / abstained / rejected | Within-locus correct abstained / rejected | Off-target errors accepted / abstained / rejected |
|---|---|---|---|
| error_0.5pct | 114 / 0 / 422 | 0 / 8 | 458 / 0 / 594 |
| error_0.2pct | 107 / 85 / 344 | 52 / 8 | 447 / 97 / 508 |
| error_0.1pct | 93 / 99 / 344 | 98 / 8 | 425 / 119 / 508 |
| error_0.05pct | 37 / 155 / 344 | 37143 / 8 | 378 / 166 / 508 |
| coverage_95pct | 3 / 191 / 342 | 4239 / 7 | 270 / 274 / 508 |
| coverage_90pct | 3 / 191 / 342 | 6156 / 7 | 267 / 277 / 508 |

## Gradient boosting against plain MAPQ thresholds at equal coverage

Difference is model minus MAPQ error, with a 95% interval resampling whole 20 kb windows. Accepted errors are split into within-locus and off-target.

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) | MAPQ accepted errors: within / off-target | Model accepted errors: within / off-target |
|---|---|---|---|---|---|---|
| 10 | 98.0% | 0.561% | 0.157% | -0.00403 (-0.00535 to -0.00265) | 355 / 693 | 6 / 288 |
| 20 | 96.2% | 0.404% | 0.155% | -0.00249 (-0.00354 to -0.00141) | 231 / 510 | 6 / 279 |
| 30 | 94.6% | 0.295% | 0.158% | -0.00137 (-0.00205 to -0.00066) | 147 / 386 | 6 / 279 |
| 40 | 93.0% | 0.242% | 0.161% | -0.00081 (-0.00130 to -0.00024) | 104 / 325 | 6 / 279 |
| 50 | 91.6% | 0.191% | 0.163% | -0.00028 (-0.00061 to +0.00016) | 70 / 264 | 6 / 279 |
| 60 | 86.7% | 0.128% | 0.172% | +0.00044 (-0.00013 to +0.00139) | 35 / 177 | 6 / 279 |

**Gradient boosting does not have reliably lower error than the plain MAPQ threshold, not at every compared coverage, so these results do not support using it instead of MAPQ.**

`minigraph` only:

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) |
|---|---|---|---|---|
| 10 | 96.3% | 0.857% | 0.207% | -0.00861 to -0.00447 |
| 20 | 92.8% | 0.565% | 0.215% | -0.00512 to -0.00206 |
| 30 | 89.8% | 0.383% | 0.222% | -0.00249 to -0.00082 |
| 40 | 87.0% | 0.289% | 0.229% | -0.00130 to +0.00006 |
| 50 | 84.6% | 0.224% | 0.235% | -0.00058 to +0.00102 |
| 60 | 82.2% | 0.173% | 0.242% | -0.00021 to +0.00212 |

`panaligner` only:

| MAPQ >= | Accepted | MAPQ error | Model error | Difference (95% interval) |
|---|---|---|---|---|
| 10 | 100.0% | 0.236% | 0.235% | -0.00004 to +0.00000 |
| 20 | 100.0% | 0.234% | 0.231% | -0.00008 to +0.00000 |
| 30 | 100.0% | 0.206% | 0.197% | -0.00019 to -0.00002 |
| 40 | 99.9% | 0.195% | 0.156% | -0.00061 to -0.00004 |
| 50 | 99.6% | 0.160% | 0.120% | -0.00114 to +0.00010 |
| 60 | 91.8% | 0.083% | 0.101% | -0.00008 to +0.00068 |

## Limits

- Simulated PBSIM2 reads at 75% accuracy on three loci; real data and whole-genome graphs are not covered.
- Errors from the same regions repeat, so the effective sample is smaller than the row counts, and there are few independent windows.
- Rankings pool probabilities from different folds' calibrators.
- Targets are compared across this one dataset; no separate untouched test set was held back for choosing among them.
