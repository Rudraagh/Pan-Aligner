# Model A: alignment reliability

> **Weak split.** These results use haplotype folds, which put the same locus on both sides because every haplotype is the same reference plus one variant. They are optimistic. See `outputs/ml_alignment/regrouped/` and `docs/17_ml_layer.md`.

Labels come from simulation truth (PBSIM2 reads from the held-out haplotypes), never from the aligners being judged.

Rows: 113295 correct, 14073 incorrect, 259920 unaligned (not used), 0 unprojectable (not used).

13977 of 14073 incorrect alignments are off-target: a read aligned to another gene's graph, where its true locus was not available to the aligner. The rest are within-locus misplacements.

Evaluation: GroupKFold over source haplotypes; isotonic calibration and thresholds on a haplotype-disjoint part of each training fold. Targets: error among accepted <= 0.01, rejection precision >= 0.9.

| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|---|---|---|---|
| mapq_only | 0.569 | 0.813 | 0.0572 | 0.0076 | 0.001 | – | 0.009 |
| logistic | 0.994 | 0.999 | 0.0047 | 0.0016 | 1.000 | 0.0038 | 0.969 |
| gradient_boosting | 0.995 | 0.999 | 0.0033 | 0.0015 | 1.000 | 0.0031 | 0.975 |

## Notes

- `mapq_only` on within-locus alignments: 96 incorrect of 113391; caught 0.010, 0 accepted as correct, 26 flagged incorrect in total.
- `logistic`: abstention band collapsed in 5 of 5 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `logistic` on within-locus alignments: 96 incorrect of 113391; caught 0.927, 7 accepted as correct, 349 flagged incorrect in total.
- `gradient_boosting`: abstention band collapsed in 5 of 5 folds (the accept and reject targets overlapped), so `insufficient_evidence` was not used there.
- `gradient_boosting` on within-locus alignments: 96 incorrect of 113391; caught 1.000, 0 accepted as correct, 232 flagged incorrect in total.

## Limits

- The incorrect class is dominated by off-target reads on single-locus graphs, so the models mainly learn that low identity and coverage mean wrong. This is evidence about detecting off-target or contaminating reads, not about whole-genome alignment.
- Simulated reads (PBSIM2, locus scale); real sequencing data and whole-genome graphs are not covered.
