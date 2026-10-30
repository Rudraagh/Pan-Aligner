# Model A: alignment reliability

Labels come from simulation truth (PBSIM2 reads from the held-out haplotypes), never from the aligners being judged.

Rows: 113295 correct, 14073 incorrect, 259920 unaligned (not used), 0 unprojectable (not used).

13977 of 14073 incorrect alignments are off-target: a read aligned to another gene's graph, where its true locus was not available to the aligner. The rest are within-locus misplacements.

Evaluation: GroupKFold over source haplotypes; isotonic calibration and thresholds on a haplotype-disjoint part of each training fold. Targets: error among accepted <= 0.001, rejection precision >= 0.99.

| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|---|---|---|---|
| mapq_only | 0.569 | 0.813 | 0.0572 | 0.0076 | 0.000 | – | 0.000 |
| logistic | 0.994 | 0.999 | 0.0047 | 0.0016 | 0.979 | 0.0017 | 0.926 |
| gradient_boosting | 0.995 | 0.999 | 0.0033 | 0.0015 | 0.982 | 0.0012 | 0.972 |

## Notes

- `mapq_only` on within-locus alignments: 96 incorrect of 113391; caught 0.000, 0 accepted as correct, 0 flagged incorrect in total.
- `logistic` on within-locus alignments: 96 incorrect of 113391; caught 0.844, 2 accepted as correct, 195 flagged incorrect in total.
- `gradient_boosting` on within-locus alignments: 96 incorrect of 113391; caught 1.000, 0 accepted as correct, 208 flagged incorrect in total.

## Limits

- The incorrect class is dominated by off-target reads on single-locus graphs, so the models mainly learn that low identity and coverage mean wrong. This is evidence about detecting off-target or contaminating reads, not about whole-genome alignment.
- Simulated reads (PBSIM2, locus scale); real sequencing data and whole-genome graphs are not covered.
