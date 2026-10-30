# 17. Machine-learning layer

The project has two separate ML tasks. Both live in `scripts/ml/` and are reachable through `main.py`. Neither is a clinical tool.

| Task | Question | Status |
|---|---|---|
| Model A, alignment reliability | Is this alignment correct, incorrect, or is there insufficient evidence? | Trained and evaluated on simulated reads with independent truth |
| Model B, healthy/disease-associated classification | Can the healthy/unhealthy labels be recovered from the variant? | Feasibility study only: the labels' source is unknown and no model recovers them |

The older healthy/unhealthy ML path (`scripts/pipeline.py`, `extract_features.py`, `train_ml.py`, `predict_ml.py`, `evaluate.py`) is legacy and leaks the label; see `docs/13_current_version_implementation_status.md`. It is not part of this layer.

## How to run

```bash
pip install -r requirements-ml.txt                       # numpy, scikit-learn >= 1.6
python3 main.py --ml-alignment-data --threads 8          # Linux/WSL only: simulate, align, label (~6 min)
python main.py --ml-alignment                            # Model A, 20 kb window folds (writes outputs/ml_alignment/regrouped/window/)
python main.py --ml-alignment --grouping leave_one_gene_out
python main.py --ml-alignment --grouping haplotype --output-dir outputs/ml_alignment/sensitivity_strict --max-accept-error 0.001 --min-reject-precision 0.99
python main.py --ml-variant                              # Model B feasibility study
```

## Model A: alignment reliability

**Data** (`scripts/ml/build_alignment_dataset.py`, run at commit `cf8c9bd`):

- PBSIM2 reads (P6C4 model) from the 46 held-out haplotypes, which are not in any graph, at read accuracy 0.95 / 0.85 / 0.75, read length 10,000 / 2,000 bp, depth 2×, seeds 1–3.
- Each read set is aligned by PanAligner v1.0 and minigraph 0.21-r606 to all three genes' combined train graphs: its own gene's graph and the other two as decoys.
- One row per read and graph: the longest alignment (the paftools mapeval rule), or `unaligned`.

**Labels** come only from simulation truth:

- Own-gene alignments are `correct` when they overlap the true interval (overlap/union > 0.1, same strand), after projecting paths through non-reference nodes onto the reference. Otherwise they are `incorrect`, or `unprojectable` when no projection exists.
- Every alignment to another gene's graph is `incorrect`: the read comes from a different chromosome.
- Strict mapeval "wrong" is never used, because on these graphs it marks correct alignments that pass a non-reference node.

| Label | Rows |
|---|---|
| correct | 113,295 |
| incorrect | 14,073 (13,977 off-target decoys, 96 within-locus misplacements) |
| unaligned (not used) | 259,920 |

**Features** are only what is visible after alignment: MAPQ, query length and coverage, matches, block length, identity, NM rate, chain minimizers (`cm`), best and second chain scores (`s1`, `s2`) and their ratio, divergence (`dv`), path segments, off-reference path, number of alignments for the read, and which aligner produced it. Simulation settings, the decoy flag and the genes are never features.

**Evaluation** (`scripts/ml/alignment_reliability.py`), with three groupings:

- `window` (default): folds over 20 kb windows of the read's true position, with reads from adjacent windows removed from training.
- `leave_one_gene_out`: test on one gene's reads, training without that gene's reads or graph. This is the hardest test.
- `haplotype`: folds over source haplotypes. **This is the grouping behind the first results, and it is too weak:** every haplotype is the same reference plus one variant, so the same locus is on both sides of a split. It is kept only for comparison.

Inside each training fold, a group-disjoint quarter calibrates the probabilities (isotonic, `CalibratedClassifierCV(FrozenEstimator(model))`) and sets two thresholds:

- accept as `correct` while error among accepted stays within the target;
- reject as `incorrect` while rejection precision stays at or above the target;
- anything between them is `insufficient_evidence`, an abstention band rather than a trained class.

There are two baselines: a MAPQ-only model calibrated the same way, and a plain MAPQ threshold sweep with no training.

**Results by grouping**, at the default targets (1% error among accepted, 90% rejection precision). Within-locus PR-AUC is measured against a prevalence of about 0.085%.

| Grouping | Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE | Within-locus PR-AUC | Error among accepted | Incorrect caught |
|---|---|---|---|---|---|---|---|---|
| Window (`regrouped/window/`) | MAPQ only | 0.564 | 0.803 | 0.0575 | 0.0081 | 0.037 | – | 0% |
| | Logistic | 0.986 | 0.996 | 0.0076 | 0.0055 | 0.465 | 0.57% | 95.4% |
| | Gradient boosting | 0.979 | 0.994 | 0.0086 | 0.0068 | 0.881 | 0.69% | 94.4% |
| Leave one gene out (`regrouped/leave_one_gene_out/`) | MAPQ only | 0.549 | 0.789 | 0.0623 | 0.0478 | 0.038 | 6.79% | 0% |
| | Logistic | 0.956 | 0.982 | 0.0111 | 0.0097 | 0.261 | 1.09% | 91.1% |
| | Gradient boosting | 0.929 | 0.957 | 0.0116 | 0.0102 | 0.502 | 1.32% | 89.3% |
| Haplotype, weak (`outputs/ml_alignment/`) | MAPQ only | 0.569 | 0.813 | 0.0572 | 0.0076 | 0.037 | – | 0.9% |
| | Logistic | 0.994 | 0.999 | 0.0047 | 0.0016 | 0.588 | 0.38% | 96.9% |
| | Gradient boosting | 0.995 | 0.999 | 0.0033 | 0.0015 | 0.856 | 0.31% | 97.5% |

**Plain MAPQ thresholds** (no training, so the same under every grouping). Accepting only MAPQ 60:

| Aligner | Coverage | Error among accepted | Incorrect caught |
|---|---|---|---|
| Both | 89.6% | 4.24% | 65.7% |
| PanAligner | 93.2% | 6.75% | 43.4% |
| Minigraph | 86.1% | 1.60% | 87.5% |

At that same 89.6% coverage, gradient boosting's error is 0.84% (window) and 1.04% (leave one gene out), against 4.24% for MAPQ 60.

What the results show:

- **The first 0.995 was optimistic.** Stricter grouping lowers gradient boosting's PR-AUC to 0.979 (window) and 0.929 (leave one gene out). The models still clearly beat MAPQ.
- **The hardest test misses the 1% target.** Leaving a gene out gives 1.09–1.32% error among accepted, and the MAPQ-only model's 6.79% shows the calibration thresholds don't transfer across genes either.
- **Logistic regression is more stable overall.** It beats gradient boosting on overall PR-AUC under both stricter groupings, while gradient boosting is ahead on within-locus PR-AUC.
- **False alarms on within-locus reads rise.** Gradient boosting still catches all 96 within-locus misplacements, but flags 619 correct within-locus alignments as incorrect under window folds and 417 when leaving a gene out, up from 136.
- **The abstention band barely opens at the defaults.** The accept and reject thresholds overlap in all folds for logistic regression, and in all but one window fold for gradient boosting. The stricter 0.1% / 99% targets (`sensitivity_strict/`) opened it for about 2% of alignments, but that run used only the weak haplotype grouping, with held-out error of 0.12–0.17%.
- **PanAligner's off-target alignments look confident:** in the features file (regenerable, not committed) they have median MAPQ 60 and median identity 0.19. Minigraph's are short, about 176 bp with MAPQ around 20. That is why MAPQ is especially weak for PanAligner.

**Limits:**

- The incorrect class is almost entirely off-target reads on single-locus graphs, where the read's true locus is not available to the aligner. The models mainly learn that low identity and coverage mean wrong. This is evidence for detecting off-target or contaminating reads, not for whole-genome alignment.
- Within-locus misplacements number only 96, too few for strong conclusions about them. The focused v2 dataset below addresses this.
- The reads are simulated at locus scale; real sequencing data and whole-genome graphs are not covered.

### Focused v2 dataset (more within-locus errors)

A pilot (depth 0.1, seed 11; scratchpad only, not kept) showed the first proposed v2 grid (accuracy 0.90/0.80/0.70, reads of 5,000/1,000/500 bp) would give roughly 120 within-locus errors, too few, and that most 70% reads do not align. So v2 uses the settings that produced the most errors in the first dataset:

```bash
python3 main.py --ml-alignment-data --accuracies 0.75 --lengths 10000 2000 --depth 10 --seeds 11 12 13 14 --threads 8 --output-dir outputs/ml_alignment_v2
python main.py --ml-alignment --features outputs/ml_alignment_v2/alignment_features.csv.gz --grouping window --output-dir outputs/ml_alignment_v2/regrouped/window
python main.py --ml-alignment --features outputs/ml_alignment_v2/alignment_features.csv.gz --grouping leave_one_gene_out --output-dir outputs/ml_alignment_v2/regrouped/leave_one_gene_out
```

**Dataset** (commit `0c22ced`, 10 minutes): 854,976 rows, of which 189,184 correct, 1,588 incorrect and 664,204 unaligned. It passed the go/no-go check set beforehand (at least 500 within-locus errors, 15 windows, all genes): 536 within-locus errors across 23 windows, with 1,052 off-target. Within-locus errors split PSEN1 364, APP 156, PSEN2 16, and minigraph 517, PanAligner 19.

**Results.** Incorrect alignments are 0.83% of rows here against 11% in the first dataset, so these numbers cannot be compared directly with those above.

| Grouping | Model | PR-AUC, all errors | Within-locus PR-AUC | Off-target PR-AUC | Error among accepted | All errors caught | Within-locus caught (of 536) | Within-locus false alarms |
|---|---|---|---|---|---|---|---|---|
| Window | MAPQ only | 0.086 | 0.035 | 0.056 | 0.76% | 0% | 0 | 0 |
| | Logistic | 0.783 | 0.751 | 0.697 | 0.37% | 55.3% | 350 | 91 |
| | Gradient boosting | 0.821 | 0.965 | 0.724 | 0.25% | 70.2% | 482 | 37 |
| Leave one gene out | MAPQ only | 0.105 | 0.038 | 0.074 | 0.83% | 0% | 0 | 0 |
| | Logistic | 0.701 | 0.525 | 0.631 | 0.53% | 37.0% | 155 | 133 |
| | Gradient boosting | 0.758 | 0.946 | 0.632 | 0.28% | 66.6% | 477 | 16 |

| Plain MAPQ rule | Accepted | Error among accepted | Errors caught | Gradient boosting error at same share (window / gene held out) |
|---|---|---|---|---|
| MAPQ >= 30 | 94.6% | 0.30% | 66.4% | 0.13% / 0.18% |
| MAPQ = 60 | 86.7% | 0.13% | 86.6% | 0.14% / 0.13% |

What v2 shows:

- **Within-locus detection works on simulated reads:** gradient boosting catches about 90% of these errors with few false alarms, including with a whole gene held out. Logistic regression does not hold up.
- **A plain MAPQ 60 cutoff is about as good at its own operating point.** The model's gain appears when fewer alignments are discarded: at 94.6% accepted, its error is well below MAPQ's. So the model is not a clear replacement for MAPQ.
- **The MAPQ-only model's low PR-AUC overstates MAPQ's weakness,** because most alignments tie at MAPQ 60. The plain threshold is the fairer baseline.
- **Off-target detection is harder at 75% accuracy** (PR-AUC 0.63–0.72).
- **Abstention still does not work at the default targets:** the band collapsed in every fold for both learned models, so `insufficient_evidence` was never used. Stricter targets have not been run on v2.

**Limits:** within-locus errors are 96% minigraph's, so PanAligner cannot be judged on them. There is one accuracy level (75%). Errors from the same regions repeat, so the effective sample is far smaller than 536, and PSEN2 has only 16.

## Model B: feasibility study (research label recovery, not diagnosis)

**The data:** each of the 220 labelled sequences (110 unhealthy) is its locus's GRCh38 reference with one inserted variant: 205 SNVs, 3 multi-SNV, 4 insertions, 8 deletions. **The source, version and criteria of the healthy/unhealthy labels are unknown.** Results therefore describe recovery of these labels, never disease association.

**The study** (`scripts/ml/variant_feasibility.py`, full results in `outputs/ml_variant/REPORT.md`):

- **Features:** each example is reduced to its variant: gene, relative position, SNV/indel, length change, transition/transversion, local GC (±50 bp) and CpG context.
- **Splits:** StratifiedGroupKFold(5) over genomic windows (1 kb, with 5 and 20 kb as a sensitivity check), repeated with seeds 0–4.
- **Interval:** 95% interval resampling whole windows, using the first seed's predictions.
- **"Uncertain":** logistic models answer `uncertain` when |p − 0.5| is below a margin chosen by inner grouped CV for 80% accuracy. Their outputs are named `healthy_label` and `unhealthy_label`, because they only predict these labels.
- **Context models** are confound diagnostics, not a classifier. No scientifically valid classifier is claimed.

Class balance is even within each gene: APP 23/23, PSEN1 78/78, PSEN2 9/9.

| Model (1 kb windows) | Balanced accuracy (mean ± SD) | 95% interval | MCC | Confident predictions (mean of 220) |
|---|---|---|---|---|
| Majority class | 0.481 ± 0.002 | 0.359–0.604 | −0.039 | – |
| Nearest training variant (same gene) | 0.628 ± 0.009 | 0.505–0.752 | 0.256 | – |
| Position vote (5 nearest, same gene) | 0.646 ± 0.014 | 0.551–0.761 | 0.338 | – |
| Position-only logistic (linear) | 0.500 ± 0.026 | 0.406–0.645 | 0.000 | 0.0 |
| Context-only logistic | 0.550 ± 0.015 | 0.502–0.649 | 0.103 | 0.2 |
| Context + position logistic | 0.555 ± 0.009 | 0.483–0.648 | 0.109 | 0.4 |

| Gene: balanced accuracy | APP | PSEN1 | PSEN2 |
|---|---|---|---|
| Nearest training variant | 0.783 | 0.577 | 0.678 |
| Position vote | 0.696 | 0.632 | 0.644 |
| Context-only logistic | 0.470 | 0.587 | 0.433 |

| Window kept out | Nearest training variant | Position vote | Context-only logistic |
|---|---|---|---|
| 1 kb | 0.628 | 0.646 | 0.550 |
| 5 kb | 0.570 | 0.602 | 0.528 |
| 20 kb | 0.450 | 0.591 | 0.453 |

**Confounding:**

- **Position:** AUC of position against the label is 0.951 for PSEN2, 0.669 for APP and 0.622 for PSEN1. APP's value understates the clustering, because AUC only measures a one-directional split: its unhealthy-label variants all lie at positions 11,228–17,076, inside the healthy-label range of 153–101,774.
- **Variant type:** deletions are 6 unhealthy to 2 healthy, and all 3 multi-SNV sequences are unhealthy; insertions are 2 and 2, and SNVs 99 and 106. Type could therefore carry some of the label, though the counts are small.
- **The original random split leaks.** On the seed-42 per-sequence split, the nearest-training-variant rule gets 33 of 46 test labels right (balanced accuracy 0.717), because neighbouring variants fall on both sides.

**Conclusion and limits:**

- **No model recovers the labels from sequence context.** Context-only reaches 0.550 at 1 kb and 0.453 at 20 kb, and it is at or below chance on APP (0.470) and PSEN2 (0.433), with wide intervals.
- **The strongest signal is position, and it weakens as windows grow.** It reflects where the labelled variants were chosen rather than anything learned about variant effect.
- **The models abstain on almost everything:** on average 0–0.4 of 220 predictions are confident, so "accuracy when confident" carries no weight.
- **Intervals are wide and regions few**, especially for PSEN2's 18 variants.
- **These are not real haplotypes.** Nothing here generalises to real genomes or individuals.

A meaningful Model B would need a documented label source and variant-level annotations (exon/intron, codon consequence, conservation, allele frequency), none of which are in the repository.
