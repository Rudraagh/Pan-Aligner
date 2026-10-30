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
python main.py --ml-alignment                            # Model A training and evaluation
python main.py --ml-alignment --output-dir outputs/ml_alignment/sensitivity_strict --max-accept-error 0.001 --min-reject-precision 0.99
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

**Evaluation** (`scripts/ml/alignment_reliability.py`):

- GroupKFold(5) over source haplotypes.
- Inside each training fold, a haplotype-disjoint quarter calibrates the probabilities (isotonic, `CalibratedClassifierCV(FrozenEstimator(model))`) and sets two thresholds:
  - accept as `correct` while error among accepted stays within the target;
  - reject as `incorrect` while rejection precision stays at or above the target;
  - anything between them is `insufficient_evidence`, an abstention band rather than a trained class.
- A MAPQ-only model calibrated the same way is the baseline.

**Results** (`outputs/ml_alignment/REPORT.md`, `model_results.json`):

| Model | PR-AUC (incorrect) | ROC-AUC | Brier | ECE |
|---|---|---|---|---|
| MAPQ only | 0.569 | 0.813 | 0.0572 | 0.0076 |
| Logistic regression | 0.994 | 0.999 | 0.0047 | 0.0016 |
| Gradient boosting | 0.995 | 0.999 | 0.0033 | 0.0015 |

| Targets (error among accepted / rejection precision) | Model | Coverage | Held-out error among accepted | Incorrect caught | Within-locus caught (of 96) |
|---|---|---|---|---|---|
| 1% / 90% | MAPQ only | 0.1% | – | 0.9% | 1 |
| | Logistic | 100% | 0.38% | 96.9% | 89 |
| | Gradient boosting | 100% | 0.31% | 97.5% | 96 |
| 0.1% / 99% (`sensitivity_strict/`) | MAPQ only | 0% | – | 0% | 0 |
| | Logistic | 97.9% | 0.17% | 92.6% | 81 |
| | Gradient boosting | 98.2% | 0.12% | 97.2% | 96 |

What the results show:

- **Defaults:** with 1% / 90% the accept and reject thresholds overlapped in all five folds, so `insufficient_evidence` was never used and the decision was binary.
- **Stricter targets:** with 0.1% / 99% the abstention band opened for about 2% of alignments. The held-out error (0.12–0.17%) was still above the 0.1% target, so thresholds set on calibration data did not fully carry over to unseen haplotypes.
- **Within-locus false alarms:** gradient boosting caught all 96 within-locus misplacements under both targets, but also flagged 136 (defaults) and 112 (strict) correct within-locus alignments as incorrect.
- **MAPQ alone met neither target**, so it abstained on almost everything.
- **PanAligner's off-target alignments look confident:** in the features file (regenerable, not committed) they have median MAPQ 60 and median identity 0.19. Minigraph's are short, about 176 bp with MAPQ around 20. This is why MAPQ alone is a weak reliability signal here.

**Limits:**

- The incorrect class is almost entirely off-target reads on single-locus graphs, where the read's true locus is not available to the aligner. The models mainly learn that low identity and coverage mean wrong. This is evidence for detecting off-target or contaminating reads, not for whole-genome alignment.
- Within-locus misplacements number only 96, too few for strong conclusions about them.
- The reads are simulated at locus scale; real sequencing data and whole-genome graphs are not covered.

## Model B: feasibility study (research label recovery, not diagnosis)

**The data:** each of the 220 labelled sequences (110 unhealthy) is its locus's GRCh38 reference with one inserted variant: 205 SNVs, 3 multi-SNV, 4 insertions, 8 deletions. **The source, version and criteria of the healthy/unhealthy labels are unknown.** Results therefore describe recovery of these labels, never disease association.

**The study** (`scripts/ml/variant_feasibility.py`, results at commit `7caa878`, generated by the code at `d3cf31a`):

- **Features:** each example is reduced to its variant: gene, relative position, SNV/indel, length change, transition/transversion, local GC (±50 bp) and CpG context.
- **Splits:** StratifiedGroupKFold(5) over 1 kb genomic windows, repeated with seeds 0–4.
- **"Uncertain":** logistic models answer `uncertain` when |p − 0.5| is below a margin chosen by inner grouped CV for 80% accuracy.

| Model | Balanced accuracy (mean ± SD) | MCC | Confident predictions (mean of 220) |
|---|---|---|---|
| Majority class | 0.481 ± 0.002 | −0.039 | – |
| Nearest training variant (same gene) | 0.628 ± 0.009 | 0.256 | – |
| Position-only logistic | 0.500 ± 0.026 | 0.000 | 0.0 |
| Context-only logistic | 0.550 ± 0.015 | 0.103 | 0.2 |
| Context + position logistic | 0.555 ± 0.009 | 0.109 | 0.4 |

**Conclusion:**

- **No model recovers the labels from sequence context.** The best signal is the nearest training variant, which reflects where the labelled variants were chosen: they cluster by position, and some clusters are wider than 1 kb.
- **The models abstain on almost everything:** on average 0–0.4 of 220 predictions are confident.
- **The original random split leaks.** On the seed-42 per-sequence split, the nearest-training-variant rule gets 33 of 46 test labels right (balanced accuracy 0.717), because neighbouring variants fall on both sides.

A meaningful Model B would need a documented label source and variant-level annotations (exon/intron, codon consequence, conservation, allele frequency), none of which are in the repository.
