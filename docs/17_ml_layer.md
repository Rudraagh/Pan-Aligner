# 17. Machine-learning layer

The project has two separate ML tasks. Both live in `scripts/ml/`. Model A and the earlier Model B study are reachable through `main.py`; the current Model B scripts run directly. Neither task is a clinical tool.

| Task | Question | Status |
|---|---|---|
| Model A, alignment reliability | Is this alignment correct, incorrect, or is there insufficient evidence? | Trained and evaluated on simulated reads with independent truth |
| Model B, variant classification experiment | Can sequence-derived features recover ClinVar classification direction in APP, PSEN1 and PSEN2? | Trained once under a fixed protocol. Negative result: the model did not beat simple baselines. Not an Alzheimer's classifier |

The older healthy/unhealthy ML path (`scripts/pipeline.py`, `extract_features.py`, `train_ml.py`, `predict_ml.py`, `evaluate.py`) is legacy and leaks the label; see `docs/13_current_version_implementation_status.md`. It is not part of this layer.

## How to run

```bash
pip install -r requirements-ml.txt                       # numpy, scikit-learn >= 1.6
python3 main.py --ml-alignment-data --threads 8          # Linux/WSL only: simulate, align, label (~6 min)
python main.py --ml-alignment                            # Model A, 20 kb window folds (writes outputs/ml_alignment/regrouped/window/)
python main.py --ml-alignment --grouping leave_one_gene_out
python main.py --ml-alignment --grouping haplotype --output-dir outputs/ml_alignment/sensitivity_strict --max-accept-error 0.001 --min-reject-precision 0.99
python main.py --ml-alignment-strict                     # Model A strict targets on the v2 dataset
python main.py --ml-variant                              # Earlier Model B study (superseded; current Model B is run directly, see below)
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

The right-hand column orders tied probabilities by row. The strict-target section below replaces it with a tie-aware comparison and intervals.

What v2 shows:

- **Within-locus detection works on simulated reads:** gradient boosting catches about 90% of these errors with few false alarms, including with a whole gene held out. Logistic regression does not hold up.
- **A plain MAPQ 60 cutoff is about as good at its own operating point.** The model's gain appears when fewer alignments are discarded. So the model is not a clear replacement for MAPQ.
- **The MAPQ-only model's low PR-AUC overstates MAPQ's weakness,** because most alignments tie at MAPQ 60. The plain threshold is the fairer baseline.
- **Off-target detection is harder at 75% accuracy** (PR-AUC 0.63–0.72).
- **Abstention still does not work at the default targets:** the band collapsed in every fold for both learned models, so `insufficient_evidence` was never used. Stricter targets have not been run on v2.

**Limits:** within-locus errors are 96% minigraph's, so PanAligner cannot be judged on them. There is one accuracy level (75%). Errors from the same regions repeat, so the effective sample is far smaller than 536, and PSEN2 has only 16.

### Strict targets on v2

```bash
python main.py --ml-alignment-strict    # writes outputs/ml_alignment_v2/strict_targets/{window,leave_one_gene_out}/
```

**Method** (`scripts/ml/strict_targets.py`):

- Each training fold is split by group into 60% fit, 20% calibration and 20% thresholds. Thresholds never see the test fold or the calibrator's rows. The committed results confirm no group overlap in any fold.
- This trains on less data than the runs above, so its figures differ slightly from them.
- Six targets were fixed before the run. Equal-coverage comparisons share tied probabilities fairly, and intervals resample whole 20 kb windows.
- A band counts as meaningful only if it opens in every fold, abstains on at least 0.1% of alignments, and its error rate lies between accepted and rejected.

**Gradient boosting by target** (190,772 alignments; requested coverage targets also require 99% rejection precision):

| Grouping | Target | Accepted | Abstentions | Rejected | Error among accepted | Error among abstained | Folds with band | Meaningful |
|---|---|---|---|---|---|---|---|---|
| Window | Error 0.5% | 99.5% | 0 | 1,024 | 0.301% | – | 0/5 | no |
| | Error 0.2% | 99.4% | 234 | 860 | 0.292% | 77.8% | 4/5 | no |
| | Error 0.1% | 99.4% | 316 | 860 | 0.273% | 69.0% | 4/5 | no |
| | Error 0.05% | 79.9% | 37,464 | 860 | 0.272% | 0.86% | 5/5 | yes |
| | Coverage 95% | 97.1% | 4,704 | 857 | 0.147% | 9.89% | 5/5 | yes |
| | Coverage 90% | 96.1% | 6,624 | 857 | 0.147% | 7.07% | 5/5 | yes |
| Leave one gene out | Error 0.5% | 100.0% | 0 | 26 | 0.819% | – | 2/3 | no |
| | Error 0.2% | 99.4% | 1,095 | 26 | 0.275% | 95.2% | 3/3 | yes |
| | Error 0.1% | 99.2% | 1,484 | 26 | 0.205% | 79.2% | 3/3 | yes |
| | Error 0.05% | 53.3% | 89,116 | 26 | 0.286% | 1.43% | 2/3 | no |
| | Coverage 95% | 99.2% | 1,585 | 26 | 0.187% | 76.3% | 3/3 | yes |
| | Coverage 90% | 99.2% | 1,585 | 26 | 0.187% | 76.3% | 3/3 | yes |

**Within-locus and off-target errors, split** (accepted / abstained / rejected):

| Grouping | Target | Within-locus errors (536) | Correct within-locus abstained | Off-target errors (1,052) |
|---|---|---|---|---|
| Window | Error 0.2% | 107 / 85 / 344 | 52 | 447 / 97 / 508 |
| | Coverage 95% | 3 / 191 / 342 | 4,239 | 270 / 274 / 508 |
| | Coverage 90% | 3 / 191 / 342 | 6,156 | 267 / 277 / 508 |
| Leave one gene out | Error 0.2% | 67 / 456 / 13 | 53 | 454 / 586 / 12 |
| | Coverage 95% | 1 / 522 / 13 | 375 | 352 / 688 / 12 |

**Against plain MAPQ at equal coverage.** Difference is model minus MAPQ error in percentage points; accepted errors are within-locus / off-target.

| MAPQ >= | Accepted | MAPQ error | Window: model error, difference (95% interval) | Gene held out: model error, difference (95% interval) | MAPQ accepted errors | Window model accepted errors |
|---|---|---|---|---|---|---|
| 10 | 98.0% | 0.561% | 0.157%, −0.40 (−0.54 to −0.27) | 0.181%, −0.38 (−0.54 to −0.25) | 355 / 693 | 6 / 288 |
| 20 | 96.2% | 0.404% | 0.155%, −0.25 (−0.35 to −0.14) | 0.171%, −0.23 (−0.37 to −0.13) | 231 / 510 | 6 / 279 |
| 30 | 94.6% | 0.295% | 0.158%, −0.14 (−0.21 to −0.07) | 0.167%, −0.13 (−0.25 to −0.05) | 147 / 386 | 6 / 279 |
| 40 | 93.0% | 0.242% | 0.161%, −0.08 (−0.13 to −0.02) | 0.167%, −0.08 (−0.15 to −0.02) | 104 / 325 | 6 / 279 |
| 50 | 91.6% | 0.191% | 0.163%, −0.03 (−0.06 to +0.02) | 0.167%, −0.03 (−0.10 to +0.03) | 70 / 264 | 6 / 279 |
| 60 | 86.7% | 0.128% | 0.172%, +0.04 (−0.01 to +0.14) | 0.167%, +0.04 (−0.03 to +0.09) | 35 / 177 | 6 / 279 |

What these results show:

- **ML does not replace MAPQ on this evidence.** The model is reliably better at MAPQ 10–40, indistinguishable at 50, and nominally worse but statistically unresolved at 60. No operating point shows the model reliably beating MAPQ 60.
- **The two kinds of error behave differently.** The model almost removes within-locus errors at every coverage compared (6 accepted, against 35–355). But at MAPQ 60 it accepts more off-target errors than the threshold (279 against 177).
- **By aligner:** the gain is mostly minigraph's. For PanAligner the model is mostly indistinguishable from MAPQ, with intervals excluding zero only at MAPQ 30 and 40 under window folds.
- **Strict error targets are not achieved.** Error among accepted stays around 0.2–0.3% on held-out data, even when validation thresholds aimed for 0.05–0.2%. Only the 0.5% target under window folds is met.
- **Why there is a floor:** over half of all rows get a predicted p(incorrect) of exactly 0, yet about 0.25% of them are wrong, nearly all off-target (window folds: 264 of 267 errors in the lowest half). These look like genuine alignments, with median identity about 0.83 and MAPQ 50–60. So overall ECE (0.0015–0.0026) looks good while calibration fails at the low end that strict targets depend on.
- **Coverage targets are not honoured either.** Because of tied probabilities, requesting 95% or 90% accepted gives 96–97% under window folds and 99.2% with a gene held out.
- **`insufficient_evidence` is usable only in a limited way.** Coverage targets open it in every fold, and abstained alignments are clearly enriched for errors (7–76%, against 0.83% overall). It catches almost all remaining within-locus errors there, at a cost of up to 6,156 correct within-locus alignments under window folds. But when a gene is held out, it holds errors that should be rejected: only 26 alignments are rejected, and 76% of abstained ones are wrong. The band is not stable enough to report as a dependable three-way classifier.
- **The meaningful criterion is lenient.** Window error-0.05% passes it, though its abstained error (0.86%) is about the base rate, so that band is not useful in practice. The criterion stays as it was fixed before the run, rather than being changed after seeing results.

**Limits:** thresholds rest on only 2–4 windows per fold. Six targets were compared on one dataset, without a further untouched test set. Figures pool different folds' calibrators, and everything remains simulated, locus-scale and 75% accuracy.

## Model B: variant classification experiment (negative result)

**Question:** can simple sequence-derived features recover the ClinVar classification direction (pathogenic or likely pathogenic, against benign or likely benign) of variants in APP, PSEN1 and PSEN2, better than simple baselines?

**Answer: no, not on this evidence.** The model was trained once under a protocol fixed beforehand. This is a valid negative finding about these features and this dataset.

**What it is not:** an Alzheimer's disease classifier, a diagnostic system or an individual-risk predictor. Labels are database classifications of variants, not patient outcomes. Benign does not mean protective. The dataset is not limited to Alzheimer's-specific assertions, and `outputs/model_b_provenance/PROVENANCE_REPORT.md` found disease assertions in these genes are not Alzheimer-specific.

### Why the dataset was replaced

The earlier study used 220 constructed sequences, each a locus reference plus one variant. The repository holds no record of where those labels came from. An audit on 2026-10-30 found 216 of 220 correspond to ClinVar records with the same alleles and no contradicted label, but that consistency does not establish origin (`outputs/model_b_provenance/`). So a new dataset was built directly from ClinVar, with documented rules.

### Dataset provenance and filtering

- **Source:** ClinVar records saved by `outputs/dataset_replacement_feasibility/20261030_audit_01/`, frozen in `outputs/model_b_task1_freeze/20261030_freeze_01/`. Retrieved 2026-10-30; label checks use the saved 2026-10-25 VCF slice. Live ClinVar will change, so only the saved records are stable.
- **Kept:** P/LP or B/LB with at least 1 star (at least 2 stars as a sensitivity tier).
- **Removed:** VUS, conflicting, not provided, below 1 star, any opposing submission, multi-gene records, and variants over 50 bp.
- **Correction:** `VCV001231596` has identical REF and ALT alleles, so it asserts the reference allele and is not a variant. The freeze has 694 rows and was left unchanged. Corrected data lives in `outputs/model_b_dataset_validation/20261030_validation_02/`, which holds the final table and splits.
- **Superseded:** `20261030_validation_01` is not the final dataset. Its APP-held-out missense validation had no benign variant.

| | APP | PSEN1 | PSEN2 | Total |
|---|---:|---:|---:|---:|
| Pathogenic or likely pathogenic | 23 | 132 | 7 | **162** |
| Benign or likely benign | 249 | 132 | 150 | **531** |
| All | 272 | 264 | 157 | **693** |
| At least 2 stars, P / B | 9 / 91 | 51 / 38 | 2 / 55 | 62 / 184 |
| Missense, P / B | 23 / 20 | 125 / 4 | 5 / 9 | 153 / 33 (**186**) |

Class and consequence are closely tied. All 241 synonymous, 213 intronic and 39 UTR variants are benign, while 153 of 162 pathogenic variants are missense.

**Real sample size:** variants were grouped into 40 genomic clusters, and only 13 hold pathogenic variants (12 among missense, with 25 missense clusters in all). That, not 162, is the effective number of independent pathogenic examples. Also, 64% of records rest on one submitter.

### Feature engineering

Built by `scripts/ml/variant_features_v2.py`, in `outputs/model_b_features/20261030_features_01/data/`. The builder reads only identity columns of the label table, and a test enforces this.

| Table | Rows | Features |
|---|---:|---|
| `features_main.csv` | 693 | Variant type, length change, transition flag, strand-neutral substitution, trinucleotide context, CpG context, GC content (±50 and ±200 bp), homopolymer run |
| `features_missense.csv` | 186 | Reference and alternate amino acid, BLOSUM62, change in hydropathy, volume and charge; proline, glycine and cysteine flags |
| `baselines.csv` | 693 | Gene, position, consequence: **baselines and stratifying only, never inputs to the sequence-derived model** |

Amino-acid features are in a separate table because filling them only for missense would let a pooled model read consequence from empty cells. Never used as inputs: classification, review status, stars, submitter, condition, accession, HGVS name, dates, cluster or variant ID.

### Leakage-aware evaluation

Code: `scripts/ml/variant_dataset_validation.py` and `scripts/ml/variant_model_v2.py`. Design: `SPLIT_DESIGN.md` and `LEAKAGE_AUDIT.md` in the validation snapshot.

- **Clusters:** variants in a gene within 1,000 bp are chained together, as are variants on the same protein residue. Clusters never cross train, validation and test.
- **Folds:** 5 folds, stratified on gene and label, seeds 0–4. Each fold is test once, the next validates and three train.
- **Fitting:** class-weighted logistic regression first, gradient boosting second. Settings and threshold are chosen on validation only, and encoders, scaling and imputation are fitted on train only.
- **Primary metric:** balanced accuracy pooled across test folds per seed, then averaged across seeds.
- **Primary comparison:** missense logistic model minus the best baseline, chosen per fold on validation.
- **Intervals:** 95%, 2,000 resamples of whole clusters, paired for differences. Seed 0 gives the headline interval.
- **Checks:** a separate verifier recomputes the primary figures from `predictions.csv`, with 0 failures. The split audit reports no cluster, residue or within-1-kb crossings in region-safe schemes, and it caught deliberately bad inputs.

The time split (`time_2024_not_region_safe`) has 31 clusters crossing roles, so it only checks against later database additions. PSEN1 and PSEN2 are paralogs, and aligned positions between them are not controlled.

### Results

All figures are from `outputs/model_b_training/20261030_training_01/`; means across 5 seeds unless stated.

| Analysis | Balanced accuracy | Sensitivity | Specificity | ROC-AUC |
|---|---:|---:|---:|---:|
| **Missense, 186 variants** | | | | |
| Logistic model | **0.734** | 0.856 | 0.612 | 0.834 |
| Logistic, without GC | 0.622 | 0.633 | 0.612 | 0.642 |
| Logistic, amino-acid features alone | 0.582 | 0.678 | 0.485 | 0.616 |
| Gradient boosting | 0.699 | 0.714 | 0.685 | 0.771 |
| Best baseline, chosen on validation | **0.848** | 0.817 | 0.879 | 0.803 |
| Gene baseline | 0.848 | 0.817 | 0.879 | 0.793 |
| Nearest labelled variant baseline | 0.833 | 0.980 | 0.685 | 0.833 |
| Position vote baseline | 0.808 | 0.889 | 0.727 | 0.859 |
| **All consequences, 693 variants** | | | | |
| Logistic model | **0.702** | 0.717 | 0.686 | 0.713 |
| Logistic, without GC | 0.578 | 0.625 | 0.532 | 0.629 |
| Gradient boosting | 0.711 | 0.728 | 0.693 | 0.774 |
| Consequence-only baseline | **0.969** | 1.000 | 0.938 | 0.965 |
| Gene baseline | 0.783 | 0.815 | 0.751 | 0.740 |

**Primary comparison:** 0.734 against 0.848, a difference of −0.114. On headline seed 0 the difference is −0.177 (95% interval −0.311 to −0.021), wholly below zero. On seeds 1–4, differences run from −0.139 to −0.078, and their intervals include zero. Seed 0 is the model's lowest seed (0.671), so it is the most negative. **Under the protocol there is no demonstrated gain.**

On all consequences, the model is 0.269 below the consequence baseline on seed 0 (−0.354 to −0.175).

**Why simple baselines score so high.** These are properties of the dataset, not biological insight:

- **Gene (chosen in 20 of 25 folds; nearest variant in the other 5):** PSEN1 missense is 125 pathogenic to 4 benign, so "pathogenic if PSEN1" scores well. This is not a biological predictor of an individual variant. It may reflect differences in gene-specific curation and reporting.
- **Consequence:** benign records here are mostly synonymous, intronic and UTR. This shows how the records are composed, not a way of classifying missense variants.
- **Sequence features see neither directly.** Gene and consequence were deliberately left out of the model.

### GC-content investigation

| Evidence | Value |
|---|---|
| Removing GC, missense | 0.734 to 0.622 |
| Removing GC, all variants | 0.702 to 0.578 |
| Headline seed paired difference, missense | 0.023 (−0.108 to 0.176): interval includes zero, and it is the seed with the smallest gap |
| `gc_200` logistic coefficient, missense | Largest in the model: mean −0.72, negative in all 25 folds |
| `gc_200` variance lying between clusters | 88% |
| `gc_200` missense AUC by gene | APP 0.05, PSEN1 0.63, PSEN2 0.78: direction reverses |
| Mean `gc_200` by gene | APP 0.443, PSEN1 0.403, PSEN2 0.573 |

| Missense, within gene (descriptive) | P / B | Model | Without GC | Amino-acid alone |
|---|---:|---:|---:|---:|
| APP | 23 / 20 | 0.66 | 0.46 | 0.49 |
| PSEN1 | 125 / 4 | 0.50 | 0.64 | 0.55 |
| PSEN2 | 5 / 9 | 0.63 | 0.73 | 0.74 |

**Interpretation:** GC content appears to contribute substantial regional or gene-associated signal, mostly through APP, whose pathogenic variants sit in two low-GC clusters. By the protocol's rule, this gain is not reported as sequence signal. This is an inference from the pattern. **It does not prove that sequence holds no biological information.** The features are simple and the data are few. Grouped folds reduce this effect but cannot remove it, since a model can still learn that low-GC regions are pathogenic. The within-gene values rest on as few as 4 benign or 14 variants, and the PSEN2 amino-acid figure is an after-the-fact observation, not a finding.

### Sensitivity analyses

| Scheme | Model | Best baseline | Note |
|---|---:|---:|---|
| At least 2 stars, 246 variants | 0.706 | 0.947 | Consequence alone 0.964 |
| PSEN1 only, all consequences | 0.631 | 0.902 | Only within-gene region test |
| Time, reference 2024-01-07 | 0.742 | 0.955 | Not region-safe; one seed |
| Hold out APP, all consequences | 0.719 | 0.960 | 23 P / 249 B to test |
| Hold out PSEN1, all consequences | 0.519 | 0.985 | Sensitivity 1.000, specificity 0.038 |
| Hold out PSEN2, all consequences | 0.497 | 0.970 | 7 P / 150 B |
| Hold out APP, missense | 0.660 | 0.500 | Baseline cannot exceed 0.5 by construction; GC included |
| Hold out PSEN1, missense | 0.467 | 0.500 | 4 benign to test |
| Hold out PSEN2, missense | 0.524 | 0.500 | 5 P / 9 B |

**Gene-held-out performance does not establish generalisation.** Three genes cannot show that, and test sets here are small and unbalanced.

### Limitations

- **Small dataset,** with only 12 pathogenic missense clusters and wide intervals.
- **Uneven genes:** PSEN1 has 125 of 153 pathogenic missense variants and just 4 benign ones. PSEN2 has 7 pathogenic variants in total.
- **Test folds have been examined.** Any future experiment on these folds is exploratory, not independent validation.
- **No independent test set exists.** Genuine evaluation requires new variants or a later frozen release. The stated minimum is 20 P/LP and 20 B/LB in new clusters.
- **Labels are not independent truth:** they share evidence and laboratories, and conservation-like signals would be partly circular with how variants get classified.
- **Scores are not calibrated probabilities,** because the models are class-weighted.
- **Seeds re-split the same clusters;** they are not independent samples.

### Future work

Each of these would be a new, separately declared analysis:

- Evaluate on a later frozen ClinVar release, with the frozen model and new clusters.
- Add more genes, to widen benign missense coverage and allow a real transfer test.
- Add region-aware features such as protein domain, or conservation, with attention to circularity.
- Compare with external predictors such as REVEL, CADD or AlphaMissense, as comparisons rather than inputs.
- Control paralog similarity between PSEN1 and PSEN2, and test holding out laboratories.

### Reproducing the record

Do not treat a rerun as confirmatory. Both scripts refuse existing folders, so nothing preserved is overwritten. `variant_model_v2.py` also needs pandas, which `requirements-ml.txt` does not list, and with Python 3.13 and pandas 3 the code sets Python-backed strings to avoid a crash seen on this laptop.

```bash
python -B scripts/ml/variant_dataset_validation.py --output <new_folder>
python -B scripts/ml/variant_model_v2.py --output <new_folder>
```

Tests: `tests/test_variant_features_v2.py`, `tests/test_variant_dataset_validation.py`, `tests/test_variant_model_v2.py`. They use synthetic data and test code behaviour, not scientific validity.

### Earlier study (superseded)

`python main.py --ml-variant` still runs the earlier study (`scripts/ml/variant_feasibility.py`, `outputs/ml_variant/REPORT.md`). It used the 220 constructed sequences, whose label source is undocumented. Under 1 kb window-grouped folds, context-only logistic reached 0.550 and the position vote 0.646, falling to 0.453 and 0.591 at 20 kb. It is kept as history and does not describe the current Model B.
