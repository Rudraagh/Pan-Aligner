# Model B follow-up protocol

Version 1, 2026-10-31. Written at commit `00f1478`, before any step below has run. **Status: draft until committed. It takes effect at that commit, which must come before any fetch of step 1 or 2 scores.**

Background: `../20261031_investigation_01/`. The original negative result (`outputs/model_b_training/20261030_training_01/`) is unchanged and remains the confirmatory result.

## Questions

| Step | Question | Status of any result |
|---|---|---|
| 1 | Do existing predictors separate these missense variants beyond gene and position? | **Exploratory** |
| 2 | Does one conservation and topology feature set add anything to a simple model? | **Exploratory** |
| 3 | Does a frozen model agree with an independent functional readout on unseen variants? | **Independent variants, different target** |

No step can show that Model B classifies disease, risk or Alzheimer's.

## General rules

1. **One run per step.** Dry runs use synthetic or permuted labels only. A failed run that produced no metrics may be repeated after a recorded fix.
2. **Log every run** in `RUN_LOG.md` in this folder, with date, commit, command, input checksums, package versions, and outcome. Include failed and abandoned runs.
3. **Changes:** any later change goes in an `AMENDMENTS.md` file here, dated, and stating whether results had been seen. A change made after seeing results makes that analysis post hoc, and it's labelled so.
4. **New folders only.** No preserved file is edited. Original splits are used read-only.
5. **Never inputs:** gnomAD frequency or absence, gene, position, consequence, REVEL, AlphaMissense, CADD, SIFT, PolyPhen, ClinVar fields.
6. **Environment:** the existing Python 3.13.7, numpy 2.3.3, pandas 3.0.2, scikit-learn 1.8.0. No new packages.
7. **Approval:** implementing and running each step needs separate approval. Step 3's data retrieval is manual.

### Inputs (SHA-256)

| File | SHA-256 |
|---|---|
| `outputs/model_b_dataset_validation/20261030_validation_02/CORRECTED_VARIANTS.csv` | `84e78413fd0ce1e85c668c1819418a6ac5bafe2bf28e0bf6aa6b3dd86f5aba92` |
| `outputs/model_b_dataset_validation/20261030_validation_02/SPLIT_MANIFEST.csv` | `df0342492c8c0990a5def2881c6c6ca08f76beb443131d7565d41d212b2b6800` |
| `outputs/model_b_features/20261030_features_01/data/features_main.csv` | `0faecfe6f0627972a952c2df4b18f40c897596e616e3e88eb510d267b3ff4345` |
| `outputs/model_b_features/20261030_features_01/data/features_missense.csv` | `e1a6de2c8e982fe6364dfba3c8e6b92e42f448bd27cace65eacf4e51fe3e51ec` |
| `outputs/model_b_features/20261030_features_01/data/baselines.csv` | `a07ed2516d881315fe5a2727ab26ae5c50e3200cdb5cf80a70ec9492fa0251ba` |
| `outputs/model_b_features/20261030_features_01/fetched/raw/ensembl_vep.json.gz` | `d1912937de248496dfe101151657cf4af1bdbb2ca7de537445b669d10a01c4e4` |
| `outputs/model_b_training/20261030_training_01/predictions.csv` | `82b097995b5f4d005bdf6c715f8eb24a2db678c4dc12f5624a82da08243764a0` |

### Shared method

Scheme `missense_primary`: seeds 0–4, 5 folds. Thresholds use `choose_threshold` from `scripts/ml/variant_model_v2.py` on the validation role only, and apply unchanged to test. Balanced accuracy is pooled per seed, then averaged. Intervals: 2,000 resamples of whole clusters, seed 20261031, paired, with seed 0 as headline.

## Step 1: predictor benchmark

- **Primary comparator: AlphaMissense.** Secondary: REVEL. Descriptive: SIFT (1 − score), PolyPhen, and CADD if available.
- **Fetch:** Ensembl VEP, canonical transcript, saving raw responses and a log. Retry failures in small batches. No score is compared with labels until coverage is recorded.
- **Coverage rule:** variants without a score are excluded from that predictor's analysis, and the baseline is recalculated on the same variants. Below 95% of missense, the predictor is descriptive only.
- **No fitting** other than thresholds.
- **Primary comparison:** AlphaMissense minus best baseline, using the preserved per-fold best-baseline predictions.
- **Also reported:** the same against the position vote baseline; ROC-AUC overall and within APP; the original model (0.734) on the same variants.

**Gate.** Step 2 runs only if AlphaMissense's mean balanced accuracy exceeds the best baseline's on the same variants. This decides what to do next. It isn't a claim of gain, because the dataset can't resolve that.

## Step 2: one feature set

Runs only if step 1's gate passes.

- **Added features, exactly two:**
    - `conservation`: VEP's Conservation output for the variant.
    - `in_transmembrane`: 1 if `protein_position` lies in a UniProt Transmembrane feature of P05067, P49768 or P49810.
- **If PSEN2's UniProt record stays unavailable** after retries, `in_transmembrane` is dropped and that is recorded as an amendment before modelling. Missing conservation is median-imputed from train.
- **Before modelling:** label-blind check of each feature's share of variance between clusters and genes. Reported whatever it shows. It doesn't remove a feature.
- **Reduced base inputs:** `variant_type`, `length_change`, `is_transition`, `substitution`, `in_cpg`, `gc_50`, `gc_200`, `homopolymer_run`, `blosum62`, `delta_hydropathy`, `delta_volume`, `delta_charge`, `to_proline`, `from_glycine`, `from_cysteine`, `to_cysteine`. This excludes `aa_ref`, `aa_alt` and `trinucleotide`.
- **Model:** class-weighted logistic regression, C in 0.01, 0.1, 1, 10, chosen on validation. No other model.
- **Arms, all four reported:** base; base plus features; each with and without GC.
- **Interpreted comparison:** no-GC base plus features, against no-GC base. Also shown against best baseline and position vote.

## Step 3: frozen model on independent data

Runs whatever step 1 shows, subject to the data rules.

**Frozen models, fixed now from `fold_choices.csv`:** two logistic models trained once on all 186 missense variants, using original inputs and class weighting.

| Model | C (most often chosen on validation) | Threshold (median across 25 folds) |
|---|---:|---:|
| Missense | 0.01 (15 of 25 folds) | 0.5245 |
| Missense without GC | 0.01 (10 of 25 folds) | 0.5822 |

Save both with checksums **before external data is opened.** Baselines are fitted on the same 186 variants: gene rate and position vote.

**Choosing the dataset, using papers only and before any model score is seen.** Candidates in order: the 2025 study of 114 VUS; the 2020 cell-based validation; Petit et al. 2022; the 2025 PSEN2 and APP study. Use the first that passes:

| Requirement | Rule |
|---|---|
| Per-variant values | Public table with a licence permitting this use |
| Readout | Measured Aβ-related experiment, not a computational score |
| Class call | **Authors' own functional classification.** If they give none, define a cut from their controls and record it as an amendment before scoring |
| Unseen variants | Missense, and not among the 186. Those in the 693 by other routes are excluded and counted |
| Minimum | At least 20 variants in each class to report as a test. 10–19 in either class: descriptive only. Fewer: stop, or try the next candidate |

Datasets may be pooled only if declared in an amendment before scoring.

**Analysis:** one scoring pass. Primary metric ROC-AUC for the frozen no-GC model, then balanced accuracy at the fixed threshold. Reported beside the GC model, both baselines and AlphaMissense. Split by gene, and by whether a variant's cluster occurs in training. Intervals resample clusters, defined as in the validation snapshot.

**Add-on, descriptive:** in APP's Aβ region only, the rank correlation between the no-GC model's score and nucleation in MaveDB `urn:mavedb:00000113-a-3`.

## How results may be described

| Outcome | Allowed statement |
|---|---|
| Steps 1–2 show no gain | "Added evidence did not clear simple baselines in exploratory analysis." |
| Steps 1–2 show a gain | "Exploratory analysis suggests possible value. This dataset can't confirm it, and predictors share evidence with the labels." |
| Step 3 agrees | "A frozen model agreed with an independent functional readout on unseen variants, for this mechanism." |
| Step 3 doesn't agree | "It did not generalise to an independent functional readout." |
| Data rules fail | "No suitable independent dataset was available." |
| **Any** | Never "improved classifier", "validated", "predicts pathogenicity" or "predicts Alzheimer's". |

## Known limits, stated in advance

- All 693 variants have been in read test folds. The primary comparison has little power: a true balanced accuracy of 0.95 clears zero in about 41% of simulations.
- AlphaMissense and other predictors are not independent of how classifications are made.
- Functional calls may have fed into ClinVar (ACMG PS3), and one assay reads one mechanism. Normal in an assay does not mean benign.
- Frozen C and thresholds were derived from folds whose test results have been read, though they were chosen on validation.
- The candidate studies are known only from search summaries so far.
