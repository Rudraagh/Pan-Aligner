# Model B label provenance and scientific validity

Investigation date and external retrieval date: 2026-10-30. Repository baseline: commit `706dd3d`, branch `wip/gpu-minimizer-stage0`. This report concerns Model B only and says nothing about Model A.

## 1. Executive summary

**Verdict: the original source of the labels is not established. The labels are strongly consistent with current ClinVar classifications, but they do not support a disease classifier.**

- The repository holds no record of where the variants or labels came from: no source file, database version, filter, or generation script.
- 216 of 220 sequences correspond to a ClinVar record with the same alleles (215 by canonical SPDI, 1 by identical coordinates and VCF alleles).
- No label is contradicted by a current ClinVar classification. 215 labels are consistent with the record's overall class, 1 is consistent only if an imprecise deletion is accepted, and 4 cannot be assessed.
- Consistency does not prove origin. ClinVar is the most likely source, or an upstream relative of it, but that is an inference, and the snapshot and selection criteria remain unknown.
- "Unhealthy" corresponds to a pathogenic or likely pathogenic variant classification. "Healthy" corresponds to benign or likely benign. Neither describes a person's health, and neither comes from patient data.
- Disease assertions are not specific to Alzheimer's disease. Only 80 of the 110 unhealthy and 56 of the 110 healthy sequences have an Alzheimer-only ClinVar condition record in the label's direction.
- Model B should stay a label-recovery feasibility study. It must not be called an Alzheimer's classifier.

### Outcome counts

| Category | Healthy | Unhealthy | Total |
|---|---:|---:|---:|
| Exact verified match | 108 | 106 | 214 |
| Likely match but representation uncertain | 0 | 3 | 3 |
| Multiple conflicting annotations | 1 | 0 | 1 |
| Match exists but unrelated to Alzheimer's disease | 0 | 1 | 1 |
| Insufficient evidence | 1 | 0 | 1 |
| No external match | 0 | 0 | 0 |
| **Total** | **110** | **110** | **220** |

"Exact verified match" means the alleles and the record's classification direction were verified. It does not mean the label's origin or any biological disease status was verified.

## 2. Dataset origin

| Finding | Evidence |
|---|---|
| 220 labelled sequences: APP 23+23, PSEN1 78+78, PSEN2 9+9 | `data/raw/*_all_labeled.fa` |
| Each is the locus reference plus one variant record | Pairwise comparison; 220 of 220 reference alleles match |
| References match NCBI exactly | efetch of NC_000021.9:25880550-26171128, NC_000014.9:73136417-73223691, NC_000001.11:226870184-226927726 |
| Headers carry only `APP_H_1`-style names | No rsID, HGVS, source, or date in any header |
| First appearance | Initial commit `2630d39`, 2026-05-09; file modification times 2026-05-08 |
| No generation route in the repository or Git history | No insertion script; no VCF or annotation table; `git log -S` finds no ClinVar, dbSNP, gnomAD, rsID, Alzforum, Benign, or Pathogenic text in the dataset's history outside outputs and tools |

These sequences are constructed, not observed. They carry one variant each on a reference background, so they are not patient haplotypes.

## 3. Original label provenance

**Not established.** Three things are consistent with a ClinVar-like source, and all three are circumstantial:

1. Exact alleles and classification direction agree for almost the whole dataset, including rare and non-SNV records.
2. PSEN2_U_1, before left-alignment, starts and ends exactly at the inner bounds of VCV000533784, `(?_226881888)_(226895599_?)del`. Its true breakpoints are unknown, yet the dataset treats those bounds as exact.
3. PSEN1_U_59 equals VCV000038296 at the coordinate level, a record with no canonical SPDI.

Two sequences suggest a possible encoding error. This is a hypothesis from two examples and is not established:

| Sequence | ClinVar record at the same site | Dataset allele | Observation |
|---|---|---|---|
| APP_U_16 | VCV000018093, c.2010_2011inv, the plus-strand change TC to GA | Plus-strand TC to CT | CT is TC reversed but not complemented |
| PSEN1_U_74 | VCV000098034, c.466_467insTTATAT | c.466_467insTATATT | TATATT is TTATAT reversed |

As encoded, these two sequences are different alleles: VEP gives p.Met671Val and p.Lys155_Tyr156insLeuTyr, not the recorded p.Lys670_Met671delinsAsnLeu and p.Lys155_Tyr156insPheIle. Their "unhealthy" label cannot be attached to the sequence actually in the FASTA.

Timing also limits the inference. 39 matched records have a submission newer than the dataset date (2026-05-08), though no matched record was created after it. Today's ClinVar is therefore not the snapshot the labels could have used.

Selection criteria are unknown. Exact label-consistent matches cover only part of what ClinVar held on 2026-10-30:

| Gene | Unhealthy matched / P/LP records in locus | Healthy matched / B/LB records in locus |
|---|---:|---:|
| APP | 22 / 28 | 23 / 283 |
| PSEN1 | 77 / 152 | 77 / 133 |
| PSEN2 | 8 / 12 | 8 / 158 |

## 4. Variant normalization

Script: `scripts/01_normalize_variants.py`.

- Each sequence was compared with its locus reference. Header coordinates are 1-based inclusive.
- Indels were left-aligned using the locus reference and written as VCF-anchored records, HGVS g. strings, and 0-based SPDI. Nine variants were shifted.
- Types: 205 SNV, 3 adjacent MNV, 4 insertions, 8 deletions.
- All records use GRCh38 (headers state GRCh38.p14). No liftover was used.

## 5. External matching

| Source | Access | Outcome |
|---|---|---|
| ClinVar | NCBI E-utilities esearch, esummary, efetch VCV XML | 1,655 locus records; full records for all matches |
| dbSNP | NCBI Variation Services, SPDI to rsID | 225 lookups, 214 with rsIDs, 2 failed |
| Ensembl VEP | REST region POST | Consequences for all but the 2 deletions over 1 kb, which were not queried |
| gnomAD r4 | GraphQL by region | Full-locus requests for APP and PSEN1 returned HTTP 500; 40 kb chunks then succeeded with no failed chunks |
| Alzforum | Page request per gene | **Failed: HTTP 403 on all three. No Alzforum data was used** |
| OMIM, ClinGen, LOVD, HGMD, PubMed | Not queried | OMIM appears only as a ClinVar submitter |

Failed dbSNP lookups: PSEN1_U_59 (4,556 bp deletion, HTTP 503) and PSEN2_U_1 (13,712 bp deletion, HTTP 414). These are not evidence of absence.

Matching rules:

- ClinVar canonical SPDI was trimmed and left-aligned to the same form as the dataset variant, then compared exactly.
- Eight locus records have no canonical SPDI. They were checked by coordinates, and the related records were fetched by site, without reference to label.
- Labels were never used to select a record. No majority vote across databases was taken. Ensembl clinical significance derives from ClinVar and was not counted as independent.

## 6. Per-gene findings

| Gene | Exact verified | Likely match | Conflicting | Unrelated to AD | Insufficient |
|---|---:|---:|---:|---:|---:|
| APP (46) | 45 | 1 (APP_U_16) | 0 | 0 | 0 |
| PSEN1 (156) | 153 | 1 (PSEN1_U_74) | 1 (PSEN1_H_12) | 1 (PSEN1_U_26) | 0 |
| PSEN2 (18) | 16 | 1 (PSEN2_U_1) | 0 | 0 | 1 (PSEN2_H_1) |

Sequences that were not exact verified matches:

| Sequence | Label | Finding |
|---|---|---|
| APP_U_16 | Unhealthy | Different allele from Swedish-site record VCV000018093. Only the c.2010G>A component has an exact record, VCV002856224, which is Likely benign. Not assessable |
| PSEN1_U_74 | Unhealthy | Different inserted sequence from VCV000098034. Not assessable |
| PSEN2_U_1 | Unhealthy | VCV000533784 (Pathogenic, Alzheimer disease 4, 1 star) has unknown outer breakpoints. Consistent only if those bounds are accepted |
| PSEN1_H_12 | Healthy | VCV001153577 is Conflicting overall. Its bundled condition record says Likely benign, and one submission is uncertain. Not assessable |
| PSEN1_U_26 | Unhealthy | VCV000018149 is Pathogenic overall, with the pathogenic assertion for Pick disease. Its Alzheimer disease 3 record says "not provided". Consistent with the overall class, but not Alzheimer-supporting |
| PSEN2_H_1 | Healthy | rs186308126, gnomAD allele frequency 0.0055, no ClinVar record at this position. No disease assertion, so not assessable |

Submission-level disagreement exists at APP_H_7 (VCV000585430, aggregate Benign/Likely benign, one uncertain submission) and PSEN1_H_12.

## 7. Alzheimer's relevance

Variant pathogenicity is not disease-specific association. These results are about assertions in condition records.

| Alzheimer evidence in the label's direction | Healthy | Unhealthy |
|---|---:|---:|
| Alzheimer-only condition record | 56 | 80 |
| Alzheimer named alongside other conditions in the record | 53 | 26 |
| None | 0 | 1 |
| No exact record | 1 | 3 |

- By condition name, in the unhealthy group: Alzheimer's disease 97, familial 6, early-onset 3, uncertain 1. In the healthy group: Alzheimer's disease 109. These are name matches only; the names "Alzheimer disease 3" and "type 1" were not reclassified.
- Conditions are not exclusive. Unhealthy records also carry Pick disease (59), Frontotemporal dementia (57), Acne inversa, familial, 3 (56), Cerebral amyloid angiopathy, APP-related (7), and Dilated cardiomyopathy 1U (3). Healthy records also carry Dilated cardiomyopathy 1U (25) and 1V (5).
- Several APP sequences are recorded with distinct amyloid phenotypes, for example APP_U_10 (Dutch type, cerebral amyloid angiopathy).
- "Healthy" is not a neutral group. APP_H_2 (p.Ala673Thr) has "Alzheimer disease, protection against = protective" and "Alzheimer disease = Benign". It is not simply a variant without effect.
- A benign assertion means a laboratory judged the variant not to cause the named condition. It does not show the carrier is healthy.
- Review strength is limited. No record has 3 or 4 stars. Healthy: 61 at 1 star, 48 at 2 stars. Unhealthy: 17 at 0 stars, 41 at 1 star, 49 at 2 stars. Last-evaluated dates range from 1996-12-01 to 2026-09-25.

## 8. Label agreement

| Current ClinVar class | Healthy | Unhealthy |
|---|---:|---:|
| Benign or likely benign | 108 | 0 |
| Pathogenic or likely pathogenic | 0 | 107 |
| Conflicting | 1 | 0 |
| No exact record | 1 | 3 |

ClinVar cannot serve as independent validation if it is, as seems likely, a source or near-source of the labels. High agreement would then be partly circular.

The per-variant table with the required columns is in `provenance_comparison_table.md`.

## 9. Validity questions

1. **Are the labels independently verifiable?** Partly. The classification direction can be checked against ClinVar today. The labels' origin and selection cannot be checked, and ClinVar may not be independent of them.
2. **Are they biologically meaningful?** Only as variant-level pathogenicity assertions. They are not phenotypes, carrier states, or patient outcomes, and the sequences are constructed.
3. **Do they correspond to Alzheimer's specifically?** No. 136 of 220 have an Alzheimer-only record in the label's direction. The rest are bundled with other conditions, lack Alzheimer support, or have no exact record. Many records also name other conditions.
4. **Is the dataset large enough?** Not for a disease claim. There are 220 variants, with PSEN2 contributing 18, and only 110 per class. Performance per gene is especially unstable.
5. **Can positional confounding explain performance?** Largely. Under window-grouped CV at 1 kb, balanced accuracy was 0.646 for position vote and 0.628 for nearest variant, against 0.550 for context only. All fell as windows grew. Nearest-variant label agrees for 87.0% of APP, 69.9% of PSEN1, and 88.9% of PSEN2 variants. APP unhealthy positions span only locus coordinates 11,228 to 17,076, against 153 to 101,774 for healthy.
6. **Does variant type differ by class?** By broad type, not significantly: SNVs are 106 healthy and 99 unhealthy (Fisher p = 0.106). By functional consequence, it differs enormously: 108 of 108 queried unhealthy are protein-altering or splice-site, against 19 of 110 healthy (Fisher p = 8.5e-43). Allele frequency also differs: 23 healthy and 0 unhealthy reach 1% in gnomAD r4, and 19 unhealthy against 103 healthy are present at all. Any model given these annotations would recover how variants were classified, not disease biology. No such model was fitted.
7. **Is there split leakage through neighbours?** Yes, in the seed-42 random split. 29 of 46 test variants are within 10 bp of a training variant, 39 within 100 bp, and 43 within 1 kb. The nearest training variant has the same label for 33 of 46, the 0.717 reported previously. Fifteen positions carry more than one variant, none with both labels. Random-split scores should not be reported as generalization.
8. **Can any result be interpreted as disease prediction?** No. There are no individuals, phenotypes, or independent test cohort, and there is no leak-free signal from sequence context.
9. **What more is needed?** See section 3 of `RECOMMENDATIONS.md`. In short: versioned sources, explicit inclusion rules, evidence filters, and real case-control or independently curated disease data.

## 10. Limitations

- All external results reflect 2026-10-30, not the label creation date. ClinVar snapshots from before May 2026 were not checked.
- Alzforum failed and other curated sources were not queried, so non-ClinVar origin is not excluded.
- Long deletions have no dbSNP or VEP results. gnomAD was checked by exact variant ID only, and MNVs were not looked up as a whole.
- Absence from gnomAD does not show a variant is absent from populations.
- Some ClinVar condition names were recorded as supplied. For example, PSEN1 record VCV000038296 lists "Alzheimer disease 4". Those names were not independently audited.
- The Alzheimer classes come from condition-name text, and the "early-onset" and "familial" classes undercount because of that.
- The reversal explanation for APP_U_16 and PSEN1_U_74 is a hypothesis.
- General statements about what classifications mean are interpretation, not retrieved data.
- No model was trained or retrained. Confound checks are descriptive, and the performance values come from the existing `outputs/ml_variant/REPORT.md`.

## 11. Final scientific verdict

| Question | Answer |
|---|---|
| Original label provenance established? | **No** |
| Labels consistent with an authoritative source today? | **Yes**, strongly, at the level of variant classification |
| Labels defensible as "ClinVar-consistent pathogenic versus benign variants"? | **Yes**, with the 6 exceptions disclosed and the source stated as inferred; no label is changed |
| Labels defensible as "healthy versus Alzheimer's" or as a disease ground truth? | **No** |
| Model B can be called an Alzheimer's classifier? | **No** |

## 12. Sources and retrieval dates

All retrieved 2026-10-30. Exact URLs and times are in `raw/` and `source_inventory.md`.

- NCBI ClinVar and Nucleotide, E-utilities: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
- ClinVar record pages: https://www.ncbi.nlm.nih.gov/clinvar/variation/ followed by the Variation ID, listed per variant in `variant_provenance.csv`
- NCBI Variation Services: https://api.ncbi.nlm.nih.gov/variation/v0/
- Ensembl VEP REST: https://rest.ensembl.org/vep/homo_sapiens/region
- gnomAD API, dataset gnomad_r4: https://gnomad.broadinstitute.org/api
- Alzforum, attempted and refused with HTTP 403: https://www.alzforum.org/mutations/app, `psen-1`, `psen-2`
