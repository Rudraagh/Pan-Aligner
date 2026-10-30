# Recommendations for Model B

The changes in section 5 were approved and applied to the working tree on 2026-10-30, uncommitted. Slides and existing reports are unchanged.

## 1. What can and cannot be claimed

### Safe claims

- "The labels' original source is undocumented. We audited them on 2026-10-30 against ClinVar: 216 of 220 sequences correspond to a record with the same alleles, and no label is contradicted by the record's current classification."
- "Unhealthy corresponds to variants currently classified pathogenic or likely pathogenic, and healthy to benign or likely benign. That agreement does not prove ClinVar was the source."
- "Model B is a feasibility study of label recovery on constructed sequences. Under window-grouped CV, sequence context does not recover the labels: balanced accuracy is 0.55 at 1 kb, below the position baselines of 0.63 to 0.65."
- "Labels cluster by position and functional consequence, so high scores on random splits would not show generalization."

### Wording to avoid

| Avoid | Reason |
|---|---|
| "Alzheimer's classifier" or "predicts Alzheimer's" | No patients or phenotypes, and no leak-free signal |
| "Healthy sequences" or "disease sequences", without qualification | Labels describe variant classifications, not people |
| "Labels are from ClinVar" | Inferred, not documented |
| "Labels were validated by ClinVar" | ClinVar may be their source, so this may be circular |
| "All variants are Alzheimer's variants" | Assertions also cover other conditions; 84 sequences lack an Alzheimer-only record in the label's direction |
| "Benign means protective", or "benign means healthy" | Different concepts |
| "Swedish mutation" for APP_U_16 | The sequence as encoded is a different allele |
| "71.7% accuracy" as a result | It reflects the random split's neighbour leakage on 46 variants |

## 2. Should the feasibility results stay in the FYP?

**Yes, as a negative result and data audit, kept separate from Model A.** This provenance audit strengthens it. Keep the original labels and original results unchanged.

Disclose these 6 sequences wherever the dataset is described: APP_U_16, PSEN1_U_74, PSEN2_U_1, PSEN1_H_12, PSEN1_U_26, PSEN2_H_1.

## 3. Is a new externally labelled dataset required?

**Yes, for any claim beyond label recovery. It is not required for the FYP as currently scoped.**

Current data cannot be repaired by relabelling. A new study would need:

1. A versioned source and frozen release, with a stored query and checksum.
2. Written inclusion rules set before modelling, covering condition, review status, and treatment of conflicts.
3. A disease question stated at the right level, for example variant pathogenicity for Alzheimer's-named conditions.
4. An independent test set from a different source or a later release, with variants grouped by position.
5. Features that do not encode the classification itself, or an explicit declaration when they do. Consequence and allele frequency often feed ClinVar classifications.
6. For any person-level disease claim: real case-control or cohort genotypes, with ethics and data-access approval. A variant-level dataset cannot support this.

## 4. Defensible future work

- Check archived ClinVar releases from before 2026-05-08 to test whether a specific release and filter reproduces the 220 selections. This could turn inference into evidence, but only if a match is unique.
- Ask whoever created the dataset for the original source.
- Retry Alzforum through an approved route, and compare with ClinGen and LOVD, reporting each one separately.
- Decide how to treat APP_U_16 and PSEN1_U_74. Do this only in a new, versioned dataset, leaving the original files as they are.
- Pathogenicity prediction with locus-held-out validation, framed as variant interpretation research rather than diagnosis.

## 5. Documentation changes (applied)

| File | Current text | Proposed text |
|---|---|---|
| `docs/17_ml_layer.md:8` | "the labels' source is unknown and no model recovers them" | "the labels' source is undocumented, though they agree with current ClinVar classes, and no model recovers them from sequence context" |
| `docs/17_ml_layer.md:205` | "The source, version and criteria of the healthy/unhealthy labels are unknown." | Keep it and add: "A 2026-10-30 audit found 216 of 220 sequences correspond to ClinVar records with the same alleles and no contradicted label, which does not establish origin. See `outputs/model_b_provenance/PROVENANCE_REPORT.md`." |
| `docs/17_ml_layer.md:252` | "none of which are in the repository" | Add that consequence and frequency differ strongly by class, so they would need treating as potential leakage |
| `README.md:264` | "the labels' source is unknown" | "the labels' source is undocumented, though an audit found they agree with current ClinVar pathogenic and benign classes" |
| `README.md:511` | "labels of unknown provenance" | "labels of undocumented provenance" |
| `VIVA_PREPARATION.md:238` | "their source is unknown, and they cluster by position" | "their source is undocumented; they agree with current ClinVar classes but cluster by position and consequence" |
| `VIVA_PREPARATION.md:409` | "The labels' source is unknown..." | "The labels' origin is undocumented, they reflect variant classifications rather than patients, and no model recovers them from sequence context under leak-free CV; it is a feasibility study, not a classifier." |

Slides are not covered here. Any slide on Model B should use the safe claims in section 1.
