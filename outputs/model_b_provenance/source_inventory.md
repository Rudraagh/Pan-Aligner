# Source inventory

## Repository inputs (read only)

| File | SHA-256 | Use |
|---|---|---|
| `data/raw/app_all_labeled.fa` | `1ba1e982c39694c164d8a5d33f19c5a9340ec313c406c5019f626307c362d8c8` | Labelled sequences |
| `data/raw/psen1_all_labeled.fa` | `735b9571c408037adbe65bd3f4853c7fe6426b2fb3d733f52a6caf76e3ebe53e` | Labelled sequences |
| `data/raw/psen2_all_labeled.fa` | `bee7507dd3644e99c9cdf94c72d0ad03a85e208616696a0cc601a523db67107b` | Labelled sequences |
| `data/raw/app_reference.fa` | `8af944d3663bcb47d1e5399ea30260efddd147eacc77681b3e583424987cd2cb` | Locus reference |
| `data/raw/psen1_reference.fa` | `cb91993f71437c5e35fe0b45931b3bc774b111fa2cf7a7577d4e825cccf12d99` | Locus reference |
| `data/raw/psen2_reference.fa` | `517347854c7a0f8d948a168f6726b10ae40059ce4acb05651e739e22c13ccf70` | Locus reference |
| `app_combined.fasta` | `e282299b119b541bf19b97dca4a8136047622d5ef40d9e107a8ad8844ebb5f59` | Origin comparison |
| `psen1_combined.fasta` | `711253422dc33eb5d152787a0931ddafcdbf17b49784fd34bafb8075aefe37f0` | Origin comparison |
| `psen2_combined.fasta` | `0160c69aa5ad099ea29949db29c93aaa5539a64d89df8e54764aab472644501d` | Origin comparison |

Also read: `scripts/preprocess.py`, `scripts/common.py`, `scripts/ml/variant_feasibility.py`, `outputs/ml_variant/REPORT.md`, `outputs/ml_variant/feasibility_results.json`, `outputs/ml_variant/variant_table.csv`, `data/metadata/`, and Git history from initial commit `2630d39`.

**Not found anywhere in the repository or its history:** a variant source table, VCF, manifest of labels, database version, selection rule, or sequence-generation script.

## Investigation scripts

| Script | Purpose |
|---|---|
| `scripts/01_normalize_variants.py` | Derive and normalize one record per sequence; verify references with NCBI |
| `scripts/02_query_external.py` | Query ClinVar, dbSNP, Ensembl VEP, gnomAD; log failures |
| `scripts/03_match_and_fetch.py` | Exact-allele matching; fetch VCV records; retry gnomAD in chunks |
| `scripts/03b_fetch_related_records.py` | Fetch site-related records; attempt Alzforum |
| `scripts/04_build_tables.py` | Build `variant_provenance.csv`, `unmatched_variants.csv`, `raw/analysis_summary.json` |
| `scripts/05_confounds.py` | Descriptive confound checks; fits no model |
| `scripts/06_build_appendices.py` | Build comparison and submission tables; print checksums |

## External databases

All retrieved 2026-10-30.

| Source | Endpoint | Cached response | Status |
|---|---|---|---|
| NCBI Nucleotide | `efetch.fcgi?db=nuccore` | URLs in `raw/normalized_variants.json` | Success; 3 of 3 references match |
| ClinVar search and summaries | `esearch.fcgi`, `esummary.fcgi`, `db=clinvar` | `raw/clinvar_locus_summaries.json` | Success: APP 639, PSEN1 637, PSEN2 379 |
| ClinVar full records | `efetch.fcgi`, `rettype=vcv` | `raw/clinvar_vcv_batch_000.xml` to `005`, `raw/clinvar_vcv_related.xml` | Success |
| dbSNP via NCBI Variation Services | `/variation/v0/spdi/{spdi}/rsids` | `raw/dbsnp_spdi_rsids.json` | 225 lookups; 2 failed |
| Ensembl VEP | `/vep/homo_sapiens/region` | `raw/ensembl_vep.json` | Success; 2 deletions over 1 kb not queried |
| gnomAD r4 | GraphQL `region` | `raw/gnomad_regions.json` | APP and PSEN1 full-locus HTTP 500, then success in 40 kb chunks; PSEN2 success |
| Alzforum | `/mutations/app`, `psen-1`, `psen-2` | `raw/related_fetch_log.json` | **Failed: HTTP 403 on all three; no data used** |

Failed requests are recorded in `raw/request_failures.json` and `raw/related_fetch_log.json`.

## Sources not queried

OMIM (directly), ClinGen, LOVD, HGMD, PubMed, and historical ClinVar archives. No literature was used as evidence. Interpretive statements in the report are marked as such.

## Output files

| File | Content |
|---|---|
| `PROVENANCE_REPORT.md` | Findings and verdict |
| `RECOMMENDATIONS.md` | Claims, wording, future work, proposed documentation changes |
| `variant_provenance.csv` | One row per variant, with identifiers, annotations, and record details |
| `unmatched_variants.csv` | Variants not confidently matched, with reasons |
| `provenance_comparison_table.md` | Required comparison table for all 220 variants |
| `clinvar_submissions.csv` | One row per ClinVar submission: submitter, date, classification, evaluation date |
