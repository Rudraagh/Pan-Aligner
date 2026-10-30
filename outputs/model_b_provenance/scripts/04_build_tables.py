"""Combine normalized variants with external records into variant_provenance.csv and unmatched_variants.csv.

Original labels are copied unchanged and never used to choose an annotation. Writes
raw/analysis_summary.json for the report. Ensembl's clinical significance is ClinVar-derived, so it
is recorded but never counted as independent support, and no votes are taken across databases.
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
RAW = OUT / "raw"
STARS = {"practice guideline": 4, "reviewed by expert panel": 3, "criteria provided, multiple submitters, no conflicts": 2,
         "criteria provided, conflicting classifications": 1, "criteria provided, single submitter": 1,
         "no assertion criteria provided": 0, "no classification provided": 0, "no classifications from unflagged records": 0}
# Records found by coordinates at sites the SPDI match could not resolve (03b_fetch_related_records.py). Chosen by site, not label.
RELATED = {
    "PSEN1_U_59": ("38296", "exact", "Deletion start, stop and VCF alleles equal the dataset's after left-alignment; the record has no canonical SPDI"),
    "PSEN2_U_1": ("533784", "imprecise", "Dataset deletion equals the record's inner bounds, but the record's outer breakpoints are unknown, so the exact allele is not defined"),
    "PSEN1_U_74": ("98034", "different_allele", "Same insertion site and length, different inserted sequence: dataset c.466_467insTATATT (VEP p.Lys155_Tyr156insLeuTyr), record c.466_467insTTATAT (p.Lys155_Tyr156insPheIle)"),
    "APP_U_16": ("18093", "different_allele", "Same two bases, different alleles: dataset c.2010_2011delinsAG (VEP p.Met671Val), record c.2010_2011inv (p.Lys670_Met671delinsAsnLeu); only c.2010G>A alone has an exact record"),
}


def condition_class(name: str) -> str:
    low = name.lower()
    if "alzheimer" in low:
        return "alzheimer"
    if low in ("not provided", "not specified", "see cases", ""):
        return "unspecified"
    if re.search(r"inborn genetic|hereditary disease|-related (condition|disorder)|related disorder", low):
        return "gene_or_generic"
    if re.search(r"dementia|frontotemporal|pick|amyloid angiopathy|cerebral|neuro|parkinson|ataxia|spastic|myoclon|lewy|prion", low):
        return "other_neurological"
    return "other_condition"


def text(element, path, default=""):
    found = element.find(path) if element is not None else None
    return (found.text or "").strip() if found is not None and found.text else default


def parse_vcv() -> dict:
    records = {}
    for path in sorted(RAW.glob("clinvar_vcv_batch_*.xml")) + [RAW / "clinvar_vcv_related.xml"]:
        root = ET.fromstring(path.read_text(encoding="utf-8"))
        for archive in root.findall("VariationArchive"):
            classified = archive.find("ClassifiedRecord")
            germline = classified.find("Classifications/GermlineClassification") if classified is not None else None
            submissions = []
            for assertion in (classified.findall("ClinicalAssertionList/ClinicalAssertion") if classified is not None else []):
                accession = assertion.find("ClinVarAccession")
                classification = assertion.find("Classification")
                traits = [name.text for name in assertion.findall("TraitSet/Trait/Name/ElementValue") if name.text]
                submissions.append({
                    "scv": accession.get("Accession") if accession is not None else "",
                    "submitter": accession.get("SubmitterName") if accession is not None else "",
                    "submitted": assertion.get("SubmissionDate") or (assertion.find("ClinVarSubmissionID").get("submitterDate") if assertion.find("ClinVarSubmissionID") is not None else ""),
                    "classification": text(classification, "GermlineClassification"),
                    "last_evaluated": classification.get("DateLastEvaluated", "") if classification is not None else "",
                    "review": text(classification, "ReviewStatus"),
                    "assertion_traits": traits,
                })
            conditions = []
            for rcv in (classified.findall("RCVList/RCVAccession") if classified is not None else []):
                names = [t.text for t in rcv.findall("ClassifiedConditionList/ClassifiedCondition") if t.text]
                rcv_class = rcv.find("RCVClassifications/GermlineClassification")
                conditions.append({"rcv": rcv.get("Accession"), "conditions": names, "classification": text(rcv_class, "Description"),
                                   "review": text(rcv_class, "ReviewStatus")})
            location = next((l.attrib for l in (classified.findall("SimpleAllele/Location/SequenceLocation") if classified is not None else []) if l.get("Assembly") == "GRCh38"), {})
            records[archive.get("VariationID")] = {
                "location": dict(location),
                "vcv": archive.get("Accession"), "version": archive.get("Version"), "name": archive.get("VariationName"), "updated": archive.get("DateLastUpdated"),
                "classification": text(germline, "Description"), "review": text(germline, "ReviewStatus"),
                "last_evaluated": germline.get("DateLastEvaluated", "") if germline is not None else "",
                "created": germline.get("DateCreated", "") if germline is not None else "", "most_recent_submission": germline.get("MostRecentSubmission", "") if germline is not None else "",
                "submissions": submissions, "rcvs": conditions,
            }
    return records


def normal_class(description: str) -> str:
    low = description.lower()
    if "conflicting" in low:
        return "conflicting"
    pathogenic = "pathogenic" in low
    benign = "benign" in low
    if pathogenic and benign:
        return "conflicting"
    if pathogenic:
        return "pathogenic_or_likely"
    if benign:
        return "benign_or_likely"
    if "uncertain" in low:
        return "uncertain"
    return "other" if low else "none"


def condition_summary(record: dict) -> tuple[list[str], str, list[str]]:
    """Conditions tied to the record, the strongest relevance class, and those with P/LP assertions."""
    names, pathogenic_names = [], []
    for rcv in record["rcvs"]:
        names.extend(rcv["conditions"])
        if normal_class(rcv["classification"]) == "pathogenic_or_likely":
            pathogenic_names.extend(rcv["conditions"])
    classes = {condition_class(n) for n in names}
    for wanted in ("alzheimer", "other_neurological", "other_condition", "gene_or_generic", "unspecified"):
        if wanted in classes:
            return sorted(set(names)), wanted, sorted(set(pathogenic_names))
    return sorted(set(names)), "none", sorted(set(pathogenic_names))


def alzheimer_evidence(record: dict, label: str) -> dict:
    """Which AD-named RCVs support the label direction. AD-only RCVs are stronger than bundled ones."""
    want = "pathogenic_or_likely" if label == "UNHEALTHY" else "benign_or_likely"
    named = [r for r in record["rcvs"] if any("alzheimer" in c.lower() for c in r["conditions"])]
    alone = [r for r in named if all("alzheimer" in c.lower() for c in r["conditions"]) and normal_class(r["classification"]) == want]
    bundled = [r for r in named if normal_class(r["classification"]) == want and r not in alone]
    tier = "AD-only RCV" if alone else ("AD named with other conditions in the RCV" if bundled else "no AD RCV in the label direction")
    support = " ".join(c.lower() for r in alone + bundled for c in r["conditions"] if "alzheimer" in c.lower())
    if not support:
        kind = "uncertain or conflicting" if named else "unrelated"
    elif "familial" in support:
        kind = "familial Alzheimer's disease"
    elif "early" in support:
        kind = "early-onset Alzheimer's disease"
    else:
        kind = "Alzheimer's disease"
    return {"tier": tier, "kind": kind, "rcvs": "; ".join(f"{','.join(r['conditions'])}={r['classification']}" for r in named)[:400]}


def categorise(row: dict, record: dict | None, partial: bool, has_rsid: bool) -> tuple[str, str]:
    """(provenance category, label agreement). Agreement never alters the original label."""
    label = row["original_label"]
    if record is None:
        if row["sequence_id"] == "PSEN2_U_1":
            return "Likely match but representation uncertain", "consistent only if breakpoints are accepted"
        if partial or row["sequence_id"] in RELATED:
            return "Likely match but representation uncertain", "not assessable"
        return ("Insufficient evidence" if has_rsid else "No external match"), "not assessable"
    cls = normal_class(record["classification"])
    _, relevance, pathogenic_names = condition_summary(record)
    if cls == "conflicting":
        return "Multiple conflicting annotations", "not assessable"
    if cls in ("uncertain", "other", "none"):
        return "Insufficient evidence", "not assessable"
    if cls == "pathogenic_or_likely":
        agreement = "consistent" if label == "UNHEALTHY" else "contradicted"
        if not any(condition_class(n) == "alzheimer" for n in pathogenic_names):
            return "Match exists but unrelated to Alzheimer's disease", agreement
        return "Exact verified match", agreement
    agreement = "consistent" if label == "HEALTHY" else "contradicted"
    return "Exact verified match", agreement


def main() -> None:
    data = json.loads((RAW / "normalized_variants.json").read_text(encoding="utf-8"))
    matches = json.loads((RAW / "clinvar_matches.json").read_text(encoding="utf-8"))
    vcv = parse_vcv()
    summaries = json.loads((RAW / "clinvar_locus_summaries.json").read_text(encoding="utf-8"))
    dbsnp = json.loads((RAW / "dbsnp_spdi_rsids.json").read_text(encoding="utf-8"))
    vep = json.loads((RAW / "ensembl_vep.json").read_text(encoding="utf-8"))["results"]
    gnomad = json.loads((RAW / "gnomad_regions.json").read_text(encoding="utf-8"))
    rows, unmatched = [], []
    for row in data["variants"]:
        match = matches[row["sequence_id"]]
        components = row["components"]
        spdi = " | ".join(c["spdi"] if len(c["spdi"]) < 80 else c["spdi"].split(":")[0] + ":" + c["spdi"].split(":")[1] + f":<{c['ref_len']} bp>:{c['spdi'].split(':')[3][:20]}" for c in components)
        uids = match["exact_clinvar_uids"]
        record = vcv.get(uids[0]) if uids else None
        related = RELATED.get(row["sequence_id"])
        related_record = vcv[related[0]] if related else None
        if related and related[1] == "exact":
            location = related_record["location"]
            assert (location["start"], location["stop"]) == tuple(components[0]["hgvs_g"].split("g.")[1].replace("del", "").split("_")) and location["alternateAlleleVCF"] == components[0]["alt_vcf"] and location["referenceAlleleVCF"] == components[0]["ref_vcf"]
            record, uids = related_record, [related[0]]
        partial = bool(not uids and any(match["component_clinvar_uids"].values()))
        rsids = sorted({str(r) for c in components for r in ((dbsnp.get(c["spdi"], {}).get("response") or {}).get("data") or {}).get("rsids", [])})
        joined = match["key"] and dbsnp.get(match["key"], {}).get("response")
        if joined and (joined.get("data") or {}).get("rsids"):
            rsids = sorted(set(rsids) | {str(r) for r in joined["data"]["rsids"]})
        dbsnp_status = "lookup failed" if any(dbsnp.get(c["spdi"], {}).get("response") is None for c in components) else ("found" if rsids else "no rsID")
        gnomad_af = {}
        if len(components) == 1:
            c = components[0]
            variant_id = f"{c['chrom']}-{c['pos_vcf']}-{c['ref_vcf']}-{c['alt_vcf']}"
            entry = gnomad[row["gene"]]["variants"].get(variant_id)
            if entry:
                for part in ("exome", "genome"):
                    if entry.get(part):
                        gnomad_af[part] = (entry[part]["ac"], entry[part]["an"], entry[part]["af"])
        ac = sum(v[0] for v in gnomad_af.values())
        an = sum(v[1] for v in gnomad_af.values())
        effect = vep.get(row["sequence_id"], {})
        category, agreement = categorise(row, record, partial, bool(rsids))
        names, relevance, pathogenic_names = condition_summary(record) if record else ([], "none", [])
        evidence = alzheimer_evidence(record, row["original_label"]) if record else {}
        submitters = sorted({s["submitter"] for s in record["submissions"]}) if record else []
        out = {
            "sequence_id": row["sequence_id"], "gene": row["gene"], "original_label": row["original_label"], "assembly": row["assembly"], "chromosome": components[0]["chrom"],
            "position_vcf": components[0]["pos_vcf"], "reference_allele_vcf": "".join(c["ref_vcf"] for c in components) if len(components) > 1 else (components[0]["ref_vcf"] if components[0]["ref_len"] < 60 else f"<{components[0]['ref_len'] + 1} bp>"),
            "alternate_allele_vcf": "".join(c["alt_vcf"] for c in components) if len(components) > 1 else components[0]["alt_vcf"][:60], "variant_type": row["kind"], "deletion_or_insertion_length": max(c["ref_len"] for c in components) - max(c["alt_len"] for c in components),
            "left_shift_applied": components[0]["shifted_by"], "reference_allele_matches_grch38": all(c["reference_allele_matches"] for c in components), "spdi": spdi,
            "dbsnp_status": dbsnp_status, "rsids": ";".join("rs" + r for r in rsids),
            "clinvar_vcv": record["vcv"] if record else "", "clinvar_match_type": ("exact coordinates and VCF alleles (no SPDI in record)" if related else "exact allele") if record else ("component of MNV only" if partial else "none"),
            "related_vcv": related_record["vcv"] if related and not record else "", "related_classification": related_record["classification"] if related and not record else "", "related_review_status": related_record["review"] if related and not record else "",
            "related_conditions": "; ".join(condition_summary(related_record)[0])[:300] if related and not record else "", "related_note": related[2] if related else "",
            "clinvar_name": record["name"] if record else "", "clinvar_classification": record["classification"] if record else "", "clinvar_review_status": record["review"] if record else "",
            "clinvar_stars": STARS.get(record["review"], "") if record else "", "clinvar_last_evaluated": record["last_evaluated"] if record else "", "clinvar_updated": record["updated"] if record else "",
            "clinvar_submission_count": len(record["submissions"]) if record else 0, "clinvar_submitters": "; ".join(s for s in submitters if s)[:300],
            "clinvar_conditions": "; ".join(names)[:400], "alzheimer_relevance": relevance if record else "no ClinVar record", "ad_evidence_tier": evidence["tier"] if record else "", "alzheimer_class": evidence["kind"] if record else "", "ad_rcv_classifications": evidence["rcvs"] if record else "",
            "submission_disagreement": ";".join(sorted({normal_class(s["classification"]) for s in record["submissions"]} - {"other", "none"})) if record and len({normal_class(s["classification"]) for s in record["submissions"]} - {"other", "none"}) > 1 else "",
            "clinvar_first_created": record["created"] if record else "", "clinvar_most_recent_submission": record["most_recent_submission"] if record else "", "clinvar_retrieved": summaries["retrieved"], "conditions_with_pathogenic_assertion": "; ".join(pathogenic_names)[:300],
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{uids[0]}/" if uids else "",
            "ensembl_consequence": effect.get("most_severe_consequence", "not queried (over 1 kb)" if not effect else ""), "gnomad_r4_ac": ac if gnomad_af else "", "gnomad_r4_an": an if gnomad_af else "", "gnomad_r4_af": round(ac / an, 8) if gnomad_af and an else "",
            "gnomad_status": "present" if gnomad_af else ("not looked up per variant (MNV components only)" if len(components) > 1 else "absent from gnomAD r4 response"),
            "provenance_category": category, "label_agreement": agreement,
        }
        rows.append(out)
        if not record:
            unmatched.append({"sequence_id": row["sequence_id"], "gene": row["gene"], "original_label": row["original_label"], "variant_type": row["kind"], "spdi": spdi, "dbsnp_status": dbsnp_status, "rsids": out["rsids"],
                              "component_clinvar_records": "; ".join(f"{k}: {','.join(v)}" for k, v in match["component_clinvar_uids"].items() if v),
                              "related_vcv": out["related_vcv"], "related_classification": out["related_classification"],
                              "reason": (related[2] if related else "No ClinVar record at this position; dbSNP rsID and gnomAD frequency exist, but no disease assertion") + ("; dbSNP lookup failed (URL too long or HTTP 503)" if dbsnp_status == "lookup failed" else "")})
    for name, table in (("variant_provenance.csv", rows), ("unmatched_variants.csv", unmatched)):
        with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(table[0]))
            writer.writeheader(); writer.writerows(table)

    labelled = defaultdict(Counter)
    for r in rows:
        labelled[r["original_label"]][normal_class(r["clinvar_classification"]) if r["clinvar_vcv"] else "no ClinVar exact match"] += 1
    def af(label):
        values = [r["gnomad_r4_af"] for r in rows if r["original_label"] == label and r["gnomad_r4_af"] != ""]
        return {"present": len(values), "median_af": statistics.median(values) if values else None}
    by_gene_cat = Counter((r["gene"], r["provenance_category"]) for r in rows)
    summary = {
        "total": len(rows), "reference_alleles_ok": sum(r["reference_allele_matches_grch38"] for r in rows),
        "category": Counter(r["provenance_category"] for r in rows), "agreement": Counter(r["label_agreement"] for r in rows),
        "category_by_label": Counter(f"{r['original_label']}|{r['provenance_category']}" for r in rows),
        "by_gene_category": {f"{g}|{c}": n for (g, c), n in by_gene_cat.items()},
        "classification_by_label": {k: dict(v) for k, v in labelled.items()},
        "relevance_by_label": Counter(f"{r['original_label']}|{r['alzheimer_relevance']}" for r in rows),
        "agreement_by_gene": Counter(f"{r['gene']}|{r['label_agreement']}" for r in rows),
        "ad_tier_by_label": Counter(f"{r['original_label']}|{r['ad_evidence_tier']}" for r in rows if r["clinvar_vcv"]),
        "alzheimer_class_by_label": Counter(f"{r['original_label']}|{r['alzheimer_class']}" for r in rows if r["clinvar_vcv"]),
        "submission_disagreement": [(r["sequence_id"], r["clinvar_vcv"], r["submission_disagreement"]) for r in rows if r["submission_disagreement"]],
        "stars_by_label": Counter(f"{r['original_label']}|{r['clinvar_stars']}" for r in rows if r["clinvar_vcv"]),
        "gnomad": {label: af(label) for label in ("HEALTHY", "UNHEALTHY")},
        "consequence_by_label": Counter(f"{r['original_label']}|{r['ensembl_consequence']}" for r in rows),
        "contradicted": [{k: r[k] for k in ("sequence_id", "gene", "original_label", "clinvar_vcv", "clinvar_classification", "clinvar_stars", "clinvar_conditions")} for r in rows if r["label_agreement"] == "contradicted"],
        "multi_exact_uids": {k: v["exact_clinvar_uids"] for k, v in matches.items() if len(v["exact_clinvar_uids"]) > 1},
        "unmatched": unmatched,
        "last_evaluated_range": (min((r["clinvar_last_evaluated"] for r in rows if r["clinvar_last_evaluated"])), max((r["clinvar_last_evaluated"] for r in rows if r["clinvar_last_evaluated"]))),
        "created_after_dataset": [r["sequence_id"] for r in rows if r["clinvar_vcv"] and vcv_field(r, vcv, "created") > "2026-05-08"],
        "submission_after_dataset": [r["sequence_id"] for r in rows if r["clinvar_vcv"] and vcv_field(r, vcv, "most_recent_submission") > "2026-05-08"],
        "locus_clinvar_totals": {g: s["count"] for g, s in summaries["searches"].items()},
        "locus_classification": locus(summaries),
    }
    (RAW / "analysis_summary.json").write_text(json.dumps(summary, indent=1, default=dict), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in summary if k not in ("unmatched",)}, indent=0, default=dict)[:9000])


def vcv_field(row: dict, vcv: dict, field: str) -> str:
    return next((r.get(field, "") for r in vcv.values() if r["vcv"] == row["clinvar_vcv"]), "")


def locus(summaries: dict) -> dict:
    counts = defaultdict(Counter)
    for record in summaries["records"].values():
        counts[record["_locus"]][normal_class(record.get("germline_classification", {}).get("description", ""))] += 1
    return {k: dict(v) for k, v in counts.items()}


if __name__ == "__main__":
    main()
