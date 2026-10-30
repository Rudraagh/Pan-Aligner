"""Exact-allele matching of the 220 variants to ClinVar, then fetch full VCV records for matches.

ClinVar canonical SPDI is reduced to the same left-aligned form as our variants, using the locus
reference. Also retries gnomAD by 40 kb chunks, because whole-locus queries failed with HTTP 500.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location("query", HERE / "02_query_external.py")
query = importlib.util.module_from_spec(spec)
spec.loader.exec_module(query)
GENE_FILES = {"APP": "app", "PSEN1": "psen1", "PSEN2": "psen2"}
OFFSET = {"NC_000021.9": ("APP", 25880550), "NC_000014.9": ("PSEN1", 73136417), "NC_000001.11": ("PSEN2", 226870184)}


def reference(gene: str) -> str:
    text = (ROOT / "data" / "raw" / f"{GENE_FILES[gene]}_reference.fa").read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in text if not line.startswith(">")).upper()


def normal_key(spdi: str, refs: dict) -> str | None:
    """Trim shared bases and left-align; None when outside the locus or not sequence-resolved."""
    parts = spdi.split(":")
    if len(parts) != 4 or parts[0] not in OFFSET:
        return None
    accession, position, deleted, inserted = parts[0], int(parts[1]), parts[2], parts[3]
    while deleted and inserted and deleted[-1] == inserted[-1]:
        deleted, inserted = deleted[:-1], inserted[:-1]
    while deleted and inserted and deleted[0] == inserted[0]:
        deleted, inserted, position = deleted[1:], inserted[1:], position + 1
    gene, start = OFFSET[accession]
    index = position - (start - 1)
    sequence = refs[gene]
    if index < 0 or index + len(deleted) > len(sequence):
        return None
    if deleted and inserted:
        return f"{accession}:{position}:{deleted}:{inserted}"
    allele = deleted or inserted
    while index > 0 and allele and sequence[index - 1] == allele[-1]:
        allele = allele[-1] + allele[:-1]
        index -= 1
    position = index + start - 1
    return f"{accession}:{position}:{allele}:" if deleted else f"{accession}:{position}::{allele}"


def main() -> None:
    refs = {gene: reference(gene) for gene in GENE_FILES}
    data = json.loads((RAW / "normalized_variants.json").read_text(encoding="utf-8"))
    clinvar = json.loads((RAW / "clinvar_locus_summaries.json").read_text(encoding="utf-8"))["records"]
    index: dict[str, list[str]] = {}
    for uid, record in clinvar.items():
        for variation in record.get("variation_set", []):
            spdi = variation.get("canonical_spdi") or ""
            key = normal_key(spdi, refs) if spdi else None
            if key:
                index.setdefault(key, []).append(uid)
    matches = {}
    for row in data["variants"]:
        keys = query.spdi_keys(row)
        whole = keys[-1] if row["kind"] == "MNV_adjacent" else keys[0]
        whole_key = normal_key(whole, refs)
        component_hits = {k: index.get(normal_key(k, refs), []) for k in keys if k != whole} if row["kind"] == "MNV_adjacent" else {}
        matches[row["sequence_id"]] = {"key": whole_key, "exact_clinvar_uids": sorted(set(index.get(whole_key, []))), "component_clinvar_uids": component_hits}
    uids = sorted({u for m in matches.values() for u in m["exact_clinvar_uids"]} | {u for m in matches.values() for ids in m["component_clinvar_uids"].values() for u in ids})
    print("matched sequences:", sum(bool(m["exact_clinvar_uids"]) for m in matches.values()), "of", len(matches), "records to fetch:", len(uids), flush=True)
    (RAW / "clinvar_matches.json").write_text(json.dumps(matches, indent=0), encoding="utf-8")

    vcv = {}
    for offset in range(0, len(uids), 40):
        batch = uids[offset:offset + 40]
        url = f"{query.EUTILS}/efetch.fcgi"
        body = query.request(url, data=urllib.parse.urlencode({"db": "clinvar", "rettype": "vcv", "is_variationid": "true", "id": ",".join(batch)}).encode())
        (RAW / f"clinvar_vcv_batch_{offset // 40:03d}.xml").write_text(body or "", encoding="utf-8")
        vcv[offset // 40] = {"ids": batch, "ok": body is not None}
    (RAW / "clinvar_vcv_index.json").write_text(json.dumps({"retrieved": query.now(), "batches": vcv}), encoding="utf-8")
    print("VCV batches:", len(vcv), "failed:", sum(not v["ok"] for v in vcv.values()), flush=True)

    gql = """query($chrom:String!,$start:Int!,$stop:Int!){region(chrom:$chrom,start:$start,stop:$stop,reference_genome:GRCh38){
      variants(dataset:gnomad_r4){variant_id rsids exome{ac an af} genome{ac an af}}}}"""
    saved = json.loads((RAW / "gnomad_regions.json").read_text(encoding="utf-8"))
    for gene in ("APP", "PSEN1"):
        chrom, start, stop = query.LOCI[gene]
        variants, errors = {}, []
        for chunk_start in range(start, stop + 1, 40000):
            chunk_stop = min(chunk_start + 39999, stop)
            body = query.request("https://gnomad.broadinstitute.org/api", data=json.dumps({"query": gql, "variables": {"chrom": chrom, "start": chunk_start, "stop": chunk_stop}}).encode(),
                                 headers={"Content-Type": "application/json"}, pause=8.0)
            parsed = json.loads(body) if body else None
            if not parsed or parsed.get("errors") or not parsed.get("data", {}).get("region"):
                errors.append({"start": chunk_start, "stop": chunk_stop, "error": parsed.get("errors") if parsed else "request failed"})
                continue
            for variant in parsed["data"]["region"]["variants"]:
                variants[variant["variant_id"]] = variant
        saved[gene] = {"dataset": "gnomad_r4", "retrieved": query.now(), "errors": errors or None, "variants": variants, "method": "40 kb chunks"}
        print(gene, "gnomAD variants:", len(variants), "failed chunks:", len(errors), flush=True)
    (RAW / "gnomad_regions.json").write_text(json.dumps(saved), encoding="utf-8")
    path = RAW / "request_failures.json"
    existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    path.write_text(json.dumps(existing + query.failures, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
