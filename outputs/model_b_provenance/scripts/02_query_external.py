"""Query ClinVar, NCBI Variation Services (dbSNP), Ensembl VEP and gnomAD; cache raw responses.

Only public variant coordinates are sent. Every response is stored under raw/ with its URL and
retrieval time, and failed requests are recorded, never treated as empty results.
"""
from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "raw"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
LOCI = {"APP": ("21", 25880550, 26171128), "PSEN1": ("14", 73136417, 73223691), "PSEN2": ("1", 226870184, 226927726)}
failures: list[dict] = []


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def request(url: str, data: bytes | None = None, headers: dict | None = None, attempts: int = 4, pause: float = 0.4) -> str | None:
    for attempt in range(attempts):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "fyp-provenance-audit", **(headers or {})})
            with urllib.request.urlopen(req, timeout=180) as response:
                body = response.read().decode("utf-8")
            time.sleep(pause)
            return body
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                time.sleep(pause)
                return json.dumps({"http_404": True, "body": exc.read().decode("utf-8", "replace")[:500]})
            last = f"HTTP {exc.code}"
        except Exception as exc:  # network errors are recorded, not hidden
            last = repr(exc)
        time.sleep(2 * (attempt + 1))
    failures.append({"url": url[:300], "error": last, "time": now()})
    return None


def clinvar() -> None:
    all_records: dict[str, dict] = {}
    searches = {}
    for gene, (chrom, start, stop) in LOCI.items():
        term = f"{chrom}[chr] AND {start}:{stop}[chrpos38]"
        url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode({"db": "clinvar", "term": term, "retmode": "json", "retmax": 100000})
        body = request(url)
        result = json.loads(body)["esearchresult"]
        ids = result["idlist"]
        searches[gene] = {"term": term, "count": int(result["count"]), "returned": len(ids), "url": url, "retrieved": now()}
        print(gene, "ClinVar records in locus:", result["count"], flush=True)
        for offset in range(0, len(ids), 200):
            batch = ids[offset:offset + 200]
            body = request(f"{EUTILS}/esummary.fcgi", data=urllib.parse.urlencode({"db": "clinvar", "id": ",".join(batch), "retmode": "json"}).encode())
            if body is None:
                continue
            result_set = json.loads(body)["result"]
            for uid in result_set["uids"]:
                all_records[uid] = {**result_set[uid], "_locus": gene}
    (RAW / "clinvar_locus_summaries.json").write_text(json.dumps({"retrieved": now(), "searches": searches, "records": all_records}), encoding="utf-8")
    print("ClinVar summaries:", len(all_records), flush=True)


def dbsnp(variants: list[dict]) -> None:
    out = {}
    for row in variants:
        for key in spdi_keys(row):
            url = f"https://api.ncbi.nlm.nih.gov/variation/v0/spdi/{urllib.parse.quote(key, safe=':')}/rsids"
            body = request(url, pause=0.35)
            out[key] = {"response": json.loads(body) if body else None, "url": url, "retrieved": now()}
    (RAW / "dbsnp_spdi_rsids.json").write_text(json.dumps(out), encoding="utf-8")
    print("dbSNP lookups:", len(out), "with rsIDs:", sum(bool(v["response"] and v["response"].get("data", {}).get("rsids")) for v in out.values()), flush=True)


def spdi_keys(row: dict) -> list[str]:
    """Each component, plus the joined delins for adjacent MNVs, which databases may hold as one record."""
    keys = [c["spdi"] for c in row["components"]]
    if row["kind"] == "MNV_adjacent":
        first = row["components"][0]["spdi"].split(":")
        refs = "".join(c["spdi"].split(":")[2] for c in row["components"])
        alts = "".join(c["spdi"].split(":")[3] for c in row["components"])
        keys.append(f"{first[0]}:{first[1]}:{refs}:{alts}")
    return keys


def ensembl(variants: list[dict]) -> None:
    lines = {}
    for row in variants:
        if row["kind"] == "MNV_adjacent":
            parts = row["components"]
            lines[row["sequence_id"]] = f"{parts[0]['chrom']} {parts[0]['pos_vcf']} {row['sequence_id']} {''.join(p['ref_vcf'] for p in parts)} {''.join(p['alt_vcf'] for p in parts)} . . ."
        else:
            c = row["components"][0]
            if c["ref_len"] > 1000:
                continue  # very long deletion: looked up by coordinates separately in the report
            lines[row["sequence_id"]] = f"{c['chrom']} {c['pos_vcf']} {row['sequence_id']} {c['ref_vcf']} {c['alt_vcf']} . . ."
    ids = list(lines)
    results = {}
    for offset in range(0, len(ids), 150):
        batch = [lines[i] for i in ids[offset:offset + 150]]
        body = request("https://rest.ensembl.org/vep/homo_sapiens/region?content-type=application/json&variant_class=1&hgvs=1&mane=1&pick=0",
                       data=json.dumps({"variants": batch}).encode(), headers={"Content-Type": "application/json", "Accept": "application/json"}, pause=1.0)
        if body:
            for item in json.loads(body):
                results[item.get("id") or item.get("input")] = item
    (RAW / "ensembl_vep.json").write_text(json.dumps({"retrieved": now(), "results": results, "inputs": lines}), encoding="utf-8")
    print("Ensembl VEP results:", len(results), "of", len(lines), flush=True)


def gnomad() -> None:
    query = """query($chrom:String!,$start:Int!,$stop:Int!){region(chrom:$chrom,start:$start,stop:$stop,reference_genome:GRCh38){
      variants(dataset:gnomad_r4){variant_id rsids exome{ac an af} genome{ac an af}}}}"""
    out = {}
    for gene, (chrom, start, stop) in LOCI.items():
        body = request("https://gnomad.broadinstitute.org/api", data=json.dumps({"query": query, "variables": {"chrom": chrom, "start": start, "stop": stop}}).encode(),
                       headers={"Content-Type": "application/json"}, pause=8.0)
        parsed = json.loads(body) if body else None
        errors = parsed.get("errors") if parsed else "request failed"
        variants = parsed["data"]["region"]["variants"] if parsed and parsed.get("data") and parsed["data"].get("region") else []
        out[gene] = {"dataset": "gnomad_r4", "retrieved": now(), "errors": errors, "variants": {v["variant_id"]: v for v in variants}}
        print(gene, "gnomAD variants:", len(variants), "errors:", str(errors)[:200], flush=True)
    (RAW / "gnomad_regions.json").write_text(json.dumps(out), encoding="utf-8")


def main() -> None:
    steps = sys.argv[1:] or ["clinvar", "dbsnp", "ensembl", "gnomad"]
    variants = json.loads((RAW / "normalized_variants.json").read_text(encoding="utf-8"))["variants"]
    if "clinvar" in steps:
        clinvar()
    if "dbsnp" in steps:
        dbsnp(variants)
    if "ensembl" in steps:
        ensembl(variants)
    if "gnomad" in steps:
        gnomad()
    path = RAW / "request_failures.json"
    existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    path.write_text(json.dumps(existing + failures, indent=1), encoding="utf-8")
    print("failed requests:", len(failures), flush=True)


if __name__ == "__main__":
    main()
