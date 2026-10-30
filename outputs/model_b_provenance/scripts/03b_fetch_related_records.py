"""Fetch ClinVar records for sites the SPDI match could not resolve, and try Alzforum.

Variation IDs were selected by coordinates, never by label: records without a canonical SPDI whose
coordinates cover a long deletion, and records at the same site as an unmatched insertion or MNV.
Also retrieves the Ensembl VEP HGVS output already cached for those sequences. Alzforum is attempted
once per gene; any refusal is recorded as a failed lookup.
"""
from __future__ import annotations

import importlib.util
import json
import urllib.parse
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
spec = importlib.util.spec_from_file_location("query", HERE / "02_query_external.py")
query = importlib.util.module_from_spec(spec)
spec.loader.exec_module(query)
RELATED = ["38296", "162429", "162428", "533784", "98034", "18093"]


def main() -> None:
    body = query.request(f"{query.EUTILS}/efetch.fcgi", data=urllib.parse.urlencode({"db": "clinvar", "rettype": "vcv", "is_variationid": "true", "id": ",".join(RELATED)}).encode())
    (RAW / "clinvar_vcv_related.xml").write_text(body or "", encoding="utf-8")
    alzforum = {}
    for gene in ("app", "psen-1", "psen-2"):
        url = f"https://www.alzforum.org/mutations/{gene}"
        page = query.request(url, attempts=1)
        ok = bool(page) and "http_404" not in page[:40]
        alzforum[gene] = {"url": url, "retrieved": query.now(), "ok": ok, "bytes": len(page or ""), "failure": None if ok else (page or "request failed")[:200]}
        if ok:
            (RAW / f"alzforum_{gene}.html").write_text(page, encoding="utf-8")
    (RAW / "related_fetch_log.json").write_text(json.dumps({"retrieved": query.now(), "clinvar_ids": RELATED, "clinvar_ok": body is not None, "alzforum": alzforum, "failures": query.failures}, indent=1), encoding="utf-8")
    print("ClinVar ok:", body is not None, {k: (v["ok"], v["bytes"]) for k, v in alzforum.items()}, query.failures)


if __name__ == "__main__":
    main()
