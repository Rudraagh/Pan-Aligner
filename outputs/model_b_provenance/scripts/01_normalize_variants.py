"""Derive one normalized record per labelled sequence by comparing it with its locus reference.

Read-only with respect to the dataset. Writes raw/normalized_variants.json inside
outputs/model_b_provenance/. Coordinates: FASTA region headers are NCBI 1-based inclusive.
Indels are left-aligned within the locus reference and written as VCF-style records with an
anchor base, as HGVS-like g. strings, and as NCBI SPDI (0-based position).
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "raw"
OUT.mkdir(parents=True, exist_ok=True)
GENES = {"APP": "app", "PSEN1": "psen1", "PSEN2": "psen2"}
CHROM = {"NC_000021.9": "21", "NC_000014.9": "14", "NC_000001.11": "1"}


def read_fasta(path: Path):
    records, name, chunks = [], None, []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    records.append((name, "".join(chunks).upper()))
    return records


def diff(reference: str, sequence: str) -> list[dict]:
    """Differences as (0-based reference start, ref allele, alt allele) on the locus reference."""
    if len(reference) == len(sequence):
        return [{"start": i, "ref": a, "alt": b} for i, (a, b) in enumerate(zip(reference, sequence)) if a != b]
    limit = min(len(reference), len(sequence))
    prefix = 0
    while prefix < limit and reference[prefix] == sequence[prefix]:
        prefix += 1
    suffix = 0
    while suffix < limit - prefix and reference[-1 - suffix] == sequence[-1 - suffix]:
        suffix += 1
    return [{"start": prefix, "ref": reference[prefix:len(reference) - suffix], "alt": sequence[prefix:len(sequence) - suffix]}]


def left_align(reference: str, start: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Shift a pure insertion or deletion left as far as the locus reference allows."""
    if ref and alt:
        return start, ref, alt
    allele = ref or alt
    while start > 0 and reference[start - 1] == allele[-1]:
        allele = allele[-1] + allele[:-1]
        start -= 1
    return (start, allele, "") if ref else (start, "", allele)


def kind(parts: list[dict]) -> str:
    if len(parts) > 1:
        adjacent = all(b["start"] == a["start"] + 1 for a, b in zip(parts, parts[1:]))
        return "MNV_adjacent" if adjacent else "multi_SNV_separated"
    ref, alt = parts[0]["ref"], parts[0]["alt"]
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if not ref:
        return "insertion"
    if not alt:
        return "deletion"
    return "complex"


def main() -> None:
    rows, checks = [], {}
    for gene, stem in GENES.items():
        records = read_fasta(ROOT / "data" / "raw" / f"{stem}_all_labeled.fa")
        header, reference = next((n, s) for n, s in records if n.startswith("ref|"))
        match = re.match(r"ref\|(NC_\d+\.\d+)\|:(\d+)-(\d+)", header)
        accession, region_start, region_end = match.group(1), int(match.group(2)), int(match.group(3))
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text&seq_start={region_start}&seq_stop={region_end}"
        with urllib.request.urlopen(url, timeout=120) as response:
            fetched = "".join(line for line in response.read().decode().splitlines() if not line.startswith(">")).upper()
        checks[gene] = {"header": header, "accession": accession, "region_start": region_start, "region_end": region_end, "length": len(reference),
                        "matches_ncbi_efetch": fetched == reference, "efetch_url": url}
        seen = {}
        for name, sequence in records:
            if name.startswith("ref|"):
                continue
            label = "HEALTHY" if "_H_" in name else "UNHEALTHY"
            parts = diff(reference, sequence)
            variant_kind = kind(parts) if parts else "identical_to_reference"
            components = []
            for part in parts:
                start, ref, alt = left_align(reference, part["start"], part["ref"], part["alt"])
                position0 = region_start - 1 + start
                if ref and alt and len(ref) == 1 and len(alt) == 1:
                    vcf_pos, vcf_ref, vcf_alt = position0 + 1, ref, alt
                    hgvs = f"{accession}:g.{vcf_pos}{ref}>{alt}"
                else:
                    anchor = reference[start - 1] if start else ""
                    vcf_pos, vcf_ref, vcf_alt = position0, anchor + ref, anchor + alt
                    if not alt:
                        hgvs = f"{accession}:g.{position0 + 1}del" if len(ref) == 1 else f"{accession}:g.{position0 + 1}_{position0 + len(ref)}del"
                    elif not ref:
                        hgvs = f"{accession}:g.{position0}_{position0 + 1}ins{alt}"
                    else:
                        hgvs = f"{accession}:g.{position0 + 1}_{position0 + len(ref)}delins{alt}"
                components.append({
                    "locus_start0": start, "shifted_by": part["start"] - start, "chrom": CHROM[accession], "pos_vcf": vcf_pos, "ref_vcf": vcf_ref, "alt_vcf": vcf_alt,
                    "hgvs_g": hgvs, "spdi": f"{accession}:{position0}:{ref}:{alt}", "ref_len": len(ref), "alt_len": len(alt),
                    "reference_allele_matches": reference[start:start + len(ref)] == ref,
                })
            key = tuple(c["spdi"] for c in components)
            seen.setdefault(key, []).append(name)
            rows.append({"sequence_id": name, "gene": gene, "original_label": label, "assembly": "GRCh38.p14", "accession": accession,
                         "kind": variant_kind, "components": components, "sequence_length": len(sequence)})
        checks[gene]["duplicate_variant_sets"] = {" ; ".join(k): v for k, v in seen.items() if len(v) > 1}
    (OUT / "normalized_variants.json").write_text(json.dumps({"reference_checks": checks, "variants": rows}, indent=1), encoding="utf-8")
    print(json.dumps({g: {k: v for k, v in c.items() if k != "efetch_url"} for g, c in checks.items()}, indent=1)[:1500])
    from collections import Counter
    print(Counter((r["gene"], r["original_label"], r["kind"]) for r in rows))
    print("shifted", [(r["sequence_id"], c["shifted_by"], c["spdi"][:60]) for r in rows for c in r["components"] if c["shifted_by"]])
    print("nonSNV", [(r["sequence_id"], r["kind"], [c["spdi"][:50] for c in r["components"]]) for r in rows if r["kind"] != "SNV"])


if __name__ == "__main__":
    main()
