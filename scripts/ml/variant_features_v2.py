"""Features for the replacement Model B (variant clinical classification), built blind to the labels.

Only identity columns of the frozen label table are read, so no feature can depend on a class. Three tables:

  features_main.csv      every variant: variant type, normalized substitution, sequence context
  features_missense.csv  missense variants only: amino-acid change properties, from Ensembl VEP
  baselines.csv          gene, consequence, position; for baselines and stratifying, never for the main model

Amino-acid properties are a separate table on purpose: filling them only for missense variants would let a
pooled model read consequence from which cells are empty. Usage: fetch, then build (see --help). No model is fitted.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FREEZE = ROOT / "outputs" / "model_b_task1_freeze" / "20261030_freeze_01" / "data" / "task1_labels.csv"
# The only columns ever read from the label table.
IDENTITY_COLUMNS = ("variation_id", "gene", "chrom", "pos", "ref", "alt")
BASELINE_INPUT_COLUMNS = ("variation_id", "consequence_group")
ACCESSION = {"21": "NC_000021.9", "14": "NC_000014.9", "1": "NC_000001.11"}
PADDING = 300
COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Kyte & Doolittle 1982, hydropathy.
HYDROPATHY = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
              "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# Zamyatnin 1972, residue volume in cubic angstroms.
VOLUME = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7,
          "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
THREE_TO_ONE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
                "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PURINES = {"A", "G"}

MAIN_COLUMNS = ["variation_id", "variant_type", "length_change", "is_transition", "substitution", "trinucleotide", "in_cpg", "gc_50", "gc_200", "homopolymer_run"]
MISSENSE_COLUMNS = ["variation_id", "aa_ref", "aa_alt", "blosum62", "delta_hydropathy", "delta_volume", "delta_charge", "to_proline", "from_glycine", "from_cysteine", "to_cysteine"]
BASELINE_COLUMNS = ["variation_id", "gene", "chrom", "pos", "relative_position", "consequence_group", "vep_consequence", "protein_position"]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "snv"
    if len(ref) == len(alt):
        return "mnv"
    return "deletion" if len(ref) > len(alt) else "insertion"


def substitution_and_context(ref: str, alt: str, before: str, after: str) -> tuple[str, str, str]:
    """(normalized substitution, normalized trinucleotide, transition flag). Strand-neutral, so the gene's strand adds nothing."""
    if len(ref) != 1 or len(alt) != 1:
        return "none", "none", ""
    transition = str(int({ref, alt} in ({"A", "G"}, {"C", "T"})))
    if ref in PURINES:
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
        before, after = after.translate(COMPLEMENT), before.translate(COMPLEMENT)
    return f"{ref}>{alt}", f"{before}{ref}{after}", transition


def in_cpg(sequence: str, start: int, end: int) -> int:
    """Whether any base of [start, end) belongs to a CG dinucleotide. Coordinates index `sequence`."""
    return int("CG" in sequence[max(0, start - 1):end + 1])


def gc_fraction(sequence: str, start: int, end: int, flank: int) -> float:
    """GC of the reference flanks around [start, end), with the variant's own bases left out."""
    window = sequence[max(0, start - flank):start] + sequence[end:end + flank]
    return round(sum(base in "GC" for base in window) / len(window), 4) if window else 0.0


def homopolymer_run(sequence: str, start: int, end: int) -> int:
    """Longest single-base run starting at the variant's first or last reference base, extending outward."""
    best = 1
    for anchor, step in ((start, -1), (end - 1, 1)):
        base = sequence[anchor]
        length, position = 1, anchor + step
        while 0 <= position < len(sequence) and sequence[position] == base:
            length += 1
            position += step
        best = max(best, length)
    return best


def main_features(row: dict, sequence: str, offset: int) -> dict:
    ref, alt = row["ref"], row["alt"]
    start = int(row["pos"]) - 1 - offset
    end = start + len(ref)
    if sequence[start:end] != ref:
        raise ValueError(f"reference mismatch for {row['variation_id']}: file {sequence[start:end]}, record {ref}")
    substitution, trinucleotide, transition = substitution_and_context(ref, alt, sequence[start - 1], sequence[end])
    return {"variation_id": row["variation_id"], "variant_type": variant_type(ref, alt), "length_change": len(alt) - len(ref),
            "is_transition": transition, "substitution": substitution, "trinucleotide": trinucleotide, "in_cpg": in_cpg(sequence, start, end),
            "gc_50": gc_fraction(sequence, start, end, 50), "gc_200": gc_fraction(sequence, start, end, 200), "homopolymer_run": homopolymer_run(sequence, start, end)}


def single_amino_change(amino_acids: str) -> tuple[str, str] | None:
    """(ref, alt) from VEP's form, such as `R/K` or `VI/VV` for a change spanning two codons; None if not one residue."""
    parts = amino_acids.split("/")
    if len(parts) != 2 or len(parts[0]) != len(parts[1]):
        return None
    changed = [(a, b) for a, b in zip(*parts) if a != b]
    if len(changed) != 1 or changed[0][0] not in VOLUME or changed[0][1] not in VOLUME:
        return None
    return changed[0]


def is_no_change(ref: str, alt: str) -> bool:
    """A record such as `c.1090+194=` asserts the reference allele, so it is not a variant."""
    return ref == alt


def missense_features(variation_id: str, amino_acids: str, blosum) -> dict:
    """`amino_acids` is VEP's form; the caller has already checked it is a single-residue change."""
    ref, alt = single_amino_change(amino_acids)
    return {"variation_id": variation_id, "aa_ref": ref, "aa_alt": alt, "blosum62": int(blosum[ref][alt]),
            "delta_hydropathy": round(HYDROPATHY[alt] - HYDROPATHY[ref], 2), "delta_volume": round(VOLUME[alt] - VOLUME[ref], 1),
            "delta_charge": CHARGE.get(alt, 0) - CHARGE.get(ref, 0), "to_proline": int(alt == "P"), "from_glycine": int(ref == "G"),
            "from_cysteine": int(ref == "C"), "to_cysteine": int(alt == "C")}


def canonical_consequence(item: dict | None, gene: str) -> dict | None:
    """The canonical-transcript consequence for this gene, if VEP returned one."""
    for transcript in (item or {}).get("transcript_consequences", []):
        if transcript.get("canonical") and transcript.get("gene_symbol") == gene:
            return transcript
    return None


def read_columns(path: Path, columns: tuple[str, ...]) -> list[dict]:
    """Read named columns only, so labels never enter memory here."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return [{column: row[column] for column in columns} for row in csv.DictReader(handle)]


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def request(url: str, data: bytes | None = None, headers: dict | None = None) -> bytes:
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "fyp-model-b-features", **(headers or {})})
            with urllib.request.urlopen(req, timeout=180) as response:
                return response.read()
        except Exception as error:  # noqa: BLE001 - retried; the last failure is raised
            if attempt == 3:
                raise RuntimeError(f"request failed: {url[:120]}") from error
            time.sleep(3 * (attempt + 1))
    raise AssertionError("unreachable")


def fetch(freeze: Path, output: Path) -> None:
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=False)
    rows = read_columns(freeze, IDENTITY_COLUMNS)
    log = {"retrieved": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "references": {}, "vep": {}}
    for chrom, accession in ACCESSION.items():
        positions = [int(row["pos"]) for row in rows if row["chrom"] == chrom]
        start, stop = min(positions) - PADDING, max(positions) + PADDING
        query = urllib.parse.urlencode({"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text", "seq_start": start, "seq_stop": stop})
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}"
        body = request(url)
        path = raw / f"reference_chr{chrom}.fa"
        path.write_bytes(body)
        log["references"][chrom] = {"accession": accession, "assembly": "GRCh38", "start_1based": start, "stop_1based": stop, "url": url, "sha256": sha256(path)}
        time.sleep(0.5)
    results = []
    strings = [f"{r['chrom']} {r['pos']} {r['variation_id']} {r['ref']} {r['alt']} . . ." for r in rows]
    for index in range(0, len(strings), 200):
        body = json.dumps({"variants": strings[index:index + 200], "canonical": 1}).encode()
        results += json.loads(request("https://rest.ensembl.org/vep/homo_sapiens/region", data=body, headers={"content-type": "application/json", "accept": "application/json"}))
        time.sleep(1)
    info = json.loads(request("https://rest.ensembl.org/info/data?content-type=application/json"))
    path = raw / "ensembl_vep.json.gz"
    path.write_bytes(gzip.compress(json.dumps(results).encode(), mtime=0))
    log["vep"] = {"endpoint": "https://rest.ensembl.org/vep/homo_sapiens/region", "release": info, "assembly": "GRCh38", "requested": len(strings), "returned": len(results), "sha256": sha256(path)}
    (raw / "fetch_log.json").write_text(json.dumps(log, indent=1), encoding="utf-8")
    print(json.dumps({"references": {k: (v["start_1based"], v["stop_1based"]) for k, v in log["references"].items()}, "vep_returned": len(results), "requested": len(strings)}))


def load_reference(raw: Path, chrom: str) -> tuple[str, int]:
    log = json.loads((raw / "fetch_log.json").read_text(encoding="utf-8"))
    lines = (raw / f"reference_chr{chrom}.fa").read_text(encoding="ascii").splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
    return sequence, log["references"][chrom]["start_1based"] - 1


def write_table(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(freeze: Path, raw: Path, output: Path) -> None:
    from Bio.Align import substitution_matrices

    blosum = substitution_matrices.load("BLOSUM62")
    output.mkdir(parents=True, exist_ok=False)
    rows = read_columns(freeze, IDENTITY_COLUMNS)
    consequence = {row["variation_id"]: row["consequence_group"] for row in read_columns(freeze, BASELINE_INPUT_COLUMNS)}
    vep = {item["id"]: item for item in json.loads(gzip.decompress((raw / "ensembl_vep.json.gz").read_bytes()))}
    references = {chrom: load_reference(raw, chrom) for chrom in ACCESSION}
    positions_by_gene: dict[str, list[int]] = {}
    for row in rows:
        positions_by_gene.setdefault(row["gene"], []).append(int(row["pos"]))

    main_rows, missense_rows, baseline_rows = [], [], []
    no_canonical, skipped_missense, excluded = [], [], []
    for row in rows:
        if is_no_change(row["ref"], row["alt"]):
            excluded.append({"variation_id": row["variation_id"], "reason": "REF equals ALT: reference-allele record, not a variant"})
            continue
        sequence, offset = references[row["chrom"]]
        main_rows.append(main_features(row, sequence, offset))
        transcript = canonical_consequence(vep.get(row["variation_id"]), row["gene"])
        if transcript is None:
            no_canonical.append(row["variation_id"])
        terms = (transcript or {}).get("consequence_terms", [])
        if "missense_variant" in terms:
            amino = transcript.get("amino_acids", "")
            if single_amino_change(amino):
                missense_rows.append(missense_features(row["variation_id"], amino, blosum))
            else:
                skipped_missense.append(row["variation_id"])
        low, high = min(positions_by_gene[row["gene"]]), max(positions_by_gene[row["gene"]])
        baseline_rows.append({"variation_id": row["variation_id"], "gene": row["gene"], "chrom": row["chrom"], "pos": row["pos"],
                              "relative_position": round((int(row["pos"]) - low) / (high - low), 6), "consequence_group": consequence[row["variation_id"]],
                              "vep_consequence": ";".join(terms), "protein_position": (transcript or {}).get("protein_start", "")})
    write_table(output / "features_main.csv", MAIN_COLUMNS, main_rows)
    write_table(output / "features_missense.csv", MISSENSE_COLUMNS, missense_rows)
    write_table(output / "baselines.csv", BASELINE_COLUMNS, baseline_rows)
    write_table(output / "excluded_variants.csv", ["variation_id", "reason"], excluded)
    import Bio

    manifest = {"dataset": "variant_pathogenicity_v1 features", "freeze": str(freeze.relative_to(ROOT)).replace("\\", "/"), "freeze_sha256": sha256(freeze),
                "label_table_columns_read": [*IDENTITY_COLUMNS, "consequence_group"], "models_fitted": 0, "label_feature_associations_computed": False,
                "python": platform.python_version(), "biopython": Bio.__version__, "rows": {"main": len(main_rows), "missense": len(missense_rows), "baselines": len(baseline_rows)},
                "excluded_from_freeze": excluded, "no_canonical_transcript": no_canonical, "missense_without_simple_amino_acids": skipped_missense,
                "outputs_sha256": {name: sha256(output / name) for name in ("features_main.csv", "features_missense.csv", "baselines.csv", "excluded_variants.csv")}}
    (output / "feature_manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps(manifest["rows"]), "no canonical:", no_canonical, "skipped missense:", skipped_missense, "excluded:", excluded)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("fetch", "build"))
    parser.add_argument("--freeze", type=Path, default=DEFAULT_FREEZE)
    parser.add_argument("--raw", type=Path, help="build: the raw folder written by fetch.")
    parser.add_argument("--output", type=Path, required=True, help="A new folder; existing folders are refused.")
    args = parser.parse_args(argv)
    if args.command == "fetch":
        fetch(args.freeze, args.output)
    else:
        if args.raw is None:
            parser.error("build needs --raw")
        build(args.freeze.resolve(), args.raw, args.output)


if __name__ == "__main__":
    main()
