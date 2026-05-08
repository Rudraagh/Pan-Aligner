from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from common import DATA_DIR, METADATA_DIR, ROOT, ensure_dir, normalize_sequence, parse_sequence_header, validate_dna, write_json, write_wrapped_fasta


def preprocess_fasta_files(input_fastas: list[Path], deduplicate: bool = True) -> dict:
    healthy_root = ensure_dir(DATA_DIR / "healthy")
    unhealthy_root = ensure_dir(DATA_DIR / "unhealthy")
    raw_root = ensure_dir(DATA_DIR / "raw")
    metadata_root = ensure_dir(METADATA_DIR)

    manifest: dict[str, dict] = {
        "input_fastas": [str(path.resolve()) for path in input_fastas],
        "genes": {},
        "summary": {
            "healthy_sequences": 0,
            "unhealthy_sequences": 0,
            "reference_sequences": 0,
            "duplicates_skipped": 0,
            "invalid_sequences_skipped": 0,
        },
    }

    per_gene_records: dict[str, dict[str, list[dict]]] = defaultdict(lambda: {"HEALTHY": [], "UNHEALTHY": [], "REFERENCE": []})
    seen_sequences: dict[tuple[str, str], set[str]] = defaultdict(set)

    for fasta_path in input_fastas:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            metadata = parse_sequence_header(record.id)
            sequence = normalize_sequence(str(record.seq))

            if not validate_dna(sequence):
                manifest["summary"]["invalid_sequences_skipped"] += 1
                continue

            dedup_key = (metadata["gene"], metadata["label"])
            if deduplicate and sequence in seen_sequences[dedup_key]:
                manifest["summary"]["duplicates_skipped"] += 1
                continue

            seen_sequences[dedup_key].add(sequence)
            per_gene_records[metadata["gene"]][metadata["label"]].append(
                {
                    "id": metadata["sequence_id"],
                    "sequence": sequence,
                    "source_file": str(fasta_path.resolve()),
                }
            )
            manifest["summary"][f"{metadata['label'].lower()}_sequences"] += 1

    for gene, buckets in per_gene_records.items():
        if len(buckets["REFERENCE"]) != 1:
            raise ValueError(f"Gene {gene} must have exactly one reference sequence. Found {len(buckets['REFERENCE'])}.")

        gene_manifest = {
            "reference": {},
            "healthy": [],
            "unhealthy": [],
        }

        reference_record = buckets["REFERENCE"][0]
        ref_path = raw_root / f"{gene.lower()}_reference.fa"
        write_wrapped_fasta([(reference_record["id"], reference_record["sequence"])], ref_path)
        gene_manifest["reference"] = {
            "id": reference_record["id"],
            "path": str(ref_path.resolve()),
            "length": len(reference_record["sequence"]),
        }

        for label_name, records, output_root in (
            ("healthy", buckets["HEALTHY"], healthy_root),
            ("unhealthy", buckets["UNHEALTHY"], unhealthy_root),
        ):
            gene_dir = ensure_dir(output_root / gene.lower())
            combined_path = output_root / f"{gene.lower()}_{label_name}.fa"
            combined_records: list[tuple[str, str]] = []

            for item in records:
                sample_path = gene_dir / f"{item['id']}.fa"
                write_wrapped_fasta([(item["id"], item["sequence"])], sample_path)
                combined_records.append((item["id"], item["sequence"]))
                gene_manifest[label_name].append(
                    {
                        "id": item["id"],
                        "path": str(sample_path.resolve()),
                        "length": len(item["sequence"]),
                    }
                )

            write_wrapped_fasta(combined_records, combined_path)
            gene_manifest[f"{label_name}_combined_fasta"] = str(combined_path.resolve())

        combined_all_path = raw_root / f"{gene.lower()}_all_labeled.fa"
        combined_records = [(reference_record["id"], reference_record["sequence"])]
        combined_records.extend((entry["id"], next(item["sequence"] for item in buckets["HEALTHY"] if item["id"] == entry["id"])) for entry in gene_manifest["healthy"])
        combined_records.extend((entry["id"], next(item["sequence"] for item in buckets["UNHEALTHY"] if item["id"] == entry["id"])) for entry in gene_manifest["unhealthy"])
        write_wrapped_fasta(combined_records, combined_all_path)
        gene_manifest["all_labeled_fasta"] = str(combined_all_path.resolve())

        manifest["genes"][gene] = gene_manifest

    manifest_path = metadata_root / "preprocess_manifest.json"
    write_json(manifest_path, manifest)
    return manifest


def default_fasta_inputs() -> list[Path]:
    return sorted(ROOT.glob("*_combined.fasta"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess AD gene FASTA files into healthy/unhealthy/reference partitions.")
    parser.add_argument(
        "--input-fastas",
        nargs="+",
        type=Path,
        default=default_fasta_inputs(),
        help="Input FASTA files. Defaults to all '*_combined.fasta' files in the project root.",
    )
    parser.add_argument("--keep-duplicates", action="store_true", help="Keep duplicate sequences instead of removing them.")
    args = parser.parse_args()

    manifest = preprocess_fasta_files([path.resolve() for path in args.input_fastas], deduplicate=not args.keep_duplicates)
    print(f"Preprocessing complete. Manifest written to {(METADATA_DIR / 'preprocess_manifest.json').resolve()}")
    print(f"Genes processed: {', '.join(sorted(manifest['genes']))}")
    print(f"Healthy sequences: {manifest['summary']['healthy_sequences']}")
    print(f"Unhealthy sequences: {manifest['summary']['unhealthy_sequences']}")


if __name__ == "__main__":
    main()

