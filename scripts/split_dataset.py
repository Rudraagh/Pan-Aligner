from __future__ import annotations

import argparse
import random
from pathlib import Path

from common import DATA_DIR, METADATA_DIR, ensure_dir, read_json, write_json, write_wrapped_fasta


def resolve_local_fasta_path(entry: dict, gene: str, payload_key: str) -> Path:
    candidate = Path(entry["path"])
    if candidate.exists():
        return candidate
    fallback = DATA_DIR / payload_key / gene.lower() / f"{entry['id']}.fa"
    if fallback.exists():
        return fallback
    raise FileNotFoundError(f"Could not locate FASTA for {entry['id']}. Checked {candidate} and {fallback}.")


def split_indices(record_count: int, test_fraction: float, rng: random.Random) -> tuple[list[int], list[int]]:
    if record_count < 2:
        raise ValueError(f"Need at least 2 sequences per class to create train/test splits. Found {record_count}.")
    indices = list(range(record_count))
    rng.shuffle(indices)
    test_count = max(1, round(record_count * test_fraction))
    test_count = min(test_count, record_count - 1)
    test_indices = sorted(indices[:test_count])
    train_indices = sorted(indices[test_count:])
    return train_indices, test_indices


def write_split_records(records: list[dict], destination_dir: Path, combined_path: Path) -> list[dict]:
    written_records: list[dict] = []
    combined_records: list[tuple[str, str]] = []
    ensure_dir(destination_dir)

    for record in records:
        sequence = record["sequence"]
        sample_path = destination_dir / f"{record['id']}.fa"
        write_wrapped_fasta([(record["id"], sequence)], sample_path)
        combined_records.append((record["id"], sequence))
        written_records.append(
            {
                "id": record["id"],
                "path": str(sample_path.resolve()),
                "length": len(sequence),
                "source_file": record["source_file"],
                "label": record["label"],
            }
        )

    write_wrapped_fasta(combined_records, combined_path)
    return written_records


def build_dataset_manifest(split_name: str, preprocess_manifest: dict, split_records: dict[str, dict[str, list[dict]]]) -> dict:
    dataset_manifest = {
        "split_name": split_name,
        "genes": {},
    }

    for gene, gene_payload in preprocess_manifest["genes"].items():
        reference_path = Path(gene_payload["reference"]["path"])
        if not reference_path.exists():
            reference_path = DATA_DIR / "raw" / f"{gene.lower()}_reference.fa"
        dataset_manifest["genes"][gene] = {
            "reference": {
                **gene_payload["reference"],
                "path": str(reference_path.resolve()),
            },
            "healthy": [],
            "unhealthy": [],
        }

        healthy_output_root = DATA_DIR / split_name / "healthy" / gene.lower()
        unhealthy_output_root = DATA_DIR / split_name / "unhealthy" / gene.lower()
        healthy_combined_path = DATA_DIR / split_name / "healthy" / f"{gene.lower()}_healthy.fa"
        unhealthy_combined_path = DATA_DIR / split_name / "unhealthy" / f"{gene.lower()}_unhealthy.fa"

        dataset_manifest["genes"][gene]["healthy"] = write_split_records(
            split_records[gene]["HEALTHY"],
            healthy_output_root,
            healthy_combined_path,
        )
        dataset_manifest["genes"][gene]["unhealthy"] = write_split_records(
            split_records[gene]["UNHEALTHY"],
            unhealthy_output_root,
            unhealthy_combined_path,
        )
        dataset_manifest["genes"][gene]["healthy_combined_fasta"] = str(healthy_combined_path.resolve())
        dataset_manifest["genes"][gene]["unhealthy_combined_fasta"] = str(unhealthy_combined_path.resolve())

    return dataset_manifest


def create_split_summary(split_manifest: dict) -> str:
    lines = [
        "PanAligner AD Pipeline Dataset Split Summary",
        "===========================================",
        f"Random seed: {split_manifest['random_seed']}",
        f"Test fraction: {split_manifest['test_fraction']:.2f}",
        f"Train fraction: {1.0 - split_manifest['test_fraction']:.2f}",
        "",
    ]

    total_train_healthy = total_train_unhealthy = 0
    total_test_healthy = total_test_unhealthy = 0

    for gene, payload in split_manifest["genes"].items():
        total = payload["total_healthy"] + payload["total_unhealthy"]
        lines.extend(
            [
                f"{gene}",
                f"  Total sequences: {total}",
                f"  Healthy train count: {payload['train_healthy']}",
                f"  Healthy test count: {payload['test_healthy']}",
                f"  Unhealthy train count: {payload['train_unhealthy']}",
                f"  Unhealthy test count: {payload['test_unhealthy']}",
                f"  Split percentage: {payload['train_fraction']:.2%} train / {payload['test_fraction']:.2%} test",
                "",
            ]
        )
        total_train_healthy += payload["train_healthy"]
        total_train_unhealthy += payload["train_unhealthy"]
        total_test_healthy += payload["test_healthy"]
        total_test_unhealthy += payload["test_unhealthy"]

    lines.extend(
        [
            "Overall",
            f"  Total sequences: {total_train_healthy + total_train_unhealthy + total_test_healthy + total_test_unhealthy}",
            f"  Healthy train count: {total_train_healthy}",
            f"  Healthy test count: {total_test_healthy}",
            f"  Unhealthy train count: {total_train_unhealthy}",
            f"  Unhealthy test count: {total_test_unhealthy}",
            f"  Split percentage: {(1.0 - split_manifest['test_fraction']):.2%} train / {split_manifest['test_fraction']:.2%} test",
        ]
    )
    return "\n".join(lines) + "\n"


def split_dataset(preprocess_manifest_path: Path, test_fraction: float, random_seed: int) -> dict:
    preprocess_manifest = read_json(preprocess_manifest_path)
    rng = random.Random(random_seed)

    split_manifest = {
        "random_seed": random_seed,
        "test_fraction": test_fraction,
        "genes": {},
    }
    train_records: dict[str, dict[str, list[dict]]] = {}
    test_records: dict[str, dict[str, list[dict]]] = {}

    ensure_dir(DATA_DIR / "train" / "healthy")
    ensure_dir(DATA_DIR / "train" / "unhealthy")
    ensure_dir(DATA_DIR / "test" / "healthy")
    ensure_dir(DATA_DIR / "test" / "unhealthy")

    for gene, payload in preprocess_manifest["genes"].items():
        train_records[gene] = {"HEALTHY": [], "UNHEALTHY": []}
        test_records[gene] = {"HEALTHY": [], "UNHEALTHY": []}
        split_manifest["genes"][gene] = {}

        for label_key, payload_key in (("HEALTHY", "healthy"), ("UNHEALTHY", "unhealthy")):
            records = []
            for entry in payload[payload_key]:
                fasta_path = resolve_local_fasta_path(entry, gene, payload_key)
                sequence_lines = fasta_path.read_text(encoding="utf-8").splitlines()[1:]
                records.append(
                    {
                        "id": entry["id"],
                        "sequence": "".join(sequence_lines),
                        "source_file": str(fasta_path.resolve()),
                        "label": label_key,
                    }
                )

            train_indices, test_indices = split_indices(len(records), test_fraction, rng)
            train_records[gene][label_key] = [records[index] for index in train_indices]
            test_records[gene][label_key] = [records[index] for index in test_indices]

        split_manifest["genes"][gene] = {
            "total_healthy": len(payload["healthy"]),
            "total_unhealthy": len(payload["unhealthy"]),
            "train_healthy": len(train_records[gene]["HEALTHY"]),
            "test_healthy": len(test_records[gene]["HEALTHY"]),
            "train_unhealthy": len(train_records[gene]["UNHEALTHY"]),
            "test_unhealthy": len(test_records[gene]["UNHEALTHY"]),
            "train_fraction": 1.0 - test_fraction,
            "test_fraction": test_fraction,
        }

    train_manifest = build_dataset_manifest("train", preprocess_manifest, train_records)
    test_manifest = build_dataset_manifest("test", preprocess_manifest, test_records)

    write_json(METADATA_DIR / "train_manifest.json", train_manifest)
    write_json(METADATA_DIR / "test_manifest.json", test_manifest)
    write_json(METADATA_DIR / "split_manifest.json", split_manifest)

    split_summary = create_split_summary(split_manifest)
    split_summary_path = METADATA_DIR / "split_summary.txt"
    split_summary_path.write_text(split_summary, encoding="utf-8")
    return split_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Split preprocessed AD gene FASTA files into reproducible train/test partitions.")
    parser.add_argument(
        "--preprocess-manifest",
        type=Path,
        default=METADATA_DIR / "preprocess_manifest.json",
        help="Path to preprocess manifest JSON generated by preprocess.py.",
    )
    parser.add_argument("--test-fraction", type=float, default=0.20, help="Fraction of each class assigned to the test set.")
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed used for reproducible splitting.")
    args = parser.parse_args()

    if not 0.0 < args.test_fraction < 1.0:
        raise ValueError("--test-fraction must be between 0 and 1.")

    split_manifest = split_dataset(args.preprocess_manifest.resolve(), args.test_fraction, args.random_seed)
    print(f"Split manifest written to {(METADATA_DIR / 'split_manifest.json').resolve()}")
    print(f"Train manifest written to {(METADATA_DIR / 'train_manifest.json').resolve()}")
    print(f"Test manifest written to {(METADATA_DIR / 'test_manifest.json').resolve()}")
    print(f"Summary written to {(METADATA_DIR / 'split_summary.txt').resolve()}")
    print(f"Genes split: {', '.join(sorted(split_manifest['genes']))}")


if __name__ == "__main__":
    main()
