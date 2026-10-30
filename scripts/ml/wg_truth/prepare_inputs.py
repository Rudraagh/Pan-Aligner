"""Extract named chromosomes from gzipped FASTAs and compare chromosomes between two references.

extract    write one unwrapped FASTA with .fai per chromosome, renamed if asked (name=newname)
compare    print whether a chromosome is identical in two FASTAs; this decides the projection route
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[2]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from simulate_truth_reads import write_fasta_with_fai  # noqa: E402


def read_records(path: Path, names: set[str]) -> dict[str, str]:
    """Upper-case sequences of the named records. Fails if any is missing."""
    opener = gzip.open if str(path).endswith(".gz") else open
    found: dict[str, list[str]] = {}
    current = None
    with opener(path, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                current = name if name in names else None
                if current is not None:
                    found[current] = []
            elif current is not None:
                found[current].append(line.strip().upper())
    missing = names - set(found)
    if missing:
        raise SystemExit(f"{path}: records not found: {sorted(missing)}")
    return {name: "".join(chunks) for name, chunks in found.items()}


def sequence_md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    extract = commands.add_parser("extract")
    extract.add_argument("--fasta", type=Path, required=True)
    extract.add_argument("--records", nargs="+", required=True, help="name or name=newname")
    extract.add_argument("--output-dir", type=Path, required=True)
    compare = commands.add_parser("compare")
    compare.add_argument("--first", type=Path, required=True)
    compare.add_argument("--second", type=Path, required=True)
    compare.add_argument("--record", required=True)
    compare.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.command == "extract":
        wanted = dict((item.split("=") + [item])[:2] for item in args.records)
        sequences = read_records(args.fasta, set(wanted))
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, new_name in wanted.items():
            path = args.output_dir / f"{new_name}.fa"
            if path.exists():
                raise SystemExit(f"{path} already exists; refusing to overwrite it")
            write_fasta_with_fai([(new_name, sequences[name])], path)
            print(new_name, len(sequences[name]), sequence_md5(sequences[name]), flush=True)
        return

    first = read_records(args.first, {args.record})[args.record]
    second = read_records(args.second, {args.record})[args.record]
    result = {"record": args.record, "first": str(args.first), "second": str(args.second), "first_length": len(first), "second_length": len(second),
              "first_md5": sequence_md5(first), "second_md5": sequence_md5(second), "identical": first == second}
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
