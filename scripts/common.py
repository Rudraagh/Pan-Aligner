from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys
from pathlib import Path, PureWindowsPath
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
GRAPHS_DIR = ROOT / "graphs"
OUTPUTS_DIR = ROOT / "outputs"
METADATA_DIR = DATA_DIR / "metadata"

DNA_PATTERN = re.compile(r"^[ACGTN]+$")
LABELED_HEADER_PATTERN = re.compile(r"^(?P<gene>[A-Za-z0-9]+)_(?P<label>[HU])_(?P<index>\d+)$")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def project_relative(path: Path) -> str:
    """Manifest form of a path: POSIX and relative to ROOT when inside the project."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return str(resolved)


def resolve_project_path(value: str | Path) -> Path:
    """Resolve a manifest path on Windows, WSL or Linux/HPC.

    Relative paths resolve against ROOT, not the working directory. Absolute paths
    written on another machine (e.g. C:\\Users\\...\\data\\raw\\x.fa) are remapped
    onto the first trailing project-relative portion that exists under ROOT.
    """
    text = str(value)
    windows_path = PureWindowsPath(text)
    candidate = Path(text)
    if not candidate.anchor and not windows_path.drive:
        return (ROOT / candidate).resolve()
    if candidate.exists():
        return candidate.resolve()
    parts = windows_path.parts if windows_path.drive else candidate.parts
    for index in range(1, len(parts)):
        remapped = ROOT.joinpath(*parts[index:])
        if remapped.exists():
            return remapped.resolve()
    return candidate


def normalize_sequence(sequence: str) -> str:
    return "".join(sequence.upper().split())


def validate_dna(sequence: str) -> bool:
    return bool(sequence) and bool(DNA_PATTERN.fullmatch(sequence))


def parse_sequence_header(header: str) -> dict:
    cleaned = header.strip()
    if cleaned.startswith("ref|"):
        return {
            "gene": infer_gene_from_reference_header(cleaned),
            "label": "REFERENCE",
            "is_reference": True,
            "sequence_id": cleaned,
        }

    match = LABELED_HEADER_PATTERN.match(cleaned)
    if not match:
        raise ValueError(
            f"Unsupported FASTA header '{cleaned}'. Expected 'GENE_H_#', 'GENE_U_#', or a reference header starting with 'ref|'."
        )

    label = "HEALTHY" if match.group("label") == "H" else "UNHEALTHY"
    return {
        "gene": match.group("gene").upper(),
        "label": label,
        "is_reference": False,
        "sequence_id": cleaned,
    }


def infer_gene_from_reference_header(header: str) -> str:
    normalized = header.upper()
    if "NC_000021" in normalized:
        return "APP"
    if "NC_000014" in normalized:
        return "PSEN1"
    if "NC_000001" in normalized:
        return "PSEN2"
    return "UNKNOWN"


def run_command(command: Sequence[str], cwd: Path | None = None, stdout_path: Path | None = None) -> subprocess.CompletedProcess:
    command_to_run = list(command)
    if command_to_run:
        launcher = Path(command_to_run[0])
        if launcher.name == "minigraph_wsl.cmd":
            command_to_run = [sys.executable, str(ROOT / "tools" / "wsl_exec.py"), "@minigraph", *command_to_run[1:]]
        elif launcher.name == "panaligner_wsl.cmd":
            command_to_run = [sys.executable, str(ROOT / "tools" / "wsl_exec.py"), "@panaligner", *command_to_run[1:]]

    if stdout_path is not None:
        ensure_dir(stdout_path.parent)
        with stdout_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                command_to_run,
                cwd=str(cwd) if cwd else None,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
    else:
        completed = subprocess.run(
            command_to_run,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if completed.returncode != 0:
        rendered = " ".join(shlex.quote(part) for part in command_to_run)
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}: {rendered}\nSTDERR:\n{completed.stderr}"
        )
    return completed


def write_json(path: Path, payload: dict | list) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_wrapped_fasta(records: Iterable[tuple[str, str]], output_path: Path, line_length: int = 80) -> None:
    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), line_length):
                handle.write(sequence[start : start + line_length] + "\n")


def resolve_binary(explicit: str | None, candidates: Sequence[Path]) -> Path:
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Binary not found: {path}")
        return path

    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.exists():
            return resolved

    candidate_text = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        f"Unable to locate required binary. Checked: {candidate_text}. Pass the path explicitly with the script argument."
    )


LEGACY_ML_NOTICE = (
    "LEGACY healthy/unhealthy ML path, kept for history only: it is not a valid classifier evaluation, "
    "so do not report its metrics. Training rows are aligned against class graphs built from those same "
    "sequences, so the healthy/unhealthy overlap features leak the label, while test rows (whose variant "
    "is in neither graph) have a different feature distribution. Each labelled sequence is the reference "
    "plus one variant, the labels correlate with genomic position, and the random per-sequence split puts "
    "neighbouring variants on both sides. See docs/13_current_version_implementation_status.md."
)


def warn_legacy_ml(module_name: str) -> None:
    print(f"[{module_name}] {LEGACY_ML_NOTICE}", file=sys.stderr)
