from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

from common import METADATA_DIR, OUTPUTS_DIR, ensure_dir, write_json


ROOT = Path(__file__).resolve().parents[1]
PAFTOOLS_WRAPPER = ROOT / "tools" / "paftools_wsl.cmd"
K8_BINARY = ROOT / "tools" / "third_party" / "k8-bin" / "k8-1.2" / "k8-x86_64-Linux"
PAFTOOLS_JS = ROOT / "tools" / "third_party" / "minimap2" / "misc" / "paftools.js"
COORDINATE_HEADER_PATTERN = re.compile(r".+[:_]\d+[-_]\d+")


def first_fasta_header(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].strip()
    return ""


def header_has_truth_coordinates(header: str) -> bool:
    return bool(COORDINATE_HEADER_PATTERN.search(header))


def run_paftools_help() -> tuple[bool, str]:
    completed = subprocess.run(
        [str(PAFTOOLS_WRAPPER)],
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = completed.stdout.strip()
    ok = "Usage: paftools.js" in output and "Commands:" in output
    return ok, output


def build_report(test_manifest_path: Path) -> dict:
    manifest = __import__("json").loads(test_manifest_path.read_text(encoding="utf-8"))
    sample_headers: list[dict[str, str | bool]] = []
    coordinate_ready = True

    for gene, payload in manifest["genes"].items():
        found_one = False
        for bucket in ("healthy", "unhealthy"):
            for entry in payload[bucket]:
                fasta_path = Path(entry["path"])
                header = first_fasta_header(fasta_path)
                has_coords = header_has_truth_coordinates(header)
                coordinate_ready = coordinate_ready and has_coords
                sample_headers.append(
                    {
                        "gene": gene,
                        "bucket": bucket,
                        "sequence_id": entry["id"],
                        "header": header,
                        "has_truth_coordinates": has_coords,
                    }
                )
                found_one = True
                break
            if found_one:
                break

    paftools_ok, paftools_output = run_paftools_help()
    exact_paper_eval_possible = paftools_ok and coordinate_ready

    return {
        "paftools_wrapper": str(PAFTOOLS_WRAPPER.resolve()),
        "k8_binary": str(K8_BINARY.resolve()),
        "paftools_js": str(PAFTOOLS_JS.resolve()),
        "paftools_runnable": paftools_ok,
        "coordinate_encoded_headers_detected_for_all_samples": coordinate_ready,
        "exact_paper_style_paftools_correctness_eval_possible": exact_paper_eval_possible,
        "why_not_possible": (
            ""
            if exact_paper_eval_possible
            else "Current dataset headers do not encode truth interval coordinates like the simulated CHM13 reads used in the paper, so paftools-based correctness evaluation cannot be reproduced exactly."
        ),
        "sample_headers": sample_headers,
        "paftools_help_snippet": paftools_output.splitlines()[:12],
    }


def write_report_files(report: dict, output_dir: Path) -> None:
    ensure_dir(output_dir)
    write_json(output_dir / "paftools_readiness.json", report)
    lines = [
        "PafTools Readiness Report",
        "=========================",
        f"Wrapper: {report['paftools_wrapper']}",
        f"k8 binary: {report['k8_binary']}",
        f"paftools.js: {report['paftools_js']}",
        f"PafTools runnable: {report['paftools_runnable']}",
        f"Coordinate-encoded headers available for all samples: {report['coordinate_encoded_headers_detected_for_all_samples']}",
        f"Exact paper-style paftools correctness eval possible: {report['exact_paper_style_paftools_correctness_eval_possible']}",
    ]
    if report["why_not_possible"]:
        lines.extend(["", "Limitation", report["why_not_possible"]])
    lines.extend(["", "Sample headers checked"])
    for sample in report["sample_headers"]:
        lines.append(
            f"- {sample['gene']} {sample['bucket']} {sample['sequence_id']}: "
            f"{sample['header']} | truth-coordinates={sample['has_truth_coordinates']}"
        )
    lines.extend(["", "PafTools help snippet"])
    lines.extend(f"- {line}" for line in report["paftools_help_snippet"])
    (output_dir / "paftools_readiness_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Check whether paftools can be used meaningfully with the current dataset.")
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "evaluation")
    args = parser.parse_args()

    report = build_report(args.test_manifest.resolve())
    write_report_files(report, args.output_dir.resolve())
    print(f"PafTools runnable={report['paftools_runnable']}")
    print(f"Exact paper-style eval possible={report['exact_paper_style_paftools_correctness_eval_possible']}")


if __name__ == "__main__":
    main()
