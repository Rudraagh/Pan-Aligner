"""Write or check scripts/hpc/pilot_manifest.sha256, the files a pilot checkout must contain.

Hashes come from the Git index, not the working tree, so a Windows checkout gives the same values as Linux.
Stage the files first. The commit hash itself is checked separately by pilot_preflight.sh.

    python scripts/hpc/make_pilot_manifest.py          write the manifest
    python scripts/hpc/make_pilot_manifest.py --check  exit 1 if it is out of date
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = "scripts/hpc/pilot_manifest.sha256"
REQUIRED = (
    ".gitattributes", ".gitignore", ".gitmodules", "requirements.txt",
    "docs/16_hpc_altair_workflow.md", "docs/18_hpc_git_deployment.md",
    "scripts/hpc/*", "scripts/ml/wg_truth/*", "scripts/ml/__init__.py", "scripts/ml/build_alignment_dataset.py",
    "scripts/common.py", "scripts/simulate_truth_reads.py", "scripts/truth_eval.py",
    "tests/test_pilot_git_deploy.py", "tests/test_pilot_hpc_scripts.py", "tests/test_wg_truth.py", "tests/test_wg_truth_report.py",
    "tools/third_party/k8-bin/k8-1.2/k8-x86_64-Linux", "wheelhouse/*",
)


def git(*args: str) -> bytes:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout


def build() -> str:
    lines = []
    for entry in git("ls-files", "-s", "-z").decode().split("\0"):
        if not entry:
            continue
        info, path = entry.split("\t", 1)
        if info.startswith("160000") or path == MANIFEST:
            continue
        if any(fnmatch.fnmatch(path, pattern) for pattern in REQUIRED):
            lines.append(f"{hashlib.sha256(git('cat-file', 'blob', info.split()[1])).hexdigest()}  {path}\n")
    return "".join(sorted(lines, key=lambda line: line.split("  ", 1)[1]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    content = build()
    path = ROOT / MANIFEST
    if args.check:
        if git("show", f":{MANIFEST}").decode() != content:
            sys.exit(f"{MANIFEST} is out of date")
        print(f"{MANIFEST}: {content.count(chr(10))} files match the index")
        return
    path.write_bytes(content.encode())
    print(f"{MANIFEST}: {content.count(chr(10))} files")


if __name__ == "__main__":
    main()
