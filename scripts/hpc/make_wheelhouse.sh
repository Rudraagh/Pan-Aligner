#!/bin/bash
# Download (not install) Linux wheels for requirements.txt so an offline cluster can build a
# venv with pbs_setup_venv.sh. Run on the laptop; needs internet, installs nothing.
#
# manylinux2014 (glibc >= 2.17) is targeted so wheels also work on older RHEL/CentOS nodes.
#
# Usage: bash scripts/hpc/make_wheelhouse.sh [OUT_DIR] [PYTHON_VERSIONS...]
#   default OUT_DIR hpc_staging/wheelhouse, default versions 3.10 3.11 3.12
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${1:-$ROOT_DIR/hpc_staging/wheelhouse}"
shift || true
VERSIONS=("$@")
[[ ${#VERSIONS[@]} -gt 0 ]] || VERSIONS=(3.10 3.11 3.12)
PIP_PYTHON="${PIP_PYTHON:-python}"

mkdir -p "$OUT_DIR"
for version in "${VERSIONS[@]}"; do
  tag="py${version/./}"
  echo "== $tag"
  "$PIP_PYTHON" -m pip download \
    --dest "$OUT_DIR/$tag" \
    --only-binary=:all: \
    --platform manylinux2014_x86_64 \
    --implementation cp \
    --python-version "$version" \
    -r "$ROOT_DIR/requirements.txt"
  (cd "$OUT_DIR/$tag" && ls -1 *.whl > WHEELS.txt)
done
"$PIP_PYTHON" -m pip --version > "$OUT_DIR/DOWNLOADED_WITH.txt"
date -Is >> "$OUT_DIR/DOWNLOADED_WITH.txt"
du -sh "$OUT_DIR"/py*
