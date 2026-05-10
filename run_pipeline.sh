#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-4}"
SPLIT_SEED="${SPLIT_SEED:-42}"
TEST_FRACTION="${TEST_FRACTION:-0.20}"
PANALIGNER_BIN="${PANALIGNER_BIN:-$ROOT_DIR/PanAligner/PanAligner}"
MINIGRAPH_BIN="${MINIGRAPH_BIN:-$ROOT_DIR/minigraph/minigraph}"

if [[ ! -d "$ROOT_DIR/PanAligner" ]]; then
  echo "PanAligner repository not found in $ROOT_DIR/PanAligner"
  exit 1
fi

if [[ ! -x "$PANALIGNER_BIN" ]]; then
  echo "PanAligner binary not found. Attempting local build..."
  (cd "$ROOT_DIR/PanAligner" && make -j"$(nproc)")
fi

if [[ ! -d "$ROOT_DIR/minigraph" ]]; then
  echo "minigraph repository not found. Clone it with:"
  echo "  git clone https://github.com/lh3/minigraph.git"
  exit 1
fi

if [[ ! -x "$MINIGRAPH_BIN" ]]; then
  echo "minigraph binary not found. Attempting local build..."
  (cd "$ROOT_DIR/minigraph" && make -j"$(nproc)")
fi

"$PYTHON_BIN" -m pip install -r "$ROOT_DIR/requirements.txt"

"$PYTHON_BIN" "$ROOT_DIR/main.py" \
  --full-pipeline \
  --input-fastas "$ROOT_DIR/app_combined.fasta" "$ROOT_DIR/psen1_combined.fasta" "$ROOT_DIR/psen2_combined.fasta" \
  --threads "$THREADS" \
  --split-seed "$SPLIT_SEED" \
  --test-fraction "$TEST_FRACTION" \
  --panaligner-bin "$PANALIGNER_BIN" \
  --minigraph-bin "$MINIGRAPH_BIN"
