#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
THREADS="${THREADS:-4}"
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

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <query_fasta> [gene]"
  echo "Example: $0 data/queries/query.fa APP"
  exit 1
fi

QUERY_FASTA="$1"
GENE_ARG="${2:-}"

if [[ ! -f "$QUERY_FASTA" ]]; then
  echo "Query FASTA not found: $QUERY_FASTA"
  echo "Use an existing FASTA file, for example:"
  echo "  $0 data/unhealthy/app/APP_U_1.fa APP"
  echo "  $0 data/healthy/psen1/PSEN1_H_1.fa PSEN1"
  exit 1
fi

if [[ -n "$GENE_ARG" ]]; then
  GENE_FLAGS=(--gene "$GENE_ARG")
else
  GENE_FLAGS=()
fi

"$PYTHON_BIN" -m pip install -r "$ROOT_DIR/requirements.txt"

"$PYTHON_BIN" "$ROOT_DIR/scripts/pipeline.py" \
  --input-fastas "$ROOT_DIR/app_combined.fasta" "$ROOT_DIR/psen1_combined.fasta" "$ROOT_DIR/psen2_combined.fasta" \
  --query-fasta "$QUERY_FASTA" \
  --threads "$THREADS" \
  --panaligner-bin "$PANALIGNER_BIN" \
  --minigraph-bin "$MINIGRAPH_BIN" \
  "${GENE_FLAGS[@]}"
