#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev python3 python3-pip

if [[ ! -d "$ROOT_DIR/PanAligner" ]]; then
  git clone https://github.com/at-cg/PanAligner "$ROOT_DIR/PanAligner"
fi

if [[ ! -d "$ROOT_DIR/minigraph" ]]; then
  git clone https://github.com/lh3/minigraph.git "$ROOT_DIR/minigraph"
fi

(cd "$ROOT_DIR/PanAligner" && make -j"$(nproc)")
(cd "$ROOT_DIR/minigraph" && make -j"$(nproc)")

python3 -m pip install -r "$ROOT_DIR/requirements.txt"

echo "PanAligner binary: $ROOT_DIR/PanAligner/PanAligner"
echo "minigraph binary: $ROOT_DIR/minigraph/minigraph"

