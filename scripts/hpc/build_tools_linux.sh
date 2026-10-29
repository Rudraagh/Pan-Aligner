#!/bin/bash
# Build PanAligner and minigraph for the HPC node's own glibc, without sudo or network.
#
# The laptop binaries were built on Ubuntu 24.04 (need GLIBC 2.34/2.38) and are unlikely
# to run on the cluster. PanAligner/, minigraph/ and zlib/ are gitignored, so their
# sources must be staged next to the repository before running this.
#
# Run it inside a PBS job on the queue that will run the experiments: PanAligner's
# Makefile uses -march=native, so a binary built on another CPU model may crash.
#
# Usage: bash scripts/hpc/build_tools_linux.sh [BUILD_ROOT]
#   Load a C++20 compiler first if needed, e.g. `module load gcc11` (name unverified).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD_ROOT="${1:-${BUILD_ROOT:-$ROOT_DIR/.hpc-build}}"
JOBS="${JOBS:-${NCPUS:-1}}"
SRC_DIR="$BUILD_ROOT/src"
BIN_DIR="$BUILD_ROOT/bin"
ZLIB_PREFIX="$BUILD_ROOT/zlib"
INFO="$BUILD_ROOT/build_info.txt"

for tool in PanAligner minigraph; do
  if [[ ! -f "$ROOT_DIR/$tool/Makefile" ]]; then
    echo "ERROR: $ROOT_DIR/$tool source not found. It is gitignored; stage it with the repository." >&2
    exit 1
  fi
done

mkdir -p "$SRC_DIR" "$BIN_DIR"
{
  echo "date: $(date -Is)"
  echo "host: $(hostname)"
  echo "uname: $(uname -srm)"
  echo "glibc: $(ldd --version 2>&1 | head -1)"
  echo "gcc: $(gcc --version 2>&1 | head -1)"
  echo "g++: $(g++ --version 2>&1 | head -1)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
  echo "jobs: $JOBS"
} > "$INFO"

# Fresh copies keep the staged sources (and any laptop .o files) untouched.
for tool in PanAligner minigraph; do
  rm -rf "${SRC_DIR:?}/$tool"
  cp -a "$ROOT_DIR/$tool" "$SRC_DIR/$tool"
  (cd "$SRC_DIR/$tool" && make clean >/dev/null 2>&1 || true)
  if [[ -d "$ROOT_DIR/$tool/.git" ]] && command -v git >/dev/null 2>&1; then
    echo "$tool source: $(git -C "$ROOT_DIR/$tool" rev-parse HEAD 2>/dev/null || echo unknown)" >> "$INFO"
  else
    echo "$tool source: no git metadata staged" >> "$INFO"
  fi
done

# Use the system zlib headers when present; otherwise build the bundled zlib statically.
ZLIB_INCLUDES=""
ZLIB_LIBDIR=""
if echo '#include <zlib.h>' | gcc -E - >/dev/null 2>&1; then
  echo "zlib: system" >> "$INFO"
else
  if [[ ! -f "$ROOT_DIR/zlib/configure" ]]; then
    echo "ERROR: zlib headers missing and bundled zlib/ source not staged." >&2
    exit 1
  fi
  rm -rf "${SRC_DIR:?}/zlib"
  cp -a "$ROOT_DIR/zlib" "$SRC_DIR/zlib"
  (cd "$SRC_DIR/zlib" && ./configure --static --prefix="$ZLIB_PREFIX" >/dev/null && make -j"$JOBS" >/dev/null && make install >/dev/null)
  ZLIB_INCLUDES="-I$ZLIB_PREFIX/include"
  ZLIB_LIBDIR="-L$ZLIB_PREFIX/lib"
  echo "zlib: bundled static ($ZLIB_PREFIX)" >> "$INFO"
fi

(cd "$SRC_DIR/minigraph" && make -j"$JOBS" INCLUDES="$ZLIB_INCLUDES" LIBS="$ZLIB_LIBDIR -lz -lpthread -lm")
(cd "$SRC_DIR/PanAligner" && make -j"$JOBS" INCLUDES="$ZLIB_INCLUDES" LIBS="$ZLIB_LIBDIR -lz -lpthread -lm -fopenmp")

cp "$SRC_DIR/minigraph/minigraph" "$SRC_DIR/PanAligner/PanAligner" "$BIN_DIR/"
{
  echo "minigraph --version: $("$BIN_DIR/minigraph" --version 2>&1 | head -1)"
  echo "PanAligner --version: $("$BIN_DIR/PanAligner" --version 2>&1 | head -1)"
  echo "minigraph max GLIBC: $(objdump -T "$BIN_DIR/minigraph" 2>/dev/null | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -1)"
  echo "PanAligner max GLIBC: $(objdump -T "$BIN_DIR/PanAligner" 2>/dev/null | grep -o 'GLIBC_[0-9.]*' | sort -uV | tail -1)"
} >> "$INFO"

cat "$INFO"
echo "Binaries: $BIN_DIR (pass with --minigraph-bin/--panaligner-bin)"
