#!/bin/bash
# Verify and extract a staging bundle on the cluster (for Altair Access, where there may
# be no interactive shell). Upload this file next to the bundle and its .sha256 file.
# It never overwrites: if the destination folder exists, it stops.
#
# UNVERIFIED cluster settings (queue, memory units): confirm with the HPC administrator.
#PBS -N pa_unpack
#PBS -q cpuq
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:15:00
#PBS -j oe
#
# Configure via environment (or edit here):
#   BUNDLE     path to panaligner_<commit>.tar.gz (default: newest one in the submit dir)
#   DEST_ROOT  where to extract (default: the submit dir)
set -euo pipefail

WORK_DIR="${PBS_O_WORKDIR:-$(pwd)}"
cd "$WORK_DIR"
BUNDLE="${BUNDLE:-$(ls -t panaligner_*.tar.gz 2>/dev/null | head -1 || true)}"
DEST_ROOT="${DEST_ROOT:-$WORK_DIR}"
LOG="$WORK_DIR/unpack_$(date +%Y%m%d_%H%M%S)_${PBS_JOBID:-local}.log"
exec > >(tee -a "$LOG") 2>&1

echo "host: $(hostname)  date: $(date -Is)  job: ${PBS_JOBID:-local}"
if [[ -z "$BUNDLE" || ! -f "$BUNDLE" ]]; then
  echo "ERROR: no bundle found (set BUNDLE). Looked in $WORK_DIR." >&2
  exit 1
fi
bundle_dir="$(cd "$(dirname "$BUNDLE")" && pwd)"
bundle_file="$(basename "$BUNDLE")"
name="${bundle_file%.tar.gz}"
dest="$DEST_ROOT/$name"

echo "== verifying upload checksum"
if [[ ! -f "$bundle_dir/$bundle_file.sha256" ]]; then
  echo "ERROR: $bundle_file.sha256 missing; upload it with the bundle." >&2
  exit 1
fi
(cd "$bundle_dir" && sha256sum -c "$bundle_file.sha256")

if [[ -e "$dest" ]]; then
  echo "ERROR: $dest already exists; not overwriting. Remove or rename it yourself." >&2
  exit 1
fi
mkdir -p "$DEST_ROOT"
echo "== extracting to $dest"
tar -xzf "$bundle_dir/$bundle_file" -C "$DEST_ROOT"

echo "== verifying extracted files"
(cd "$dest" && sha256sum -c --quiet STAGED_FILES.sha256)
echo "files verified: $(wc -l < "$dest/STAGED_FILES.sha256")"
echo "commit: $(head -1 "$dest/scripts/hpc/STAGED_COMMIT")"
cat "$dest/STAGED_TOOLS.txt"
echo "OK: next, submit $dest/scripts/hpc/pbs_cpu_smoke.sh from inside $dest"
