#!/bin/bash
# Tiny CPU smoke test for the PBS cluster: environment report, checksums and the
# focused unit tests. No graph building, alignment or benchmarks.
#
# UNVERIFIED cluster settings: queue access, memory units (the Altair form's
# "Amount of Memory (MB)" field has shown inconsistent values) and the account.
# Confirm them with the HPC administrator; values set in the Altair form may
# override these directives.
#PBS -N pa_cpu_smoke
#PBS -q cpuq
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=00:15:00
#PBS -j oe
#
# Configure via environment (or edit here):
#   PROJECT_DIR  staged repository on the cluster (default: PBS submit directory)
#   MODULES      space-separated modules to load first, e.g. "gcc11" (names from pbs_env_probe.sh)
#   PYTHON_BIN   Python >= 3.10 (default: .hpc-venv from pbs_setup_venv.sh, else the first
#                python3.13..python3.10/python3 on PATH that qualifies)
#   RUN_BUILD=1  also run scripts/hpc/build_tools_linux.sh inside this job
set -uo pipefail

PROJECT_DIR="${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}"
MODULES="${MODULES:-}"
RUN_BUILD="${RUN_BUILD:-0}"
VENV_PYTHON="$PROJECT_DIR/.hpc-venv/bin/python"

# A venv built from a module Python needs that module loaded again to run.
if [[ -z "${PYTHON_BIN:-}" && -x "$VENV_PYTHON" ]]; then
  PYTHON_BIN="$VENV_PYTHON"
  MODULES="$(cat "$PROJECT_DIR/.hpc-venv/MODULES_USED" 2>/dev/null) $MODULES"
  MODULES="${MODULES# }"
  MODULES="${MODULES% }"
fi

if [[ -n "$MODULES" ]]; then
  if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
      [[ -f "$init" ]] && source "$init" && break
    done
  fi
  for mod in $MODULES; do
    module load "$mod" || { echo "ERROR: module load $mod failed" >&2; exit 2; }
  done
fi

if [[ -z "${PYTHON_BIN:-}" ]]; then
  for candidate in python3.13 python3.12 python3.11 python3.10 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; then
      PYTHON_BIN="$candidate"
      break
    fi
  done
  PYTHON_BIN="${PYTHON_BIN:-python3}"
fi
JOB_TAG="${PBS_JOBID:-local}"
RUN_DIR="$PROJECT_DIR/hpc_runs/$(date +%Y%m%d_%H%M%S)_cpu_smoke_${JOB_TAG%%.*}"

cd "$PROJECT_DIR" || { echo "ERROR: cannot cd to PROJECT_DIR=$PROJECT_DIR" >&2; exit 2; }
mkdir -p "$RUN_DIR"
exec > >(tee -a "$RUN_DIR/job.log") 2>&1

status=0
step() {
  local name="$1"; shift
  echo "=== $name: $*"
  "$@" > "$RUN_DIR/$name.log" 2>&1
  local code=$?
  cat "$RUN_DIR/$name.log"
  echo "=== $name exit=$code"
  echo "$name $code" >> "$RUN_DIR/step_status.txt"
  [[ $code -eq 0 ]] || status=1
}

# Commit identity: a git checkout, or the STAGED_COMMIT file filled in by `git archive`.
if command -v git >/dev/null 2>&1 && git -C "$PROJECT_DIR" rev-parse HEAD >/dev/null 2>&1; then
  COMMIT="$(git -C "$PROJECT_DIR" rev-parse HEAD)"
else
  COMMIT="$(head -1 "$PROJECT_DIR/scripts/hpc/STAGED_COMMIT" 2>/dev/null || echo unknown)"
fi

{
  echo "commit: $COMMIT"
  echo "date: $(date -Is)"
  echo "host: $(hostname)"
  echo "user: $(id -un)"
  echo "uname: $(uname -srm)"
  echo "glibc: $(ldd --version 2>&1 | head -1)"
  echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
  echo "nproc (visible): $(nproc)"
  echo "PBS_JOBID=${PBS_JOBID:-} PBS_QUEUE=${PBS_QUEUE:-} NCPUS=${NCPUS:-} PBS_O_WORKDIR=${PBS_O_WORKDIR:-}"
  echo "modules requested: ${MODULES:-none}"
  echo "python: $(command -v "$PYTHON_BIN") $("$PYTHON_BIN" --version 2>&1)"
  echo "gcc: $(gcc --version 2>&1 | head -1)"
  echo "nvidia-smi: $(command -v nvidia-smi || echo 'not found')"
} > "$RUN_DIR/environment.txt"
cat "$RUN_DIR/environment.txt"
if [[ -n "${PBS_JOBID:-}" ]] && command -v qstat >/dev/null 2>&1; then
  qstat -f "$PBS_JOBID" > "$RUN_DIR/qstat.txt" 2>&1 || true
fi
(module list) > "$RUN_DIR/modules.txt" 2>&1 || true

if ! "$PYTHON_BIN" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
  echo "ERROR: $PYTHON_BIN is $("$PYTHON_BIN" --version 2>&1); Python >= 3.10 is required." | tee "$RUN_DIR/FAILED"
  exit 3
fi
"$PYTHON_BIN" -m pip freeze > "$RUN_DIR/pip_freeze.txt" 2>&1 || true

(sha256sum data/raw/*.fa data/metadata/*.json *_combined.fasta > "$RUN_DIR/checksums.sha256") 2>>"$RUN_DIR/job.log"

export PYTHONDONTWRITEBYTECODE=1
step imports "$PYTHON_BIN" -c 'import Bio, matplotlib, networkx; print("biopython", Bio.__version__, "matplotlib", matplotlib.__version__, "networkx", networkx.__version__)'
step unit_tests "$PYTHON_BIN" -B -m unittest tests/test_proposed_minimizer_engine.py tests/test_path_portability.py -v

if [[ "$RUN_BUILD" == "1" ]]; then
  step build_tools bash scripts/hpc/build_tools_linux.sh
fi

"$PYTHON_BIN" - "$RUN_DIR" "$COMMIT" "$status" <<'EOF'
import json, os, sys
from pathlib import Path
run_dir, commit, status = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
steps = dict(line.split() for line in (run_dir / "step_status.txt").read_text().splitlines() if line.strip())
manifest = {
    "experiment": "cpu_smoke",
    "commit": commit,
    "pbs_jobid": os.environ.get("PBS_JOBID"),
    "queue": os.environ.get("PBS_QUEUE"),
    "ncpus": os.environ.get("NCPUS"),
    "command": "scripts/hpc/pbs_cpu_smoke.sh",
    "steps": {name: int(code) for name, code in steps.items()},
    "overall_exit": status,
    "gpu_validated": False,
    "files": sorted(p.name for p in run_dir.iterdir()),
}
(run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(json.dumps(manifest, indent=2))
EOF

echo "Results: $RUN_DIR (copy this folder back to the laptop)"
exit "$status"
