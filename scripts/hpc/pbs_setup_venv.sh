#!/bin/bash
# Offline Python environment for the staged project: create .hpc-venv/ from the bundled
# wheelhouse/ (no network). pbs_cpu_smoke.sh then uses it automatically.
#
# Submit from inside the extracted panaligner_<commit>/ folder, after the unpack job.
# An existing working .hpc-venv is reused; a broken one is renamed, never deleted.
#
# UNVERIFIED cluster settings (queue, memory units): confirm with the HPC administrator.
#PBS -N pa_setup_venv
#PBS -q cpuq
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=00:20:00
#PBS -j oe
#
# Configure via environment (or edit here):
#   PROJECT_DIR  extracted bundle folder (default: PBS submit directory)
#   MODULES      space-separated modules that provide Python >= 3.10 (from pbs_env_probe.sh)
#   PYTHON_BIN   base Python >= 3.10 (default: first python3.13..python3.10/python3 that qualifies)
set -uo pipefail

PROJECT_DIR="${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}"
MODULES="${MODULES:-}"
JOB_TAG="${PBS_JOBID:-local}"
VENV="$PROJECT_DIR/.hpc-venv"
IMPORT_CHECK='import Bio, matplotlib, networkx; print("biopython", Bio.__version__, "matplotlib", matplotlib.__version__, "networkx", networkx.__version__)'

cd "$PROJECT_DIR" || { echo "ERROR: cannot cd to PROJECT_DIR=$PROJECT_DIR" >&2; exit 2; }
RUN_DIR="$PROJECT_DIR/hpc_runs/$(date +%Y%m%d_%H%M%S)_setup_venv_${JOB_TAG%%.*}"
mkdir -p "$RUN_DIR"
exec > >(tee -a "$RUN_DIR/setup.log") 2>&1
echo "host: $(hostname)  date: $(date -Is)  job: $JOB_TAG  project: $PROJECT_DIR"

fail() { echo "ERROR: $*" | tee "$RUN_DIR/FAILED"; exit "${2:-1}"; }

if [[ -x "$VENV/bin/python" ]] && "$VENV/bin/python" -c "$IMPORT_CHECK"; then
  echo "Existing $VENV already works; nothing to do."
  cp "$VENV/ENV_READY.txt" "$RUN_DIR/" 2>/dev/null || true
  exit 0
fi
if [[ -e "$VENV" ]]; then
  mv "$VENV" "$VENV.broken_$(date +%Y%m%d_%H%M%S)"
  echo "Renamed unusable $VENV"
fi

if [[ -n "$MODULES" ]]; then
  if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
      [[ -f "$init" ]] && source "$init" && break
    done
  fi
  for mod in $MODULES; do
    module load "$mod" || fail "module load $mod failed" 2
  done
fi

if [[ -z "${PYTHON_BIN:-}" ]]; then
  for candidate in python3.13 python3.12 python3.11 python3.10 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; then
      PYTHON_BIN="$candidate"
      break
    fi
  done
fi
[[ -n "${PYTHON_BIN:-}" ]] || fail "no Python >= 3.10 found on PATH; set MODULES or PYTHON_BIN (see pbs_env_probe.sh output)" 3
"$PYTHON_BIN" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' || fail "$PYTHON_BIN is $("$PYTHON_BIN" --version 2>&1); need >= 3.10" 3

tag="$("$PYTHON_BIN" -c 'import sys; print("py%d%d" % sys.version_info[:2])')"
WHEELS="$PROJECT_DIR/wheelhouse/$tag"
echo "base python: $(command -v "$PYTHON_BIN") $("$PYTHON_BIN" --version 2>&1)  wheels: $WHEELS"
[[ -d "$WHEELS" ]] || fail "no wheelhouse for $tag (have: $(ls "$PROJECT_DIR/wheelhouse" 2>/dev/null | tr '\n' ' ')); choose a matching Python via MODULES/PYTHON_BIN" 4

"$PYTHON_BIN" -m venv "$VENV" || fail "python -m venv failed (venv/ensurepip missing in this Python?)" 5
"$VENV/bin/python" -m pip install --no-index --find-links "$WHEELS" -r "$PROJECT_DIR/requirements.txt" || fail "offline pip install failed" 6
"$VENV/bin/python" -c "$IMPORT_CHECK" || fail "imports failed after install" 7

echo "$MODULES" > "$VENV/MODULES_USED"
{
  echo "python: $VENV/bin/python ($("$VENV/bin/python" --version 2>&1))"
  echo "base: $(command -v "$PYTHON_BIN")"
  echo "modules: ${MODULES:-none}"
  echo "wheelhouse: $WHEELS"
  "$VENV/bin/python" -c "$IMPORT_CHECK"
} > "$VENV/ENV_READY.txt"
"$VENV/bin/python" -m pip freeze > "$RUN_DIR/pip_freeze.txt"
cp "$VENV/ENV_READY.txt" "$RUN_DIR/"
cat "$VENV/ENV_READY.txt"
echo "OK: now submit scripts/hpc/pbs_cpu_smoke.sh from $PROJECT_DIR (it uses .hpc-venv automatically)"
