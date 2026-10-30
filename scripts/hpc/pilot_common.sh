# Sourced by the pbs_pilot_*.sh jobs. Sets up the environment and a run folder with a manifest.
#
# Configure via environment (or edit here):
#   PROJECT_DIR  extracted bundle folder (default: PBS submit directory)
#   PILOT_DIR    pilot data and results (default: $PROJECT_DIR/pilot_data); put this on scratch
#   MODULES      space-separated modules to load, e.g. "gcc11" (names from pbs_env_probe.sh)
#   PYTHON_BIN   Python >= 3.10 (default: .hpc-venv from pbs_setup_venv.sh)

PROJECT_DIR="${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}"
PILOT_DIR="${PILOT_DIR:-$PROJECT_DIR/pilot_data}"
MODULES="${MODULES:-}"
BIN_DIR="${BIN_DIR:-$PROJECT_DIR/.hpc-build/bin}"
MODEL="${MODEL:-$PROJECT_DIR/.hpc-build/pbsim2_data/P6C4.model}"
WG="$PROJECT_DIR/scripts/ml/wg_truth"
THREADS="${NCPUS:-1}"
export PYTHONDONTWRITEBYTECODE=1

if [[ -z "${PYTHON_BIN:-}" && -x "$PROJECT_DIR/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$PROJECT_DIR/.hpc-venv/bin/python"
  MODULES="$(cat "$PROJECT_DIR/.hpc-venv/MODULES_USED" 2>/dev/null) $MODULES"
fi
if [[ -n "${MODULES// /}" ]]; then
  if ! type module >/dev/null 2>&1; then
    for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
      [[ -f "$init" ]] && source "$init" && break
    done
  fi
  for mod in $MODULES; do
    module load "$mod" || { echo "ERROR: module load $mod failed" >&2; exit 2; }
  done
fi
PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! "$PYTHON_BIN" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)'; then
  echo "ERROR: $PYTHON_BIN is not Python >= 3.10." >&2
  exit 3
fi

# pilot_start NAME: create the run folder, send all output to its job.log, record the environment.
pilot_start() {
  JOB_TAG="${PBS_JOBID:-local}"
  RUN_DIR="$PILOT_DIR/runs/$(date +%Y%m%d_%H%M%S)_$1_${JOB_TAG%%.*}"
  STAGE="$1"
  cd "$PROJECT_DIR" || exit 2
  mkdir -p "$RUN_DIR"
  exec > >(tee -a "$RUN_DIR/job.log") 2>&1
  COMMIT="$(head -1 "$PROJECT_DIR/scripts/hpc/STAGED_COMMIT" 2>/dev/null || echo unknown)"
  {
    echo "stage: $STAGE"
    echo "commit: $COMMIT"
    echo "date: $(date -Is)"
    echo "host: $(hostname)"
    echo "user: $(id -un)"
    echo "cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
    echo "memory total: $(grep -m1 MemTotal /proc/meminfo)"
    echo "PBS_JOBID=${PBS_JOBID:-} PBS_QUEUE=${PBS_QUEUE:-} NCPUS=${NCPUS:-}"
    echo "modules: ${MODULES:-none}"
    echo "python: $("$PYTHON_BIN" --version 2>&1)"
    for tool in PanAligner minigraph; do
      echo "$tool: $("$BIN_DIR/$tool" --version 2>&1 | head -1)"
    done
    echo "minimap2: $("$BIN_DIR/minimap2" --version 2>&1 | head -1)"
    echo "gpu_used: no; this pilot is CPU only"
  } > "$RUN_DIR/environment.txt"
  cat "$RUN_DIR/environment.txt"
  if [[ -n "${PBS_JOBID:-}" ]] && command -v qstat >/dev/null 2>&1; then
    qstat -f "$PBS_JOBID" > "$RUN_DIR/qstat.txt" 2>&1 || true
  fi
}

# timed LABEL COMMAND...: run with wall time and, when available, peak memory recorded in RUN_DIR/timings.tsv.
timed() {
  local label="$1"; shift
  local start end code rss=""
  start="$(date +%s)"
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "$RUN_DIR/$label.time" "$@"
    code=$?
    rss="$(grep 'Maximum resident' "$RUN_DIR/$label.time" 2>/dev/null | awk '{print $NF}')"
  else
    "$@"
    code=$?
  fi
  end="$(date +%s)"
  printf '%s\t%s\t%s\t%s\n' "$label" "$code" "$((end - start))" "${rss:-unknown}" >> "$RUN_DIR/timings.tsv"
  echo "=== $label exit=$code seconds=$((end - start)) peak_kb=${rss:-unknown}"
  return "$code"
}

pilot_finish() {
  "$PYTHON_BIN" - "$RUN_DIR" "$STAGE" "$COMMIT" "$1" <<'EOF'
import hashlib, json, os, sys
from pathlib import Path
run_dir, stage, commit, status = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
timings = []
path = run_dir / "timings.tsv"
if path.exists():
    for line in path.read_text().splitlines():
        label, code, seconds, peak = line.split("\t")
        timings.append({"label": label, "exit": int(code), "seconds": int(seconds), "peak_kb": peak})
manifest = {"stage": stage, "commit": commit, "pbs_jobid": os.environ.get("PBS_JOBID"), "queue": os.environ.get("PBS_QUEUE"),
            "ncpus": os.environ.get("NCPUS"), "overall_exit": status, "gpu_validated": False, "timings": timings,
            "files": sorted(p.name for p in run_dir.iterdir())}
(run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
EOF
  echo "Results: $RUN_DIR (copy this folder back to the laptop)"
  exit "$1"
}
