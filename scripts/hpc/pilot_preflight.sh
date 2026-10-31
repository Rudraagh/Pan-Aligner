#!/bin/bash
# Read-only checks before Job 0. It never submits a job, downloads data, or changes a shared environment.
# It writes only a temporary compile test under ${TMPDIR:-/tmp}, which it removes.
#
# Usage, from the checkout's top folder:
#   bash scripts/hpc/pilot_preflight.sh --expected-commit SHA [--before-venv] [--no-pbs] [--branch NAME]
#     --before-venv  skip the checks on .hpc-venv, for use before pbs_setup_venv.sh
#     --no-pbs       skip the PBS command check, for testing away from the cluster
# Exit: 0 all required checks passed, 1 at least one failed, 2 bad usage.
set -uo pipefail

EXPECTED=""
BRANCH="hpc/model-a-giab-pilot"
BEFORE_VENV=0
NO_PBS=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --expected-commit) EXPECTED="${2:-}"; shift 2 ;;
    --branch) BRANCH="${2:-}"; shift 2 ;;
    --before-venv) BEFORE_VENV=1; shift ;;
    --no-pbs) NO_PBS=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
if [[ -z "$EXPECTED" ]]; then
  echo "usage: $0 --expected-commit SHA [--before-venv] [--no-pbs] [--branch NAME]" >&2
  exit 2
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT" || exit 2
source scripts/hpc/hpc_env.sh
failures=0
warnings=0
pass() { echo "PASS  $*"; }
fail() { echo "FAIL  $*"; failures=$((failures + 1)); }
warn() { echo "WARN  $*"; warnings=$((warnings + 1)); }
info() { echo "INFO  $*"; }
section() { echo; echo "== $*"; }
need() { if command -v "$1" >/dev/null 2>&1; then pass "$1: $(command -v "$1")"; else fail "$1 is not on PATH"; fi; }

section "Git"
need git
head_commit="$(git rev-parse HEAD 2>/dev/null || echo none)"
if [[ "$head_commit" == "$EXPECTED" ]]; then pass "HEAD is the expected commit $EXPECTED"; else fail "HEAD is $head_commit, expected commit $EXPECTED"; fi
branch="$(git branch --show-current 2>/dev/null)"
if [[ "$branch" == "$BRANCH" ]]; then pass "branch $branch"; else fail "branch is '${branch:-detached}', expected $BRANCH"; fi
changed="$(git status --porcelain --untracked-files=no 2>/dev/null)"
if [[ -z "$changed" ]]; then pass "no tracked file differs from the commit"; else fail "tracked files differ from the commit:"; echo "$changed" | head -20; fi

section "Tool sources pinned by the commit"
for path in PanAligner minigraph zlib tools/third_party/pbsim2 tools/third_party/minimap2; do
  pinned="$(git ls-tree HEAD -- "$path" | awk '{print $3}')"
  actual="$(git -C "$path" rev-parse HEAD 2>/dev/null || echo none)"
  if [[ ! -e "$path/.git" ]]; then
    fail "$path is not checked out; run: git submodule update --init -- $path"
  elif [[ -n "$pinned" && "$pinned" == "$actual" ]]; then
    pass "$path at $pinned"
  else
    fail "$path is at $actual, pinned $pinned"
  fi
done

section "Required files and checksums"
if [[ -s scripts/hpc/pilot_manifest.sha256 ]]; then
  if output="$(sha256sum --quiet -c scripts/hpc/pilot_manifest.sha256 2>&1)"; then
    pass "$(wc -l < scripts/hpc/pilot_manifest.sha256) files match pilot_manifest.sha256"
  else
    fail "checksum mismatch or missing file:"; echo "$output" | head -20
  fi
else
  fail "scripts/hpc/pilot_manifest.sha256 is missing"
fi

section "Python"
if [[ -n "${PYTHON_BIN:-}" ]]; then info "PYTHON_BIN overrides the base interpreter: $PYTHON_BIN"; HPC_PYTHON="$PYTHON_BIN"; fi
if [[ -x "$(command -v "$HPC_PYTHON" 2>/dev/null)" ]]; then
  version="$("$HPC_PYTHON" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
  if [[ "$version" == "$HPC_PYTHON_VERSION" || "$version" == "$HPC_PYTHON_VERSION".* ]]; then pass "$HPC_PYTHON is Python $version"; else fail "$HPC_PYTHON is Python $version, expected $HPC_PYTHON_VERSION"; fi
  if "$HPC_PYTHON" -c 'import venv, ensurepip' 2>/dev/null; then pass "venv and ensurepip are available"; else fail "venv or ensurepip is missing in $HPC_PYTHON"; fi
  tag="$("$HPC_PYTHON" -c 'import sys; print("py%d%d" % sys.version_info[:2])')"
else
  fail "base interpreter not found: $HPC_PYTHON"
  tag="py${HPC_PYTHON_VERSION/./}"
fi
info "default python3 here: $(command -v python3 || echo none) $(python3 --version 2>&1). Not used."

section "Dependencies"
wheels="wheelhouse/$tag"
if [[ -s "$wheels/WHEELS.txt" ]]; then
  missing=0
  while read -r name; do [[ -f "$wheels/$name" ]] || { fail "missing wheel $wheels/$name"; missing=1; }; done < "$wheels/WHEELS.txt"
  [[ $missing -eq 0 ]] && pass "offline wheels present in $wheels"
else
  fail "no offline wheels for $tag in wheelhouse/"
fi
if [[ $BEFORE_VENV -eq 1 ]]; then
  info "--before-venv: .hpc-venv not checked"
elif [[ -x .hpc-venv/bin/python ]]; then
  if .hpc-venv/bin/python -c 'import numpy, Bio, matplotlib, networkx; print("numpy", numpy.__version__, "biopython", Bio.__version__, "matplotlib", matplotlib.__version__, "networkx", networkx.__version__)'; then
    pass ".hpc-venv imports the required packages"
  else
    fail ".hpc-venv cannot import the required packages"
  fi
else
  fail ".hpc-venv is missing; run scripts/hpc/pbs_setup_venv.sh"
fi

section "System tools"
for tool in gcc g++ make tar gzip gunzip awk sed sort comm cut sha256sum md5sum; do need "$tool"; done
for tool in objdump curl samtools; do
  if command -v "$tool" >/dev/null 2>&1; then pass "$tool: $(command -v "$tool")"; else warn "$tool not found: objdump is informational; curl is for a later job; samtools only for the giab route"; fi
done
if [[ -x /usr/bin/time ]]; then pass "/usr/bin/time"; else warn "/usr/bin/time is missing: peak memory would be reported as unknown"; fi
info "glibc: $(getconf GNU_LIBC_VERSION 2>/dev/null || echo unknown). Wheels need 2.17 or newer."
info "gcc: $(gcc --version 2>&1 | head -1)"
info "bash: $BASH_VERSION. pbs_pilot_5_label.sh needs 4.4 or newer if its k-mer check fails."

section "Compile check"
scratch="$(mktemp -d "${TMPDIR:-/tmp}/pilot_preflight.XXXXXX")"
trap 'rm -rf "${scratch:?}"' EXIT
printf 'int main() { auto f = []<typename T>(T x) { return x; }; return f(0); }\n' > "$scratch/cxx.cpp"
if g++ -std=c++2a "$scratch/cxx.cpp" -o "$scratch/cxx" 2> "$scratch/cxx.log"; then pass "g++ accepts -std=c++2a"; else fail "g++ cannot compile C++20; load a newer compiler"; head -5 "$scratch/cxx.log"; fi
printf '#include <zlib.h>\nint main(void) { return zlibVersion() == 0; }\n' > "$scratch/z.c"
if gcc "$scratch/z.c" -lz -o "$scratch/z" 2>/dev/null; then pass "zlib links"; else warn "zlib does not link; the build falls back to the zlib source"; fi

section "k8 and read simulation inputs"
K8=tools/third_party/k8-bin/k8-1.2/k8-x86_64-Linux
PAFTOOLS=tools/third_party/minimap2/misc/paftools.js
if [[ -x "$K8" ]]; then
  pass "$K8 is executable"
  if version="$("$K8" "$PAFTOOLS" version 2>&1)"; then pass "k8 starts; paftools $version"; else fail "k8 does not start: $version"; fi
else
  fail "$K8 is not executable or is missing"
fi
for file in tools/third_party/pbsim2/data/P6C4.model tools/third_party/pbsim2/configure scripts/ml/wg_truth/synthetic_smoke.sh scripts/ml/wg_truth/synthetic_fixture.py; do
  if [[ -s "$file" ]]; then pass "$file"; else fail "missing $file"; fi
done

section "Jobs"
for job in scripts/hpc/*.sh scripts/ml/wg_truth/synthetic_smoke.sh; do
  bash -n "$job" 2>/dev/null || fail "syntax error in $job"
done
pass "shell syntax checked"
while read -r job cpus mem wall; do
  file="scripts/hpc/$job"
  if grep -qx "#PBS -l select=1:ncpus=$cpus:mem=$mem" "$file" && grep -qx "#PBS -l walltime=$wall" "$file"; then pass "$job: $cpus CPUs, $mem, $wall"; else fail "$job does not request $cpus CPUs, $mem, $wall"; fi
  if grep -Eq 'ngpus|gpuq' "$file"; then fail "$job refers to a GPU"; fi
done <<'EOF'
pbs_pilot_0_build_smoke.sh 4 8gb 01:00:00
pbs_pilot_1_inputs.sh 16 32gb 04:00:00
pbs_pilot_2_loadtest.sh 32 64gb 01:00:00
pbs_pilot_3_simulate.sh 4 8gb 02:00:00
pbs_pilot_4_align.sh 32 64gb 08:00:00
pbs_pilot_5_label.sh 4 32gb 04:00:00
EOF
if [[ $NO_PBS -eq 1 ]]; then
  info "--no-pbs: PBS commands not checked"
else
  for tool in qsub qstat; do need "$tool"; done
fi

section "Storage"
info "free space for ${PILOT_DIR:-$ROOT}: $(df -h "${PILOT_DIR:-$ROOT}" 2>/dev/null | awk 'NR==2 {print $4}'). Job 0 needs little; later jobs are estimated at 25 GB."

echo
echo "Result: $failures failed, $warnings warnings. This is not cluster validation; only Job 0 tests the compute nodes."
[[ $failures -eq 0 ]]
