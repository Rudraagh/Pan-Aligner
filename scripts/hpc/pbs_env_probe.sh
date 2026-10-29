#!/bin/bash
# Read-only environment probe: which Python, compilers, modules and Conda envs a job sees.
# Needs no bundle; it answers "what is PYTHON_BIN / which module to load" without a shell.
# It installs nothing and makes no network connections.
#
# UNVERIFIED cluster settings (queue, memory units): confirm with the HPC administrator.
#PBS -N pa_env_probe
#PBS -q cpuq
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=00:10:00
#PBS -j oe
set -uo pipefail

WORK_DIR="${PBS_O_WORKDIR:-$(pwd)}"
OUT="$WORK_DIR/env_probe_$(date +%Y%m%d_%H%M%S)_${PBS_JOBID:-local}.txt"
exec > >(tee -a "$OUT") 2>&1

section() { echo; echo "===== $*"; }

section identity
echo "host: $(hostname)  user: $(id -un)  date: $(date -Is)"
echo "PBS_JOBID=${PBS_JOBID:-} PBS_QUEUE=${PBS_QUEUE:-} NCPUS=${NCPUS:-} TMPDIR=${TMPDIR:-}"
[[ -n "${PBS_JOBID:-}" ]] && command -v qstat >/dev/null 2>&1 && qstat -f "$PBS_JOBID" 2>&1 | grep -Ei 'resource|queue|walltime|mem|ncpus|exec_host'

section os
cat /etc/os-release 2>/dev/null | grep -E '^(NAME|VERSION)='
uname -srm
ldd --version 2>&1 | head -1
grep -m1 'model name' /proc/cpuinfo
echo "nproc visible: $(nproc)"

section storage
df -h "$HOME" "$WORK_DIR" 2>&1 | sort -u
for candidate in /scratch "/scratch/$(id -un)" /tmp "${TMPDIR:-}"; do
  [[ -n "$candidate" && -d "$candidate" ]] && echo "exists: $candidate"
done

section compilers
for c in gcc g++ cc make cmake git; do
  printf '%-6s %s\n' "$c" "$(command -v "$c" >/dev/null 2>&1 && "$c" --version 2>&1 | head -1 || echo 'not found')"
done
echo '#include <zlib.h>' | gcc -E - >/dev/null 2>&1 && echo "zlib headers: yes" || echo "zlib headers: no"

section python
for p in python3.13 python3.12 python3.11 python3.10 python3.9 python3 python; do
  command -v "$p" >/dev/null 2>&1 && echo "$p -> $(command -v "$p") $("$p" --version 2>&1)"
done

section modules
if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/modules.sh /usr/share/Modules/init/bash /etc/profile.d/lmod.sh; do
    [[ -f "$init" ]] && source "$init" && break
  done
fi
if type module >/dev/null 2>&1; then
  module list 2>&1
  module avail 2>&1 | head -200
else
  echo "module command not available"
fi

section conda
if command -v conda >/dev/null 2>&1; then
  conda --version 2>&1
  conda env list 2>&1
else
  echo "conda not on PATH"
fi

section gpu
echo "nvidia-smi: $(command -v nvidia-smi || echo 'not found')"
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi 2>&1 | head -20
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"

echo
echo "Probe written to $OUT"
