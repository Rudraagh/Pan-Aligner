#!/bin/bash
# Pilot stage 0d: build the tools on this node type, run the unit tests and the synthetic smoke test.
# No real genome data. UNVERIFIED cluster settings (queue, memory units): confirm with the administrator.
#PBS -N pa_pilot_smoke
#PBS -q cpuq
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=01:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
JOBS="$THREADS" bash "$PROJECT_DIR/scripts/hpc/build_tools_linux.sh" || { echo "ERROR: tool build failed" >&2; exit 4; }
pilot_start build_smoke
status=0
timed unit_tests "$PYTHON_BIN" -B -m unittest tests/test_wg_truth.py tests/test_wg_truth_report.py -v || status=1
timed synthetic bash "$WG/synthetic_smoke.sh" "$RUN_DIR/synthetic" "$PYTHON_BIN" "$BIN_DIR/PanAligner" "$BIN_DIR/minigraph" "$BIN_DIR/minimap2" "$BIN_DIR/pbsim" "$MODEL" || status=1
pilot_finish "$status"
