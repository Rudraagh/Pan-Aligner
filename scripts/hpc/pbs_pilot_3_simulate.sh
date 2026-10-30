#!/bin/bash
# Pilot stage 3: simulate all 6 settings at 1x per haplotype, and repeat the smallest setting to test reproducibility.
# UNVERIFIED cluster settings: confirm with the administrator.
#PBS -N pa_pilot_sim
#PBS -q cpuq
#PBS -l select=1:ncpus=4:mem=8gb
#PBS -l walltime=02:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
pilot_start simulate
PREPARED="$PILOT_DIR/prepared"
SIM="$PILOT_DIR/sim"
[[ -e "$SIM" ]] && { echo "ERROR: $SIM already exists; refusing to overwrite it" >&2; pilot_finish 1; }
HAPS=(); for f in "$PREPARED"/haplotypes/*.fa; do name="$(basename "$f" .fa)"; case "$name" in *MAT*) HAPS+=("mat=$f");; *PAT*) HAPS+=("pat=$f");; esac; done
[[ ${#HAPS[@]} -eq 2 ]] || { echo "ERROR: expected 2 haplotypes in $PREPARED/haplotypes" >&2; pilot_finish 1; }
timed simulate "$PYTHON_BIN" "$WG/simulate_haplotypes.py" --haplotypes "${HAPS[@]}" --depth 1 --pbsim "$BIN_DIR/pbsim" --model "$MODEL" --output-dir "$SIM" || pilot_finish 1
timed repeat "$PYTHON_BIN" "$WG/simulate_haplotypes.py" --haplotypes "${HAPS[@]}" --depth 1 --accuracies 0.95 --lengths 10000 --pbsim "$BIN_DIR/pbsim" --model "$MODEL" --output-dir "$RUN_DIR/repeat" || pilot_finish 1
identical=true
for hap in mat pat; do
  cmp -s "$SIM/acc95_len10000/$hap.reads.fa" "$RUN_DIR/repeat/acc95_len10000/$hap.reads.fa" || identical=false
done
echo "{\"identical\": $identical}" | tee "$PILOT_DIR/repro.json"
cp "$PILOT_DIR/repro.json" "$SIM/simulation_manifest.json" "$RUN_DIR/"
rm -rf "${RUN_DIR:?}/repeat"
pilot_finish 0
