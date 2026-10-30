#!/bin/bash
# Pilot stage 2: load the whole graph and align about 1,000 reads per haplotype, measuring peak memory and speed.
# Its result sets the walltime for stage 4. UNVERIFIED cluster settings: confirm with the administrator.
#PBS -N pa_pilot_load
#PBS -q cpuq
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=01:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
pilot_start loadtest
PREPARED="$PILOT_DIR/prepared"
HAPS=(); for f in "$PREPARED"/haplotypes/*.fa; do name="$(basename "$f" .fa)"; case "$name" in *MAT*) HAPS+=("mat=$f");; *PAT*) HAPS+=("pat=$f");; esac; done
status=0
timed simulate "$PYTHON_BIN" "$WG/simulate_haplotypes.py" --haplotypes "${HAPS[@]}" --accuracies 0.95 --lengths 10000 --depth 0.22 --pbsim "$BIN_DIR/pbsim" --model "$MODEL" --output-dir "$RUN_DIR/sim" || pilot_finish 1
READS="$RUN_DIR/sim/acc95_len10000/all.reads.fa"
cat "$RUN_DIR"/sim/acc95_len10000/*.reads.fa > "$READS"
grep -c '^>' "$READS" | tee "$RUN_DIR/read_count.txt"
for aligner in minigraph PanAligner; do
  timed "align_$aligner" "$BIN_DIR/$aligner" -t "$THREADS" -cx lr "$PREPARED/graph.gfa" "$READS" > "$RUN_DIR/$aligner.gaf" 2> "$RUN_DIR/$aligner.log" || status=1
done
# Real path naming, needed to confirm the projection assumptions.
for aligner in minigraph PanAligner; do
  cut -f6 "$RUN_DIR/$aligner.gaf" | awk 'length($0) < 400' | head -200 > "$RUN_DIR/$aligner.path_examples.txt"
done
pilot_finish "$status"
