#!/bin/bash
# Pilot stage 4: align every setting with one aligner. Set ALIGNER=minigraph or ALIGNER=PanAligner, one job each.
# Resumable: a setting with a finished GAF is skipped, so a job that hits its limit can be resubmitted.
# Set the walltime from the stage 2 measurement. UNVERIFIED cluster settings: confirm with the administrator.
#PBS -N pa_pilot_align
#PBS -q cpuq
#PBS -l select=1:ncpus=32:mem=64gb
#PBS -l walltime=08:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
ALIGNER="${ALIGNER:?set ALIGNER=minigraph or ALIGNER=PanAligner}"
pilot_start "align_$ALIGNER"
status=0
for dir in "$PILOT_DIR"/sim/acc*; do
  setting="$(basename "$dir")"
  gaf="$dir/$ALIGNER.gaf"
  [[ -s "$gaf" ]] && { echo "done already: $gaf"; continue; }
  [[ -s "$dir/all.reads.fa" ]] || cat "$dir"/mat.reads.fa "$dir"/pat.reads.fa > "$dir/all.reads.fa"
  if timed "${ALIGNER}_$setting" "$BIN_DIR/$ALIGNER" -t "$THREADS" -cx lr "$PILOT_DIR/prepared/graph.gfa" "$dir/all.reads.fa" > "$gaf.partial" 2> "$gaf.log"; then
    mv "$gaf.partial" "$gaf"
  else
    status=1
    echo "ERROR: $ALIGNER failed on $setting; see $gaf.log"
  fi
done
pilot_finish "$status"
