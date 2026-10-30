#!/bin/bash
# Pilot stage 5: label every GAF, run the k-mer independence check, and compute the go/no-go gates.
# The k-mer check scans the whole reference. UNVERIFIED cluster settings: confirm with the administrator.
#PBS -N pa_pilot_label
#PBS -q cpuq
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=04:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
pilot_start label
PREPARED="$PILOT_DIR/prepared"
INPUTS="$PILOT_DIR/inputs"
LABELS="$PILOT_DIR/labels"
[[ -e "$LABELS" ]] && { echo "ERROR: $LABELS already exists; refusing to overwrite it" >&2; pilot_finish 1; }
mkdir -p "$LABELS"
complete=yes
features=()
for dir in "$PILOT_DIR"/sim/acc*; do
  setting="$(basename "$dir")"
  accuracy="$(echo "$setting" | sed -E 's/acc([0-9]+)_.*/0.\1/')"
  length="${setting##*_len}"
  for aligner in PanAligner minigraph; do
    gaf="$dir/$aligner.gaf"
    if [[ ! -s "$gaf" ]]; then
      echo "MISSING: $gaf"
      complete=no
      continue
    fi
    out="$LABELS/$setting.$aligner.csv.gz"
    timed "label_${setting}_$aligner" "$PYTHON_BIN" "$WG/label_alignments.py" --reads "$dir/all.reads.fa" --gaf "$gaf" \
      --alignments "$PREPARED"/haplotype_alignments/* --reference-names "$PREPARED/reference_names.txt" \
      --ambiguous-bed "$INPUTS"/HG002.*_alllowmapandsegdupregions.bed.gz --difficult-bed "$INPUTS"/HG002.*_alldifficultregions.bed.gz \
      --setting "$setting" --accuracy "$accuracy" --read-length "$length" --seed 11 --aligner "$aligner" --output "$out" || { complete=no; continue; }
    features+=("$out")
  done
done
[[ ${#features[@]} -gt 0 ]] || { echo "ERROR: nothing was labelled" >&2; pilot_finish 1; }
timed kmer_check "$PYTHON_BIN" "$WG/kmer_check.py" --features "${features[@]}" --haplotypes "$PREPARED"/haplotypes/*.fa --reference "$PREPARED/reference.fa" --output "$LABELS/kmer_check.csv.gz"
kmer_args=(); [[ -s "$LABELS/kmer_check.csv.gz" ]] && kmer_args=(--kmer-check "$LABELS/kmer_check.csv.gz")
repro_args=(); [[ -s "$PILOT_DIR/repro.json" ]] && repro_args=(--repro-json "$PILOT_DIR/repro.json")
timed report "$PYTHON_BIN" "$WG/pilot_report.py" --features "${features[@]}" "${kmer_args[@]}" "${repro_args[@]}" --runs-complete "$complete" --output "$RUN_DIR/pilot_report.json"
cp "$PILOT_DIR"/sim/simulation_manifest.json "$RUN_DIR/" 2>/dev/null
pilot_finish 0
