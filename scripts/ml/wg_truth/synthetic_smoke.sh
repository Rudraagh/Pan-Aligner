#!/bin/bash
# Run the whole truth pipeline on synthetic data with the real tools. Checks plumbing only: the numbers
# say nothing about real genomes, and the gate verdict here is expected to be NO_GO (too few errors).
#
# Usage: synthetic_smoke.sh WORK_DIR PYTHON PANALIGNER MINIGRAPH MINIMAP2 PBSIM MODEL
set -euo pipefail

WORK="$1"; PYTHON="$2"; PANALIGNER="$3"; MINIGRAPH="$4"; MINIMAP2="$5"; PBSIM="$6"; MODEL="$7"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -e "$WORK" ]]; then
  echo "ERROR: $WORK already exists; choose a new folder." >&2
  exit 1
fi
mkdir -p "$WORK"
export PYTHONDONTWRITEBYTECODE=1

"$PYTHON" "$HERE/synthetic_fixture.py" "$WORK"
"$MINIGRAPH" -cxggs -t2 "$WORK/reference.fa" "$WORK/graph_hap.fa" > "$WORK/graph.gfa" 2> "$WORK/graph.log"
"$MINIMAP2" -cx asm5 -t2 --secondary=no "$WORK/reference.fa" "$WORK/source.fa" > "$WORK/source_to_ref.paf" 2> "$WORK/minimap2.log"
"$PYTHON" "$HERE/simulate_haplotypes.py" --haplotypes "mat=$WORK/source.fa" --depth 8 --pbsim "$PBSIM" --model "$MODEL" --output-dir "$WORK/sim"

features=()
for setting_dir in "$WORK"/sim/acc*; do
  setting="$(basename "$setting_dir")"
  accuracy="$(echo "$setting" | sed -E 's/acc([0-9]+)_.*/0.\1/')"
  length="${setting##*_len}"
  for aligner in panaligner minigraph; do
    binary="$PANALIGNER"; [[ "$aligner" == minigraph ]] && binary="$MINIGRAPH"
    gaf="$setting_dir/mat.$aligner.gaf"
    "$binary" -t2 -cx lr "$WORK/graph.gfa" "$setting_dir/mat.reads.fa" > "$gaf" 2> "$gaf.log"
    out="$setting_dir/mat.$aligner.features.csv.gz"
    "$PYTHON" "$HERE/label_alignments.py" --reads "$setting_dir/mat.reads.fa" --gaf "$gaf" \
      --alignments "$WORK/source_to_ref.paf" --reference-names "$WORK/reference.fa.fai" \
      --ambiguous-bed "$WORK/ambiguous.bed" --setting "$setting" --accuracy "$accuracy" --read-length "$length" \
      --seed 11 --aligner "$aligner" --output "$out"
    features+=("$out")
  done
done

"$PYTHON" "$HERE/kmer_check.py" --features "${features[@]}" --haplotypes "$WORK/source.fa" --reference "$WORK/reference.fa" --sample 500 --output "$WORK/kmer_check.csv.gz"
"$PYTHON" "$HERE/simulate_haplotypes.py" --haplotypes "mat=$WORK/source.fa" --depth 8 --accuracies 0.95 --lengths 2000 --pbsim "$PBSIM" --model "$MODEL" --output-dir "$WORK/repro"
"$PYTHON" -c "
import json, sys, hashlib
a, b, out = sys.argv[1:4]
digest = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
json.dump({'identical': digest(a) == digest(b)}, open(out, 'x'))
" "$WORK/sim/acc95_len2000/mat.reads.fa" "$WORK/repro/acc95_len2000/mat.reads.fa" "$WORK/repro.json"
"$PYTHON" "$HERE/pilot_report.py" --features "${features[@]}" --kmer-check "$WORK/kmer_check.csv.gz" --runs-complete yes --repro-json "$WORK/repro.json" --output "$WORK/pilot_report.json"
echo "Synthetic smoke finished: $WORK/pilot_report.json"
