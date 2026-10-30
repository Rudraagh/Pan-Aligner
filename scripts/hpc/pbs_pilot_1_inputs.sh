#!/bin/bash
# Pilot stage 1: fetch inputs (only with FETCH=1), record checksums, extract chr21, compare CHM13 versions,
# and produce the haplotype-to-reference alignments. Downloads about 7 GB: needs explicit approval first.
#
# Without FETCH=1, the files must already be in $PILOT_DIR/inputs under the names below.
# ROUTE: auto (default), giab (dipcall BAMs, needs samtools and identical chr21) or minimap2.
# UNVERIFIED cluster settings: confirm queue, memory and internet access with the administrator.
#PBS -N pa_pilot_inputs
#PBS -q cpuq
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -l walltime=04:00:00
#PBS -j oe
set -uo pipefail
source "${PROJECT_DIR:-${PBS_O_WORKDIR:-$(pwd)}}/scripts/hpc/pilot_common.sh"
pilot_start inputs
FETCH="${FETCH:-0}"
ROUTE="${ROUTE:-auto}"
INPUTS="$PILOT_DIR/inputs"
PREPARED="$PILOT_DIR/prepared"
mkdir -p "$INPUTS" "$PREPARED"

HP=https://s3-us-west-2.amazonaws.com/human-pangenomics
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
DIP="$GIAB/AshkenazimTrio/HG002_NA24385_son/v5.0q/dipcall_output/CHM13v2.0_HG2-T2TQ100-V1.1_dipcall-z2k"
STRAT="$GIAB/genome-stratifications/v3.6"
URLS=(
  "$HP/T2T/HG002/assemblies/hg002v1.1.mat.fasta.gz"
  "$HP/T2T/HG002/assemblies/hg002v1.1.pat.fasta.gz"
  "$HP/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-chm13.gfa.gz"
  "$HP/pangenomes/freeze/freeze1/minigraph/CHM13v11Y.fa.gz"
  "$HP/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
  "$DIP.hap1.bam" "$DIP.hap1.bam.bai" "$DIP.hap2.bam" "$DIP.hap2.bam.bai"
  "$STRAT/HG002.mat@all/Union/HG002.mat_alllowmapandsegdupregions.bed.gz"
  "$STRAT/HG002.pat@all/Union/HG002.pat_alllowmapandsegdupregions.bed.gz"
  "$STRAT/HG002.mat@all/Union/HG002.mat_alldifficultregions.bed.gz"
  "$STRAT/HG002.pat@all/Union/HG002.pat_alldifficultregions.bed.gz"
)
printf '%s\n' "${URLS[@]}" > "$RUN_DIR/input_urls.txt"

fail() { echo "ERROR: $*" >&2; pilot_finish 1; }

if [[ "$FETCH" == "1" ]]; then
  for url in "${URLS[@]}"; do
    file="$INPUTS/${url##*/}"
    [[ -s "$file" ]] && { echo "present: $file"; continue; }
    timed "fetch_${url##*/}" curl -fL --retry 3 -o "$file.partial" "$url" || fail "download failed: $url"
    mv "$file.partial" "$file"
  done
fi
for name in hg002v1.1.mat.fasta.gz hg002v1.1.pat.fasta.gz hprc-v1.0-minigraph-chm13.gfa.gz CHM13v11Y.fa.gz HG002.mat_alllowmapandsegdupregions.bed.gz HG002.pat_alllowmapandsegdupregions.bed.gz; do
  [[ -s "$INPUTS/$name" ]] || fail "missing input $INPUTS/$name"
done
(cd "$INPUTS" && sha256sum ./* > "$RUN_DIR/inputs.sha256") || fail "checksums failed"

# Published in the GIAB v3.6 HG002.mat README, for the uncompressed FASTA.
mat_md5="$(gunzip -c "$INPUTS/hg002v1.1.mat.fasta.gz" | md5sum | cut -d' ' -f1)"
echo "hg002v1.1.mat md5 $mat_md5"
[[ "$mat_md5" == "d26bf77757a98e69fdcd43ad0894ec4c" ]] || fail "maternal FASTA does not match the published md5"

gunzip -c "$INPUTS/hprc-v1.0-minigraph-chm13.gfa.gz" > "$PREPARED/graph.gfa.partial" && mv "$PREPARED/graph.gfa.partial" "$PREPARED/graph.gfa" || fail "graph decompression failed"
gunzip -c "$INPUTS/CHM13v11Y.fa.gz" > "$PREPARED/reference.fa.partial" && mv "$PREPARED/reference.fa.partial" "$PREPARED/reference.fa" || fail "reference decompression failed"
awk '/^>/{print substr($1,2)}' "$PREPARED/reference.fa" > "$PREPARED/reference_names.txt"
# Rank-0 names in the graph must match the reference names, or every path would look off-reference.
awk -F'\t' '$1=="S"{for(i=4;i<=NF;i++){if($i=="SR:i:0"){r=1}; if($i ~ /^SN:Z:/){n=substr($i,6)}}; if(r)print n; r=0}' "$PREPARED/graph.gfa" | sort -u > "$RUN_DIR/graph_rank0_names.txt"
comm -13 <(sort -u "$PREPARED/reference_names.txt") "$RUN_DIR/graph_rank0_names.txt" > "$RUN_DIR/graph_rank0_not_in_reference.txt"
if [[ -s "$RUN_DIR/graph_rank0_not_in_reference.txt" ]]; then
  fail "graph rank-0 names are not all in the reference; see graph_rank0_not_in_reference.txt"
fi

timed extract_mat "$PYTHON_BIN" "$WG/prepare_inputs.py" extract --fasta "$INPUTS/hg002v1.1.mat.fasta.gz" --records "${MAT_RECORD:-chr21_MATERNAL}" --output-dir "$PREPARED/haplotypes" || fail "maternal chr21 extraction failed"
timed extract_pat "$PYTHON_BIN" "$WG/prepare_inputs.py" extract --fasta "$INPUTS/hg002v1.1.pat.fasta.gz" --records "${PAT_RECORD:-chr21_PATERNAL}" --output-dir "$PREPARED/haplotypes" || fail "paternal chr21 extraction failed"

identical=false
if [[ -s "$INPUTS/chm13v2.0.fa.gz" ]]; then
  timed compare "$PYTHON_BIN" "$WG/prepare_inputs.py" compare --first "$PREPARED/reference.fa" --second "$INPUTS/chm13v2.0.fa.gz" --record chr21 --output "$RUN_DIR/chr21_comparison.json" || fail "chr21 comparison failed"
  grep -q '"identical": true' "$RUN_DIR/chr21_comparison.json" && identical=true
fi

have_bams=false
[[ -s "$INPUTS/CHM13v2.0_HG2-T2TQ100-V1.1_dipcall-z2k.hap1.bam" && -s "$INPUTS/CHM13v2.0_HG2-T2TQ100-V1.1_dipcall-z2k.hap2.bam" ]] && have_bams=true
if [[ "$ROUTE" == "auto" ]]; then
  if [[ "$identical" == "true" && "$have_bams" == "true" ]] && command -v samtools >/dev/null 2>&1; then ROUTE=giab; else ROUTE=minimap2; fi
fi
echo "projection route: $ROUTE (chr21 identical: $identical; BAMs: $have_bams; samtools: $(command -v samtools || echo none))" | tee "$RUN_DIR/projection_route.txt"

ALIGN_DIR="$PREPARED/haplotype_alignments"
[[ -e "$ALIGN_DIR" ]] && fail "$ALIGN_DIR already exists; refusing to overwrite it"
mkdir -p "$ALIGN_DIR"
if [[ "$ROUTE" == "giab" ]]; then
  [[ "$identical" == "true" ]] || fail "route giab needs identical chr21"
  for hap in hap1 hap2; do
    # Keep records from either chr21 haplotype; label_alignments.py stops if a haplotype has none.
    samtools view -@ "$THREADS" "$INPUTS/CHM13v2.0_HG2-T2TQ100-V1.1_dipcall-z2k.$hap.bam" | awk -F'\t' '$1 ~ /^chr21_/' | gzip -c > "$ALIGN_DIR/$hap.sam.gz" || fail "samtools failed on $hap"
  done
elif [[ "$ROUTE" == "minimap2" ]]; then
  cat "$PREPARED"/haplotypes/*.fa > "$PREPARED/haplotypes_all.fa"
  timed minimap2 "$BIN_DIR/minimap2" -cx asm5 -t "$THREADS" --secondary=no -o "$ALIGN_DIR/haplotypes.paf" "$PREPARED/reference.fa" "$PREPARED/haplotypes_all.fa" || fail "minimap2 failed"
else
  fail "unknown ROUTE=$ROUTE"
fi
ls -la "$ALIGN_DIR" > "$RUN_DIR/alignment_files.txt"
pilot_finish 0
