#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'Usage: bash run.sh HG19_BARCODED_BAM HG19_FASTA NEW_OUTPUT_DIR' >&2; exit 2; }
bam=$1; ref=$2; out=$3
command -v conda >/dev/null
for file in "$bam" "$ref"; do [[ -s $file ]] || { echo "Missing input: $file" >&2; exit 1; }; done
[[ -s ${ref}.fai ]] || { echo 'Missing FASTA index' >&2; exit 1; }
[[ -s ${bam}.bai || -s ${bam%.bam}.bai || -s ${bam}.csi ]] || { echo 'Missing BAM index' >&2; exit 1; }
[[ ! -e $out ]] || { echo 'Use a new output directory' >&2; exit 1; }
bam_dir=$(cd "$(dirname "$bam")" && pwd); bam="$bam_dir/$(basename "$bam")"
ref_dir=$(cd "$(dirname "$ref")" && pwd); ref="$ref_dir/$(basename "$ref")"
mkdir -p "$out"
out=$(cd "$out" && pwd)
run_env() {
    local env=$1; shift
    if [[ $env == */* ]]; then conda run --no-capture-output -p "$env" "$@"
    else conda run --no-capture-output -n "$env" "$@"; fi
}
cd "$out"
run_env "${CHISEL_ENV:-chisel}" chisel_rdr -t "$bam" -r "$ref" \
    -b 50000 -m 100000 -j 8 2>&1 | tee chisel_rdr.log
[[ -s rdr/rdr.tsv ]] || { echo 'Missing CHISEL RDR output' >&2; exit 1; }
gzip -c rdr/rdr.tsv > sprinter.input.tsv.gz
# Keep native intermediates; the original broad rm -rf is omitted.
run_env "${SPRINTER_ENV:-sprinter}" sprinter "$out/sprinter.input.tsv.gz" \
    --refgenome "$ref" --minreads 100000 --rtreads 200 --cnreads 1000 \
    --minnumcells 15 -j 8 2>&1 | tee sprinter.log
