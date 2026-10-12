#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'Usage: bash run.sh HG19_BED_DIR GINKGO_REPO NEW_OUTPUT_DIR' >&2; exit 2; }
input=$1; repo=$2; out=$3
input=$(cd "$input" && pwd); repo=$(cd "$repo" && pwd)
[[ -s $repo/cli/ginkgo.sh && -s $repo/genomes/hg19/original/variable_500000_76_bowtie ]] || {
    echo 'Missing Ginkgo CLI or hg19 variable_500000_76_bowtie bins' >&2; exit 1;
}
[[ ! -e $out ]] || { echo 'Use a new output directory' >&2; exit 1; }
shopt -s nullglob
files=("$input"/*_hg19_rmdup.bed.gz)
(( ${#files[@]} > 0 )) || { echo 'No *_hg19_rmdup.bed.gz files found' >&2; exit 1; }
mkdir -p "$out"
out=$(cd "$out" && pwd)
[[ $out != *[[:space:]]* ]] || { echo 'Ginkgo upstream requires paths without whitespace' >&2; exit 1; }
# The native CLI may rewrite BED files and writes results beside them.
# Use copies rather than modifying the original input dataset.
for file in "${files[@]}"; do cp "$file" "$out/"; done
printf '%s\n' "$out"/*_hg19_rmdup.bed.gz > "$out/cell.list"
cd "$out"
bash "$repo/cli/ginkgo.sh" --input "$out" --genome hg19 \
    --binning variable_500000_76_bowtie --cells "$out/cell.list" \
    2>&1 | tee "$out/ginkgo.log"
[[ -s $out/SegCopy ]] || { echo 'Ginkgo did not produce SegCopy; inspect ginkgo.log' >&2; exit 1; }
