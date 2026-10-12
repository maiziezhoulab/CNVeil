#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 6 ]] || {
    echo 'Usage: bash run.sh CELL_BAM_DIR HG19_FASTA HG19_MAPPABILITY PHASED_VCF HG19_NORMAL_BAM NEW_OUTPUT_DIR' >&2; exit 2;
}
input=$1; ref=$2; map=$3; vcf=$4; normal=$5; out=$6
command -v seacon >/dev/null
[[ -d $input ]] || { echo 'Missing cell BAM directory' >&2; exit 1; }
input=$(cd "$input" && pwd)
for file in "$ref" "$map" "$vcf" "$normal"; do [[ -s $file ]] || { echo "Missing input: $file" >&2; exit 1; }; done
[[ -s ${ref}.fai && ( -s ${vcf}.tbi || -s ${vcf}.csi ) ]] || { echo 'Missing FASTA or phased-VCF index' >&2; exit 1; }
[[ -s ${normal}.bai || -s ${normal%.bam}.bai || -s ${normal}.csi ]] || { echo 'Missing normal BAM index' >&2; exit 1; }
[[ ! -e $out ]] || { echo 'Use a new output directory' >&2; exit 1; }
mkdir -p "$out"
seacon prep_readcount -o "$out" -i "$input" -r "$ref" -b 500000 \
    --map-file "$map" -P 12 --normal-bam "$normal" 2>&1 | tee "$out/prep_readcount.log"
seacon prep_baf -o "$out" -i "$input" --vcf "$vcf" --block-size 50 \
    --normal-bam "$normal" 2>&1 | tee "$out/prep_baf.log"
seacon call -o "$out" --tolerance 0.2 --upper-filter 5 --max-wgd 2 \
    2>&1 | tee "$out/call.log"
