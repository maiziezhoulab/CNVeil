#!/usr/bin/env bash
set -euo pipefail
SAMPLE_ID="TN1"
cb_bam="/path/to/merged_single"
outbam="$cb_bam/all_tumor_cells_norg.bam"
REF="/path/to/hg19/genome.fa"
chisel_result="/path/to/CHISEL/TN1/calls"
VARTRIX="/path/to/vartrix_linux"
tail -n +2 "$chisel_result/calls.tsv" | awk '{print $4}' | sort -u | \
    sed 's/$/-0/' > "$chisel_result/barcodes.tsv"
mkdir -p "$cb_bam/chr"
for i in {1..22}; do
    vcf_in="$cb_bam/${SAMPLE_ID}_phased_chr${i}.vcf"
    vcf_chr="$cb_bam/${SAMPLE_ID}_phased_chr${i}_chrprefix.vcf"
    awk 'BEGIN {OFS="\t"}
         /^##contig=<ID=/ {if ($0 !~ /ID=chr/) sub("ID=", "ID=chr"); print; next}
         /^#/ {print; next}
         {if ($1 !~ /^chr/) $1="chr"$1; print}' "$vcf_in" > "$vcf_chr"
    # Original filename convention retained; see README before interpreting REF/ALT.
    "$VARTRIX" -v "$vcf_chr" -s coverage -b "$outbam" -f "$REF" \
        -c "$chisel_result/barcodes.tsv" \
        --out-matrix "$cb_bam/chr/chr${i}_matrix_ref" \
        --ref-matrix "$cb_bam/chr/chr${i}_matrix_alt"
done
