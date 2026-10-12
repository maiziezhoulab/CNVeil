#!/usr/bin/env bash
set -euo pipefail
# Activate the original CHISEL environment before running.
patient="TN1"
cb_bam="/path/to/merged_single"
output_path="/path/to/CHISEL/TN1"
reference="/path/to/hg19/genome.fa"
FINAL_OUTPUT="$cb_bam/all_chromosomes.tsv"
# Original VCF-to-TSV conversion; no additional variant/genotype filtering.
> "$FINAL_OUTPUT"
for chr in {1..22}; do
    bcftools query -f '%CHROM\t%POS\t[%GT]\n' "$cb_bam/${patient}_phased_chr${chr}.vcf" | \
        awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3}' >> "$FINAL_OUTPUT"
done
mkdir -p "$output_path"
cd "$output_path"
chisel_nonormal -x "$output_path" -t "$cb_bam/all_tumor_cells.bam" \
    -r "$reference" -l "$FINAL_OUTPUT" --seed 12 -b 5Mb
