#!/usr/bin/env bash
set -euo pipefail
SAMPLE_ID="TN1"
OUTDIR="/path/to/merged_single"
REFERENCE_DIR="/path/to/1000Genome"
GENETIC_MAP="/path/to/Eagle/genetic_map_hg19_withX.txt"
EAGLE="/path/to/Eagle/eagle"
RENAME_FILE="/path/to/rename.txt"
THREADS=10
# Starts from the existing GATK VCF, as did the active original script.
VCF_OUT="$OUTDIR/${SAMPLE_ID}_gatk.vcf"
VCF_GZ="$OUTDIR/${SAMPLE_ID}_gatk.vcf.gz"
VCF_NOCHR="$OUTDIR/${SAMPLE_ID}_gatk.nochr.vcf.gz"
bgzip -c "$VCF_OUT" > "$VCF_GZ"
tabix -f -p vcf "$VCF_GZ"
bcftools annotate --rename-chrs "$RENAME_FILE" -Oz -o "$VCF_NOCHR" "$VCF_GZ"
tabix -f -p vcf "$VCF_NOCHR"
for chr in {1..22}; do
    target="$OUTDIR/${SAMPLE_ID}_gatk.chr${chr}.vcf.gz"
    bcftools view -r "$chr" "$VCF_NOCHR" -Oz -o "$target"
    tabix -f -p vcf "$target"
done
for chr in {1..22}; do
    "$EAGLE" \
        --vcfRef "$REFERENCE_DIR/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
        --geneticMapFile "$GENETIC_MAP" --vcfOutFormat v \
        --vcfTarget "$OUTDIR/${SAMPLE_ID}_gatk.chr${chr}.vcf.gz" \
        --outPrefix "$OUTDIR/${SAMPLE_ID}_phased_chr${chr}" --numThreads "$THREADS"
done
