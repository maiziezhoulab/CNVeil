#!/usr/bin/env bash
set -euo pipefail
# Activate CNRein before running.
BAM_FILE_RAW="/path/to/all_tumor_cells.bam"
BAM_FILE_INT="/path/to/all_tumor_cells_rgcb.bam"
BAM_FILE="/path/to/all_tumor_cells_rgcb_sorted.bam"
REF_FOLDER="/path/to/refNew/ACT"
OUTPUT_DIR="/path/to/CNRein/TN2"
CNREIN_REPO="/path/to/CNRein"
FIX_BAM_HEADER="/path/to/fix_bam_header.py"
REF_GENOME="hg19"
python "$FIX_BAM_HEADER" "$BAM_FILE_RAW" "$BAM_FILE_INT" --threads 20
samtools sort -@8 -o "$BAM_FILE" "$BAM_FILE_INT"
samtools index "$BAM_FILE"
bamLoc="$BAM_FILE"
refLoc="$REF_FOLDER"
outLoc="$OUTPUT_DIR"
refGenome="$REF_GENOME"
mkdir -p "$outLoc"/{temp,counts,info,phased,phasedCounts,readCounts/pos,initial,binScale,model}
for chr in {1..22}; do mkdir -p "$outLoc/readCounts/pos/$chr"; done
all_done=true
for chr in {1..22}; do
    if [[ ! -f $outLoc/counts/ignore_chr${chr}.vcf.gz || ! -f $outLoc/counts/ignore_chr${chr}.vcf.gz.tbi ]]; then all_done=false; break; fi
done
if ! $all_done; then
    parallel -j 18 "bcftools mpileup --ignore-RG -Ou -R '$refLoc/vcf_${refGenome}/chr{}.vcf.gz' -f '$refLoc/${refGenome}.fa' '$bamLoc' | bcftools call -vmO z -o '$outLoc/counts/ignore_chr{}.vcf.gz' && bcftools index -t '$outLoc/counts/ignore_chr{}.vcf.gz'" ::: {1..22}
fi
all_done=true
for chr in {1..22}; do
    if [[ ! -f $outLoc/phased/phased_chr${chr}.bcf ]]; then all_done=false; break; fi
done
if ! $all_done; then
    parallel -j 18 "shapeit4 --input '$outLoc/counts/ignore_chr{}.vcf.gz' --reference '$refLoc/vcf_${refGenome}/chr{}.vcf.gz' --region chr{} --output '$outLoc/phased/phased_chr{}.bcf' --thread 8 && bcftools index '$outLoc/phased/phased_chr{}.bcf'" ::: {1..22}
fi
all_done=true
for chr in {1..22}; do
    if [[ ! -f $outLoc/phased/restricted_chr${chr}.vcf.gz ]]; then all_done=false; break; fi
done
if ! $all_done; then
    parallel -j 18 "bcftools view -i 'GT=\"0|1\" || GT=\"1|0\"' '$outLoc/phased/phased_chr{}.bcf' -Oz -o '$outLoc/phased/restricted_chr{}.vcf.gz' && bcftools index -t '$outLoc/phased/restricted_chr{}.vcf.gz'" ::: {1..22}
fi
all_done=true
for chr in {1..22}; do
    if [[ ! -f $outLoc/counts/seperates_chr${chr}.vcf.gz ]]; then all_done=false; break; fi
done
if ! $all_done; then
    parallel -j 18 "bcftools mpileup --annotate FORMAT/AD -Ou -R '$outLoc/phased/restricted_chr{}.vcf.gz' -f '$refLoc/${refGenome}.fa' '$bamLoc' | bcftools call -vmO z -o '$outLoc/counts/seperates_chr{}.vcf.gz' && bcftools index -t '$outLoc/counts/seperates_chr{}.vcf.gz'" ::: {1..22}
fi
cd "$CNREIN_REPO"
python script.py -input "$BAM_FILE" -ref "$REF_FOLDER" \
    -output "$OUTPUT_DIR" -refGenome "$REF_GENOME"
