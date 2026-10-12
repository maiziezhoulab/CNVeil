#!/usr/bin/env bash
# Inputs are existing indexed BAM/FASTA and phased, uncompressed VCFs.
set -euo pipefail
if [[ $# -ne 6 ]]; then
    echo "Usage: bash prepare_vartrix.sh SAMPLE BAM FASTA CHISEL_CALLS VCF_DIR OUTPUT_DIR" >&2
    exit 2
fi
sample=$1; bam=$2; fasta=$3; calls=$4; vcf_dir=$5; out=$6
vartrix=${VARTRIX:-vartrix}
threads=${SLURM_CPUS_PER_TASK:-1}
command -v "$vartrix" >/dev/null || { echo "Vartrix executable unavailable: $vartrix" >&2; exit 1; }
for path in "$bam" "$fasta" "$calls"; do
    [[ -s "$path" ]] || { echo "Missing input: $path" >&2; exit 1; }
done
[[ -s "${fasta}.fai" ]] || { echo "Missing FASTA index" >&2; exit 1; }
[[ -s "${bam}.bai" || -s "${bam%.bam}.bai" || -s "${bam}.csi" ]] || { echo "Missing BAM index" >&2; exit 1; }
for chr in {1..22}; do
    [[ -s "$vcf_dir/${sample}_phased_chr${chr}.vcf" ]] || { echo "Missing chromosome $chr VCF" >&2; exit 1; }
done
# A new directory prevents mixing old, reversed matrices with corrected output.
mkdir "$out"
mkdir "$out/chr" "$out/vcf"
# TN1 convention: CHISEL column 4 plus '-0'. Must match BAM CB tags exactly.
awk 'NR > 1 {print $4}' "$calls" | LC_ALL=C sort -u | sed 's/$/-0/' > "$out/barcodes.tsv"
[[ -s "$out/barcodes.tsv" ]] || { echo "No cell barcodes found" >&2; exit 1; }
for chr in {1..22}; do
    vcf="$out/vcf/${sample}_phased_chr${chr}_chrprefix.vcf"
    awk 'BEGIN {OFS="\t"}
        /^##contig=<ID=/ {if ($0 !~ /^##contig=<ID=chr/) sub("ID=", "ID=chr"); print; next}
        /^#/ {print; next}
        {if ($1 !~ /^chr/) $1="chr"$1; print}' \
        "$vcf_dir/${sample}_phased_chr${chr}.vcf" > "$vcf"
    "$vartrix" -v "$vcf" -s coverage -b "$bam" -f "$fasta" \
        -c "$out/barcodes.tsv" --threads "$threads" \
        --out-matrix "$out/chr/chr${chr}_matrix_alt" \
        --ref-matrix "$out/chr/chr${chr}_matrix_ref"
    echo "Chromosome $chr completed"
done
