#!/usr/bin/env bash
set -euo pipefail

patient=""
bam=""
ref=""
barcodes=""
vcf_dir=""
out_dir=""
vartrix="vartrix_linux"
chromosomes=$(seq 1 22)

usage() {
    echo "Usage:"
    echo "  bash run_vartrix.sh \\"
    echo "    --patient <sample_id> \\"
    echo "    --bam <possorted_bam.bam> \\"
    echo "    --ref <genome.fa> \\"
    echo "    --barcodes <barcodes.txt> \\"
    echo "    --vcf-dir <vcf_folder> \\"
    echo "    --out-dir <output_folder> \\"
    echo "    [--vartrix <path_to_vartrix>]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --patient) patient="$2"; shift 2 ;;
        --bam) bam="$2"; shift 2 ;;
        --ref) ref="$2"; shift 2 ;;
        --barcodes) barcodes="$2"; shift 2 ;;
        --vcf-dir) vcf_dir="$2"; shift 2 ;;
        --out-dir) out_dir="$2"; shift 2 ;;
        --vartrix) vartrix="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$patient" || -z "$bam" || -z "$ref" || -z "$barcodes" || -z "$vcf_dir" || -z "$out_dir" ]]; then
    echo "Missing required arguments."
    usage
    exit 1
fi

mkdir -p "$out_dir"

echo "Running VarTrix for sample: $patient"

for chr in $chromosomes; do
    vcf="${vcf_dir}/${patient}_phased_chr${chr}.vcf.gz"

    if [[ ! -f "$vcf" ]]; then
        echo "Skipping chr${chr}, VCF not found."
        continue
    fi

    echo "Processing chr${chr}"

    "$vartrix" \
        -v "$vcf" \
        -s coverage \
        -b "$bam" \
        -f "$ref" \
        -c "$barcodes" \
        --out-matrix "${out_dir}/chr${chr}_matrix_alt" \
        --ref-matrix "${out_dir}/chr${chr}_matrix_ref"
done

echo "Done."

if [[ -n "${SLURM_JOBID:-}" ]]; then
    echo "SLURM_JOBID: $SLURM_JOBID"
fi
