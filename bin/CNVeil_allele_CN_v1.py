from smooth import smooth_total_cn
from EM import allele_CN
from subprocess import Popen
import os
import glob


def load_snp(vcf_dir, out_dir, sample, num_threads):
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    cmd = f"Rscript {SCRIPT_DIR}/load_snp.r {vcf_dir} {out_dir} {sample} {num_threads}"
    print(cmd)
    Popen(cmd, shell=True).wait()


#!/usr/bin/env python3

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ASCN/BAF phasing pipeline."
    )

    parser.add_argument(
        "--work_dir",
        required=True,
        help="Working directory for CNVeil pipeline."
    )

    parser.add_argument(
        "--ref_type",
        required=True,
        choices=["hg38", "hg19"],
        help="Reference type."
    )

    parser.add_argument(
        "--convert",
        default=None,
        help="Optional cell name/barcode converter file."
    )

    parser.add_argument(
        "--snp_source",
        choices=["tumor", "tumor-normal"],
        required=True,
        help="SNP source type."
    )

    parser.add_argument(
        "--snp_file",
        default=None,
        help="Optional SNP genotype file."
    )

    parser.add_argument(
        "--num_threads",
        "-t",
        type=int,
        default=22,
        help="Number of parallel threads."
    )

    # -----------------------------
    # New arguments for updated EM.py
    # -----------------------------
    parser.add_argument(
        "--threshold_file",
        default=None,
        help=(
            "Threshold model TSV, e.g. best_mBAF_threshold_by_CN.tsv. "
            "Required/recommended when --snp_source tumor."
        )
    )

    parser.add_argument(
        "--threshold_feature",
        default="median_mBAF",
        choices=["median_mBAF", "mean_mBAF"],
        help="Feature used for threshold-based LOH prediction."
    )

    parser.add_argument(
        "--min_depth",
        type=int,
        default=10,
        help="Minimum depth for precomputed position-level BAF SNPs."
    )

    parser.add_argument(
        "--min_snps",
        type=int,
        default=10,
        help="Minimum number of SNP positions required per segment."
    )

    return parser.parse_args()


args = parse_args()

work_dir = os.path.abspath(args.work_dir)
ref_type = args.ref_type[2:]
convert = args.convert
snp_source = args.snp_source
snp_file = args.snp_file
num_threads = args.num_threads

threshold_file = args.threshold_file
threshold_feature = args.threshold_feature
min_depth = args.min_depth
min_snps = args.min_snps

total_cn_dir = work_dir + "/Total_CN/"
allele_cn_dir = work_dir + "/Allele_CN/"
baf_dir = allele_cn_dir + "/allele_count_by_cell/"

sample = (
    glob.glob(f"{total_cn_dir}/CNV_*.csv")[0]
    .split("/")[-1]
    .replace("CNV_", "")
    .split(".")[0]
)

# smooth_total_cn(total_cn_dir, ref_type)

# -------- deprecated
# load_snp(vcf_dir, baf_dir, sample, num_threads)
# -------- deprecated

# Optional safety check
if snp_source == "tumor" and threshold_file is None:
    raise ValueError(
        "--threshold_file is required when --snp_source tumor "
        "because updated EM.py uses the threshold model to predecide LOH segments."
    )

allele_CN(
    total_cn_dir,
    allele_cn_dir,
    baf_dir,
    sample,
    convert,
    snp_source,
    snp_file,
    num_threads,
    threshold_file=threshold_file,
    threshold_feature=threshold_feature,
    min_depth=min_depth,
    min_snps=min_snps,
    loh_only= False,
    # only_chrom= 'chr7'
)