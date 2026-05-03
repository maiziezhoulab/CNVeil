from smooth import smooth_total_cn
from EM import allele_CN
from subprocess import Popen
import os
import glob
def load_snp(vcf_dir, out_dir, sample, num_threads):
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    cmd = f"Rscript {SCRIPT_DIR}/load_snp.r {vcf_dir}  {out_dir} {sample} {num_threads}"
    Popen(cmd, shell= True).wait()



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
        "--var_dir",
        required=True,
        help="Directory containing input VCF/RDS-derived files."
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
        choices=["tumor", "tumor-normal",],
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

    return parser.parse_args()


args = parse_args()

work_dir = os.path.abspath(args.work_dir)
vcf_dir = os.path.abspath(args.var_dir)
ref_type = args.ref_type[2:]
convert = args.convert
snp_source = args.snp_source
snp_file = args.snp_file
num_threads = args.num_threads

total_cn_dir = work_dir + "/Total_CN/"
allele_cn_dir = work_dir + "/Allele_CN/"
baf_dir = allele_cn_dir +"/allele_count_by_cell/"
sample = glob.glob(f"{total_cn_dir}/CNV_*.csv").split('_')[1].split('.')[0]

smooth_total_cn(total_cn_dir, ref_type)
load_snp(vcf_dir, baf_dir, sample, num_threads)
allele_CN(total_cn_dir, allele_cn_dir, baf_dir, sample, convert, snp_source, snp_file, num_threads )
