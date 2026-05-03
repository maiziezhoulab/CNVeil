from subprocess import Popen
import os

def bin_reads_GC(merged_bam, split_bam_folder, refv, out_dir,cell_node_file, num_threads):
    ''' refv : hg19|hg38 '''
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    cmd = f"Rscript {SCRIPT_DIR}/Preprocess/preprocess.R {merged_bam} {split_bam_folder} {refv} {out_dir}/bin_reads/  {num_threads}"
    Popen(cmd, shell= True).wait()

    if cell_node_file is None:
        cmd = f"Rscript {SCRIPT_DIR}/Preprocess/reformat_bin_reads.R  -i {out_dir}/bin_reads/binned-GC  -o {out_dir}/RDmatrix.csv"
        Popen(cmd, shell= True).wait()
    else:
        cmd = f"Rscript {SCRIPT_DIR}/Preprocess/reformat_bin_reads.R  -i {out_dir}/bin_reads/binned-GC  -o {out_dir}/RDmatrix.csv -t {cell_node_file}"
        Popen(cmd, shell= True).wait()


import argparse
from argparse import ArgumentParser
parser = ArgumentParser(description="",
    usage='use "python3 %(prog)s --help" for more information',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--merged_bam_file','-mgbam')
parser.add_argument('--split_bam_dir','-spbam')
parser.add_argument('--work_dir','-w')
parser.add_argument('--reftype','-rt', choices=['hg19','hg38'])
parser.add_argument('--num_threads','-t', type = int, default= 20)
parser.add_argument('--cell_node','-cn',help = "this is designed for the simulation data, not meaningful for real data. cell node file; optional, if not given, will not transform cell to node in the end")
args = parser.parse_args()

merged_bam_file = args.merged_bam_file
split_bam_dir = args.split_bam_dir
output_dir = args.work_dir + "/Total_CN/"
cell_node = args.cell_node
num_threads = args.num_threads
reftype = args.reftype
if not os.path.exists(output_dir):
    os.system("mkdir -p "+output_dir)


bin_reads_GC(merged_bam_file, split_bam_dir, reftype, 
             output_dir, cell_node, num_threads)