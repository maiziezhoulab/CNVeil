from pyfaidx import Fasta
import pyranges as pr
import os
import pandas as pd
import pysam
import numpy as np
from scipy.stats import median_abs_deviation
import glob 
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse

global code_dir
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

def get_bam_bed(  reference_fasta, hgref="hg19", resolution=500, sex=False):
	valid_references = ["hg19", "hg38", "mm10"]

	if hgref not in valid_references:
		raise ValueError("Reference genome should be one of: hg19, hg38, or mm10.")

	# reference_fasta = f"path/to/{hgref}_genome.fasta"  # Replace with the actual path to your genome FASTA file

	if not os.path.exists(reference_fasta):
		raise FileNotFoundError(f"Reference genome FASTA file {reference_fasta} not found for {hgref}.")

	genome = Fasta(reference_fasta)

	if resolution <= 0:
		raise ValueError("Invalid fixed bin length.")

	bins = {}
	
	for chrom in genome.keys():
		chrom_length = len(genome[chrom])
		bin_size = resolution * 1000
		num_bins = (chrom_length + bin_size - 1) // bin_size
		# bin_positions = [(i * bin_size+1, min((i + 1) * bin_size, chrom_length)) for i in range(num_bins)]
		bin_positions = [(i * bin_size+1, (i + 1) * bin_size) for i in range(num_bins)]
		bins[chrom] = bin_positions
		


	autochr = 22 if hgref != "mm10" else 19

	if sex:
		sex_chromosomes = [f"chr{str(i)}" for i in range(1, autochr + 1)] + ["chrX", "chrY"]
		ref = {chrom: positions for chrom, positions in bins.items() if chrom in sex_chromosomes}



	else:
		autosomes = [f"chr{str(i)}" for i in range(1, autochr + 1)]
		ref = {chrom: positions for chrom, positions in bins.items() if chrom in autosomes}

	ref = {f"chr{chrom}" if not chrom.startswith("chr") else chrom: positions for chrom, positions in ref.items()}

	n_bin = 0
	interval_list = []
	bin_str_list = []
	for chrom,poss in ref.items():
		for start, end in poss:
			interval = {
				'Chromosome': chrom,
				'Start': start,
				'End': end
			}
			interval_list.append(interval)
			bin_str_list.append(f"{chrom}:{start}-{end}")

		n_bin += len(poss)
	# print(interval_list)
	# Step 2: Convert the list of dictionaries into a pyranges object
	intervals = pr.from_dict(interval_list)

	# Display the pyranges object
	# print(intervals)

	print(f"num of raw bins by resolution of {resolution}k: {n_bin}")
	return {  "ref": intervals, "bin": bin_str_list}

def get_map(ref_type,use_bin_index):
    map_file = code_dir+f"/{ref_type}_ref_file/genome_mappability.tab"

    with open(map_file,'r') as f:
        map_list = [ eval(x) for x in f.read().split('\n')[:-1]]
    map_list = [map_list[i] for i in use_bin_index]

    return np.array(map_list)

def get_GC(ref_type,use_bin_index):
	gc_file = code_dir+f"/{ref_type}_ref_file/genome_gc.bed"

	with open(gc_file,'r') as f:
		gc_list = [ eval(x.split()[-1]) for x in f.read().split('\n')[:-1]]
	gc_list = [gc_list[i] for i in use_bin_index]
	return np.array(gc_list)


def MAD(data):
    # Assuming data is your NumPy array
    median = np.median(data)
    mad = 1.4826 * np.median(np.abs(data - median))
    return mad

def filter_bin_cell(Y_raw, samples,ref_raw,  mapp, gc,
               gc_thresh=(20, 80), mapp_thresh=0.9, nMAD=3, filter_by_Nmat=False, keep_low_ploidy_cell = False):
    
	# Check if input dimensions match
	if len(ref_raw) != Y_raw.shape[0]:
		raise ValueError("Invalid inputs: length of ref and # of rows in read count matrix must be the same")

	# Apply filters on bins based on GC content and mappability
	binfilter1 = (gc < gc_thresh[0]) | (gc > gc_thresh[1])
	print(f"Excluded {sum(binfilter1)} bins due to extreme GC content.")

	binfilter2 = (mapp < mapp_thresh)
	print(f"Excluded {sum(binfilter2)} bins due to low mappability.")

	# Exclude bins that failed any of the above filters
	bin_filter_mask = binfilter1 | binfilter2
	if bin_filter_mask.any():
		ref = ref_raw.loc[~bin_filter_mask, :]
		Y = Y_raw[~bin_filter_mask, :]
	else:
		ref = ref_raw
		Y = Y_raw

	# Compute normalization factors
	filtered_rows = np.apply_along_axis(lambda x: not np.any(x == 0), axis=1, arr=Y)
	# filtered_rows = np.apply_along_axis(lambda x: np.sum(x == 0) < 0.1 * len(x), axis=1, arr=Y)

	Y_nonzero = Y[filtered_rows]

	print("Y.shape", Y.shape)
	print("Y_nonzero.shape", Y_nonzero.shape)
	# exit()

	if Y_nonzero.shape[0] <= 10:
		print("Adopt arithmetic mean instead of geometric mean")
		# exit()
		pseudo_sample = Y.mean(axis=1)
		N = np.median(Y / pseudo_sample[:, np.newaxis], axis=0)
	else:
		pseudo_sample = np.array([np.exp(np.sum(np.log(row)) / len(row)) for row in Y_nonzero]).reshape(-1, 1)
		# pseudo_sample = np.array([row.mean() for row in Y_nonzero]).reshape(-1, 1)
		# print(pseudo_sample)
		N = np.nanmedian(Y_nonzero / pseudo_sample, axis=0)

	print("N shape:", N.shape)

	# Normalize and filter bins based on deviation from median
	n_rows, n_cols = Y.shape
	Nmat = np.tile(N, (n_rows, 1))
	# print("======",Y.shape, Nmat.shape)
	# filter_by_Nmat = 0
	if filter_by_Nmat:
		print("filter by Nmat...")
		bin_sum = (Y / Nmat).sum(axis=1)
		# print(bin_sum.mean())

		# exit()
		binfilter3 = (bin_sum >= (np.median(bin_sum) - nMAD * MAD(bin_sum))) & \
						(bin_sum <= (np.median(bin_sum) + 5 * MAD(bin_sum)))
		# binfilter3 = (bin_sum >= (np.median(bin_sum) - nMAD * MAD(bin_sum))) 
		## ---------- use nmad bin filter
		Y = Y[binfilter3, :]
		ref = ref.loc[binfilter3, :].reset_index(drop=True)

		print(f"Filtered {n_rows - Y.shape[0]} bins based on deviation from median.")
	# exit()

	#------- filter cells
	if not keep_low_ploidy_cell:
		cell_filter = (Y==0).mean(0)<0.05
		Y = Y [:, cell_filter]
		print(f"Filtered {n_cols - Y.shape[1]} cells based on number of 0 cov bins.")

	print(f"There are {Y.shape[1]} samples and {Y.shape[0]} bins after QC step.")

	# 
	# df1 = pd.concat([ref, df], axis=1)
	# df1.columns = ['chr', 'start', 'end'] + [f"Sample_{i}" for i in range(Y.shape[1])]
	ref.columns = ['chr', 'start', 'end']
	ref['bin'] = ref.apply(lambda row: f"{row['chr']}:{row['start']}-{row['end']}", axis=1)
	if not keep_low_ploidy_cell:
		df = pd.DataFrame(Y, columns=np.array(samples)[cell_filter])
	else:
		df = pd.DataFrame(Y, columns=np.array(samples))
	df.insert(0,  'bin',ref['bin'].tolist())

	return df




def perform_qc(rdmatrix, reftype, outfile, filter_by_Nmat,keep_low_ploidy_cell):
	assert reftype in ['hg19','hg38']
	if reftype == 'hg38':
		reference_fasta = "/data/maiziezhou_lab/Softwares/refdata-GRCh38-2.1.0/fasta/genome.fa" 
		
	else:
		reference_fasta = "/data/maiziezhou_lab/Softwares/refdata-hg19-2.1.0/fasta/genome.fa"

	df = pd.read_csv(rdmatrix)
	af_bins = set(df.iloc[:,0].values)
	Y_raw = df.iloc[:,1:].values
	samples = df.columns.tolist()[1:]


	bambed_dc = get_bam_bed( reference_fasta,
							hgref=reftype, resolution=500, sex=False)
	print(f"num of bins in raw RDmatrix:{len(af_bins)}")
	all_bins = bambed_dc['bin']
	use_bin_index = [ i for i in range(len(all_bins)) if all_bins[i] in af_bins]
	map_list = get_map(reftype, use_bin_index)
	gc_list = get_GC(reftype, use_bin_index)

	ref_list = []
	for _,df in bambed_dc['ref'].items():
		ref_list.append(df)
	ref_raw = pd.concat(ref_list)
	ref_raw = ref_raw.iloc[use_bin_index,:].reset_index(drop = True)
	print("num of bins:",ref_raw.shape[0])

	df = filter_bin_cell(Y_raw, samples, ref_raw,  map_list, gc_list, filter_by_Nmat = filter_by_Nmat , keep_low_ploidy_cell = keep_low_ploidy_cell)
	df.to_csv(outfile, index = False)

if __name__ == '__main__':
	reftype = 'hg38'
	rdmatrix = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/Use_AF_matrix/T10/RDmatrix.csv"
	outfile = "/data/maiziezhou_lab/CanLuo/Single_Cell_Project/DEV/Use_AF_matrix/T10/RDmatrix_qc.csv"
	perform_qc(rdmatrix, reftype, outfile, keep_low_ploidy_cell)