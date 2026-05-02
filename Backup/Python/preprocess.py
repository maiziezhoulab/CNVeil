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


global code_dir
code_dir = os.path.dirname(os.path.realpath(__file__))+'/'

def get_bam_bed(bamdir, sampname, reference_fasta, hgref="hg19", resolution=500, sex=False):
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
        bin_positions = [(i * bin_size+1, min((i + 1) * bin_size, chrom_length)) for i in range(num_bins)]
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
    for chrom,poss in ref.items():
        for start, end in poss:
            interval = {
                'Chromosome': chrom,
                'Start': start,
                'End': end
            }
            interval_list.append(interval)

        n_bin += len(poss)

    # Step 2: Convert the list of dictionaries into a pyranges object
    intervals = pr.from_dict(interval_list)

    # Display the pyranges object
    print(intervals)

    print(f"num of raw bins by resolution of {resolution}k: {n_bin}")
    return {"bamdir": bamdir, "sampname": sampname, "ref": intervals}


def get_map(ref_type):
    map_file = code_dir+f"/{ref_type}_ref_file/genome_mappability.tab"

    with open(map_file,'r') as f:
        map_list = [ eval(x) for x in f.read().split('\n')[:-1]]

    return np.array(map_list)

def get_GC(ref_type):
    gc_file = code_dir+f"/{ref_type}_ref_file/genome_gc.bed"

    with open(gc_file,'r') as f:
        gc_list = [ eval(x.split()[-1]) for x in f.read().split('\n')[:-1]]

    return np.array(gc_list)


def get_masked_ref(ref_type):
    mask_file = code_dir+f"/{ref_type}_ref_file/mask.ref.{ref_type}"

    df=pd.read_csv(mask_file,sep = "\t", header = None)
    df.columns = ['chrom','start','end','a','b','c']

    masked_ref = pr.PyRanges(chromosomes=df["chrom"], starts=df["start"] + 1, ends=df["end"] )
    

    print("masked ref:")
    print(masked_ref)
    return masked_ref


def get_cov_one_bam(bamurl,seq,mask_ref,ref,mapqthresh, temp_dir,n_bin):
    if seq == "paired-end":

        # filter_flag = 1024 optical duplicate

        bam_ref = pr.readers.read_bam(bamurl, sparse=True, as_df=False, mapq=int(mapqthresh),
                                required_flag = 65,  )
    else:
        # print("process single end")

        bam_ref = pr.readers.read_bam(bamurl, sparse=True, as_df=False, mapq=int(mapqthresh),
                                      
                                          )

    if len(bam_ref) == 0:
        # Failed library preparation
        print(0, file= open(temp_dir+"/"+bamurl.split('/')[-1].split('_')[0],'w'))
        return [0]*n_bin
        
    else:
        # remove reads that overlap with DUP and GAP mask
        # print(1)s
        cov_mask = bam_ref.count_overlaps(mask_ref)
        # print(2)
        bam_ref = bam_ref[cov_mask.NumberOverlaps==0]

        # count reads in bins
        # print(3)
        coverage = ref.count_overlaps(bam_ref)
        # print(4)
        val = coverage.NumberOverlaps.values
        print(val, file= open(temp_dir+"/"+bamurl.split('/')[-1].split('_')[0],'w'))

        return val
        
def get_coverage_scDNA(bambedObj, mapq_thresh, seq, hgref, n_thread, outdir,n_bin):
    valid_references = ["hg19", "hg38", "mm10"]

    if hgref not in valid_references:
        raise ValueError("Reference genome should be one of: hg19, hg38, or mm10.")
    
    if seq == "paired-end":
        flag = 1  # Set the flag value for paired-end reads
        # flag = pysam.build_flag(is_paired=True, is_unmapped=False, is_duplicate=False, is_qcfail=False, is_firstmate=True)
    elif seq == "single-end":
        flag = 0  # Set the flag value for single-end reads
        # flag = pysam.build_flag(is_paired=False, is_unmapped=False, is_duplicate=False, is_qcfail=False)
    else:
        print("seq type can only be paired-end or single-end")
        exit()


    ref = bambedObj["ref"]
    bamdir = bambedObj["bamdir"]
    sampname = bambedObj["sampname"]

    mask_ref = get_masked_ref(hgref)

    # Convert the PyRanges object to a Pandas DataFrame
    interval_dfs = ref.dfs
    formatted_intervals = []
    # Iterate through the DataFrame rows
    for _,interval_df in interval_dfs.items(): 
        # print(interval_df[1])
        for index, row in interval_df.iterrows():
            chrom = row["Chromosome"]
            start = row["Start"]
            end = row["End"]
            formatted_interval = f"{chrom}:{start}-{end}"
            formatted_intervals.append(formatted_interval)


    Y = pd.DataFrame(index=formatted_intervals, columns=sampname)

    temp_dir = outdir+"/temp/"
    os.system("mkdir -p " + temp_dir)

    if n_thread == 1:
        sequences = []

        for i in range(len(sampname)):
            bamurl = bamdir[i]
            print(f"Getting coverage for sample {i + 1}: {sampname[i]}...")
            val = get_cov_one_bam(bamurl,seq,mask_ref,ref, mapq_thresh, temp_dir,n_bin)
            print(val)
            sequences.append(val)
            # Y[sampname[i]] = val
    
    else:
        sequences = Parallel(n_jobs=n_thread)(delayed(get_cov_one_bam)\
                                              (bamurl,seq,mask_ref,ref,mapq_thresh, temp_dir,n_bin) \
                                              for bamurl in tqdm(bamdir, desc = 'Get coverage'))
    # print(sequences)
    Y = pd.DataFrame(np.array(sequences).T, columns = sampname, index = formatted_intervals)
        # for i in range(len(sampname)):
        #     Y[sampname[i]] = sequences[i]

    return Y

def qc_one_bam(bam_file_path):
    # Open BAM file
    samfile = pysam.AlignmentFile(bam_file_path, "rb")

    # Initialize variables to calculate metrics
    total_reads = 0
    mapped_reads = 0
    non_duplicated_mapped_reads = 0
    mapq20_reads = 0
    total_read_length = 0

    # Iterate through reads in the BAM file
    for read in samfile:
        total_reads += 1
        total_read_length += read.query_length
        if not (read.reference_name is None):
            mapped_reads += 1
            if not read.is_duplicate:
                non_duplicated_mapped_reads += 1
            if read.mapping_quality >= 20:
                mapq20_reads += 1
    


    # Calculate QC metrics
    read_length_avg = total_read_length / total_reads if total_reads > 0 else 0
    mapped_prop = mapped_reads / total_reads if total_reads > 0 else 0
    non_dup_prop = non_duplicated_mapped_reads / total_reads if total_reads > 0 else 0
    mapq20_prop = mapq20_reads / total_reads if total_reads > 0 else 0

    # Close BAM file
    samfile.close()
    dat = [read_length_avg, total_reads, mapped_reads, mapped_prop,
                        non_duplicated_mapped_reads, non_dup_prop, mapq20_reads, mapq20_prop]
    return dat


def get_samp_QC(bambedObj, n_thread):
    # Extracting data from the bambedObj
    ref = bambedObj["ref"]
    bamdir = bambedObj["bamdir"]
    sampname = bambedObj["sampname"]

    # Initialize QC metric matrix
    num_samples = len(sampname)
    QCmetric = np.zeros((num_samples, 8), dtype=float)
    colnames = ["readlength", "total", "mapped", "mapped_prop", "non_dup", "non_dup_prop", "mapq20", "mapq20_prop"]

    if n_thread == 1:
        # Iterate through each sample to compute QC metrics
        for i, bam_file_path in enumerate(bamdir):
            print(f"Getting sample QC metric for sample {i + 1}")
            dat = qc_one_bam(bam_file_path)
            # Populate the QC metric matrix
            QCmetric[i, :] = dat
    else:
        sequences = Parallel(n_jobs=n_thread)(delayed(qc_one_bam)(bam_file_path) \
                                              for bam_file_path in tqdm(bamdir, desc = 'QC'))
        for i, bam_file_path in enumerate(bamdir):
            QCmetric[i, :] = sequences[i]


    # Create a dictionary with column names
    result = {}
    for i, colname in enumerate(colnames):
        result[colname] = QCmetric[:, i]

    df = pd.DataFrame(result)
    return df


def MAD(data):
    # Assuming data is your NumPy array
    median = np.median(data)
    mad = 1.4826 * np.median(np.abs(data - median))
    return mad

def perform_qc(Y_raw, sampname_raw, ref_raw, QCmetric_raw,mapp,gc,
               cov_thresh=0, minCountQC=20, mapq20_thresh=0.3,
               mapp_thresh=0.9, gc_thresh=(20, 80), nMAD=3):
    
    # Check if input dimensions match
    if len(ref_raw) != Y_raw.shape[0]:
        raise ValueError("Invalid inputs: length of ref and # of rows in read count matrix must be the same")
    if len(sampname_raw) != Y_raw.shape[1]:
        raise ValueError("Invalid inputs: length of sample names and # of cols in read count matrix must be the same")
    if QCmetric_raw.shape[0] != Y_raw.shape[1]:
        raise ValueError("Invalid inputs: # of rows in QC metric and # of cols in read count matrix must be the same")

    # Extract mapping and GC content information
    # mapp = ref_raw["mapp"]  # From reference file
    # gc = ref_raw["gc"]  # From reference file

    # Apply various filters based on thresholds
    # Filter samples based on coverage threshold
    sampfilter1 = (Y_raw.sum(axis=0) <= cov_thresh)
    print(f"Removed {sum(sampfilter1)} samples due to failed library preparation.")

    # Filter samples based on minimum coverage requirement
    sampfilter2 = (Y_raw.mean(axis=0) <= minCountQC)
    print(f"Removed {sum(sampfilter2)} samples due to failure to meet min coverage requirement.")

    # Filter samples based on proportion of mapped reads
    sampfilter3 = (QCmetric_raw["mapped_prop"] < mapq20_thresh)
    print(f"Removed {sum(sampfilter3)} samples due to low proportion of mapped reads.")

    # Exclude samples that failed any of the above filters
    filter_mask = sampfilter1 | sampfilter2 | sampfilter3
    if filter_mask.any():
        Y = Y_raw[:, ~filter_mask]
        sampname = [sampname_raw[i] for i, keep in enumerate(~filter_mask) if keep]
        QCmetric = QCmetric_raw.loc[~filter_mask, :]
    else:
        Y = Y_raw
        sampname = sampname_raw
        QCmetric = QCmetric_raw
    # print("Y1",Y)
    # Apply filters on bins based on GC content and mappability
    # Filter bins based on extreme GC content
    binfilter1 = (gc < gc_thresh[0]) | (gc > gc_thresh[1])
    print(f"Excluded {sum(binfilter1)} bins due to extreme GC content.")

    # Filter bins based on mappability
    binfilter2 = (mapp < mapp_thresh)
    print(f"Excluded {sum(binfilter2)} bins due to low mappability.")

    # Exclude bins that failed any of the above filters
    bin_filter_mask = binfilter1 | binfilter2
    if bin_filter_mask.any():
        ref = ref_raw.loc[~bin_filter_mask, :]
        Y = Y[~bin_filter_mask, :]
    else:
        ref = ref_raw
        Y = Y
    # print("Y2",Y)

    # Compute normalization factors
    # Calculate a normalization factor based on non-zero values
    # Y_nonzero = Y[Y.sum(axis=1) != 0, :]

    # Assuming Y is your NumPy array
    filtered_rows = np.apply_along_axis(lambda x: not np.any(x == 0), axis=1, arr=Y)
    Y_nonzero = Y[filtered_rows]


    print("Y.shape",Y.shape)
    print("Y_nonzero.shape",Y_nonzero.shape)
    if Y_nonzero.shape[0] <= 10:
        print("Adopt arithmetic mean instead of geometric mean")
        pseudo_sample = Y.mean(axis=1)
        N = (Y / pseudo_sample).median(axis=1)
    else:
        # pseudo_sample = Y_nonzero.apply(lambda x: np.exp(np.sum(np.log(x)) / len(x)), axis=1)
        pseudo_sample = np.array([np.exp(np.sum(np.log(row)) / len(row)) for row in Y_nonzero]).reshape(-1, 1)
        print("pseudo_sample.shape",pseudo_sample.shape)
        print(Y_nonzero / pseudo_sample)
        N = np.nanmedian(Y_nonzero / pseudo_sample,axis=0)
        # N = np.median(Y_nonzero / pseudo_sample[:, np.newaxis], axis=0)
        print("N dim", N.shape)

    # Filter samples based on library size calculation
    sampfilter4 = (N == 0)
    print(f"Removed {sum(sampfilter4)} samples due to excessive zero read counts in library size calculation.")
    if sampfilter4.any():
        Y = Y[:, ~sampfilter4]
        sampname = [sampname[i] for i, keep in enumerate(~sampfilter4) if keep]
        QCmetric = QCmetric.loc[~sampfilter4, :]
        N = N.loc[~sampfilter4]
    # print("Y3",Y)
    # Normalize and filter bins based on deviation from median
        # Normalize and filter bins based on deviation from median
    print("N shape:",N.shape)
    # print("N:",N)

    n_rows, n_cols = Y.shape
    Nmat = np.tile(N, (n_rows, 1))
    print("Nmat shape:",Nmat.shape)
    print("Nmat:",Nmat)

    # print(Y)
    bin_sum = (Y / Nmat).sum(axis=1)
    print(bin_sum)
    print(np.median(bin_sum))
    print(MAD(bin_sum))
    print(sum(bin_sum))
    print(nMAD)
    binfilter3 = (bin_sum >= (np.median(bin_sum) - nMAD * MAD(bin_sum))) & \
                 (bin_sum <= (np.median(bin_sum) + nMAD * MAD(bin_sum)))
    Y = Y[binfilter3, :]
    # print("Y4",Y)
    ref = ref.loc[binfilter3, :].reset_index(drop = True)

    print(f"Filtered {Y_raw.shape[0] - Y.shape[0]} bins based on deviation from median.")

    # Convert QC metric to a DataFrame and return results
    QCmetric = QCmetric.reset_index(drop=True)
    print(f"There are {Y.shape[1]} samples and {Y.shape[0]} bins after QC step.")

    result = {
        "Y": Y,
        "sampname": sampname,
        "ref": ref,
        "QCmetric": QCmetric
    }

    df = pd.DataFrame(Y, columns = sampname)
    # print(df)
    # print(ref)
    df1 = pd.concat([ref,df],axis = 1)
    df1.columns = ['chr','start','end'] + sampname

    return df1

def preprocess(bam_dirs, pre_bam_list, reference_file, outdir, reftype, seq_type,
               n_thread = 50 ,mapq_thresh = 0,resolution=500):
    if pre_bam_list is None:
        raw_bam_list = []
        for bam_dir in bam_dirs:
            raw_bam_list.extend(glob.glob(bam_dir+"/*.bam"))
    else:
        raw_bam_list = pre_bam_list


    samp_list = [ ]
    bam_list = []
    for file in raw_bam_list:
        samp = file.split('/')[-1].split(',')[0].split('_')[0] 
        if samp!='all':
            samp_list.append(samp)
            bam_list.append(file)
    

    # samp_list=samp_list[:20]
    # bam_list=bam_list[:20]
    print("num of samples:", len(samp_list))
    print(samp_list)
    # print(bam_list)
    # exit()



    if not os.path.exists(outdir):
        os.system("mkdir -p "+ outdir)

    bambed_dc = get_bam_bed(bam_list, samp_list, reference_fasta= reference_file ,hgref=reftype, resolution=resolution, sex=False)
    map_list = get_map(reftype)
    gc_list = get_GC(reftype)

    ref_list = []
    for _,df in bambed_dc['ref'].items():
        ref_list.append(df)
    ref_raw = pd.concat(ref_list)
    print(ref_raw)


    print(len(map_list))
    print(len(gc_list))


    print("num of bins:",ref_raw.shape[0])
    Y_raw = get_coverage_scDNA(bambed_dc, mapq_thresh, 
                                    seq = seq_type, hgref = reftype,
                                    n_thread=n_thread, outdir = outdir, n_bin= ref_raw.shape[0])


    print(Y_raw)
    Y_raw.to_csv(outdir+"/Y_raw.csv",index = False)


    QCmetric_raw = get_samp_QC(bambed_dc, n_thread)
    QCmetric_raw.to_csv(outdir + "/QCmetric_raw.csv", index = False)
    print(QCmetric_raw)

    QCmetric_raw = pd.read_csv(outdir+"/QCmetric_raw.csv")
    Y_raw = pd.read_csv(outdir+"/Y_raw.csv")
    Y_raw = Y_raw.values


    print(Y_raw.shape)
    print(QCmetric_raw.shape)

    # Example usage:
    # Replace Y_raw, sampname_raw, ref_raw, and QCmetric_raw with your data
    df = perform_qc(Y_raw, samp_list, ref_raw, QCmetric_raw, map_list, gc_list)
    df.to_csv(outdir+"/genome_cov.bed",index = False, sep = '\t')
    return 











if __name__ == "__main__":
    # Example usage:
    bam_dirs = ["/data/maiziezhou_lab/Datasets/singlecell_data/TNBC/DNA_bowtie2_rmdup/DNA3020_rmdup/",
                "/data/maiziezhou_lab/Datasets/singlecell_data/TNBC/DNA_bowtie2_rmdup/DNA3022_rmdup/"]

    n_thread = 50
    reference_file = "/data/maiziezhou_lab/Softwares/refdata-GRCh38-2.1.0/fasta/genome.fa"
    reftype = "hg38"
    seq_type = 'paired-end'
    mapq_thresh = 0
    print("min mapq:", mapq_thresh)
    resolution = 500
    # outfile = "test.csv"
    outdir = "./out_v3/"
    preprocess(bam_dirs, reference_file, outdir, reftype, seq_type,
               n_thread = 50 ,mapq_thresh = 0,resolution=500)





















