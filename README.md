# CNVeil

## Table of content
- [CNVeil-Python](#CNVeil-Python)
  - [Install python version CNVeil](#install-python-version-cnveil)
  - [Run CNVeil-Python](#Run-CNVeil-Python)
- [CNVeil-R](#CNVeil-R)
  - [Install R version CNVeil](#install-R-version-cnveil)
  - [Run CNVeil-R](#Run-CNVeil-R)
- Benchmarking
  - [Prepare scDNA-seq data](#prepare-scDNA-seq-data)
  - [Run Tools]

## CNVeil-Python
### Install python version CNVeil

You can install python version CNVeil through github
```
git clone git@github.com:maiziezhoulab/CNVeil.git
```
To set up the environment, you should have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed on your machine. Then simple run
```
conda env create -f CNVeil/Python/environment.yml
```
You will have a virtual environment called CNVeil installed. You must load this virtual environment before you run the CNVeil-python pipeline.
```
conda activate CNVeil
```

### Run CNVeil-Python

The main code for CNVeil-python is 'CNVeil/python/CNVeil.py`. The input arguments for this code are listed below:

```
  --bam_dirs BAM_DIRS [BAM_DIRS ...], -bam BAM_DIRS [BAM_DIRS ...]
                        Input directory that contains BAM files; if you want to provide multiple directories, separate them by blank
                        space (default: None)
  --pre_bam_list PRE_BAM_LIST [PRE_BAM_LIST ...], -baml PRE_BAM_LIST [PRE_BAM_LIST ...]
                        either bam_dirs or pre_bam_list need to be provided; provide multiple bam files, separate them by blank space
                        (default: None)
  --reference REFERENCE, -ref REFERENCE reference genome path
  --reftype {hg19,hg38}, -rt {hg19,hg38}
  --seq_type {paired-end,single-end}, -st {paired-end,single-end}
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
  --prefix PREFIX, -px PREFIX the output file's prefix, default = Sample
  --n_thread N_THREAD, -t N_THREAD  number of threads the pipeline will be executed under
  --cell_node CELL_NODE, -cn CELL_NODE
                        cell node file (every line should be cell name , tab and its corresponding node name); optional, if not given, will not transform cell to node in the end (default: None)
  --change_rate_mode {q,h}, -crm {q,h}
                        q -> use 0.84 quantile for change rate cut-off; h -> use 0.1 as change rate cut-off (default: q)
  --ploidy_mode {gl,cl}, -pm {gl,cl}
                        estimate ploidy globally(gl) or by cluster(cl) (default: cl)

```

An example command to run the pipeline is provided below:

```
python3 CNVeil/Python/CNVeil.py \
-bam <bam_dir> \
-ref <refernce_file> \
-rt <hg19 or hg38> \
-st <paired-end or single-end> \
-o CNVeil_output \
-t 10 \
-px <sample_name> \
-crm q -pm cl
```
After the pipeline is finished, you will see the final CNV matrix in CNVeil_output/CNV_<sample_name>.csv.

## CNVeil-R
### Install R version CNVeil

You can install R version CNVeil by loading all the necessary packages described below first, make sure the R version is above 4.2.1:
```
optparse/1.7.3
data.table/1.14.2
ggplot2/3.4.4
cluster/2.1.3
stats/4.2.1
dplyr/1.0.9
matrixStats/1.2.0
kernlab/0.9.31
pracma/2.3.8
gtools/3.9.5
DescTools/0.99.45
GenomicRanges1.50.2
Rsamtools/2.14.0
BSgenome.Hsapiens.UCSC.hg38/1.4.5
SCOPE/1.10.0
WGSmapp/1.10.0
```
### Run CNVeil-R

To establish read depth matrix from bam file, please run the preprocess.r first. The input arguments for this code are listed below:

```
1. work_path: Path to the working directory.
2. ref: Type of the reference.
3. bam_path: Directory containing BAM files.
4. pattern: Pattern to match the BAM files (if needed).
5. bin_size: Size of the bins for the read depth matrix.
```
Run the script as follows:
```
Rscript preprocess.r [work_path] [ref] [bam_path] [pattern] [bin_size] 
```
The major result for this step is named as "genome_cov.bed". 

To correct the bias by normal cell, please use normal_cell_correction.r for to obtain RDmatrix.rds. 
Run the script as follows:
```
Rscript normal_cell_correction.r /path/to/work/directory /path/to/normal/cell/file(optional) /path/to/save/output/
```

Finally, we could call CNVs from the standardized matrix. Please run the call_cn.r. The input arguments for this code are listed below:

```
Rscript call_cn.r --input_dir /path/to/input --output_dir /path/to/output --prefix Sample --cell_node /path/to/cell_node.txt --change_rate_mode q --ploidy_mode gl
```


