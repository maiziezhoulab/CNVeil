# CNVeil

## Table of content
- [CNVeil-Python](#CNVeil-Python)
  - [Install python version CNVeil](#install-python-version-cnveil)
  - [Run CNVeil-Python](#Run-CNVeil-Python)
- CNVeil-R(#CNVeil-R)
  - [Install R version CNVeil](#install-R-version-cnveil)
  - [Run CNVeil-R](#Run-CNVeil-R)
- Benchmarking
  - Prepare scDNA-seq data(#prepare-scDNA-seq-data)
  - Run Tools

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




