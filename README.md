# CNVeil

## Table of content
- [CNVeil](#CNVeil)
  - [Install CNVeil](#install-cnveil)
  - [Run CNVeil](#Run-CNVeil)
    - [Step 1: Extract read depth matrix](#Step-1-Extract-read-depth-matrix)
    - [Step 2: Infer total copy number](#Step-2-Infer-total-copy-number)
    - [Step 3: Infer allele-specific copy number](#Step-3-Infer-allele-specific-copy-number)
    - [Step 4: Infer haplotype-specific copy number](#Step-3-Infer-allele-specific-copy-number)
- [Benchmark](#Benchmark)
  - [Prepare scDNA-seq data](#prepare-scDNA-seq-data)
  - [Run Tools](#run-tools)

## CNVeil
### Install CNVeil

You can install CNVeil through github
```
git clone git@github.com:maiziezhoulab/CNVeil.git
```
To set up the environment, you should have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed on your machine. Then simply run
```
# create conda env CNVeil_py3
conda env create -f CNVeil/env_yamls/environment_py3.yml 
# create conda env CNVeil_py2
conda env create -f CNVeil/env_yamls/environment_py2.yml
# create conda env CNVeil_R451
conda env create -f CNVeil/env_yamls/environment_R451.yml
```

You will have 3 virtual environments installed (CNVeil_py3, CNVeil_py2 and CNVeil_R451). And they will be used in diferent steps. You will need to load the corresponding virtual environment before you run each step in CNVeil pipeline.

### Run CNVeil
#### Step 1 Extract read depth matrix

This step bins reads from the input BAM files, performs GC correction, and reformats the corrected bin-level read depth into a matrix for downstream CNVeil analysis.

Script:

```bash
CNVeil_extract_read_depth.py
```

Required conda environment:

```bash
conda activate CNVeil_R451
```

### Arguments

| Argument | Short | Required | Description |
|---|---|---|---|
| `--merged_bam_file` | `-mgbam` | Yes | Path to the merged BAM file containing reads from all cells. This is used for binning and GC correction. |
| `--split_bam_dir` | `-spbam` | Yes | Directory containing split BAM files for individual cells. |
| `--work_dir` | `-w` | Yes | Working directory for CNVeil output. The script will create a `Total_CN/` subdirectory inside the working folder. |
| `--reftype` | `-rt` | Yes | Reference genome version. Must be either `hg19` or `hg38`. |
| `--num_threads` | `-t` | No | Number of threads to use. Default: `20`. |
| `--cell_node` | `-cn` | No | Optional cell-node mapping file. This is mainly designed for simulation data generated from SimscTree. If provided, cell IDs will be transformed to node IDs in the final matrix. For real data, this is not needed. |

### Example usage

For real data:

```bash
python CNVeil_extract_read_depth.py \
  --merged_bam_file /path/to/merged.bam \
  --split_bam_dir /path/to/split_bams \
  --work_dir /path/to/work_dir \
  --reftype hg38 \
  --num_threads 20
```

For simulation data with a cell-node mapping file:

```bash
python CNVeil_extract_read_depth.py \
  --merged_bam_file /path/to/merged.bam \
  --split_bam_dir /path/to/split_bams \
  --work_dir /path/to/work_dir \
  --reftype hg38 \
  --num_threads 20 \
  --cell_node /path/to/cell_node.tsv
```

### Output files

The script creates the following output directory:

```text
<work_dir>/Total_CN/
```

Main output:

```text
<work_dir>/Total_CN/RDmatrix.csv
```

This is the final read-depth matrix used by downstream CNVeil steps.


#### Step 2 Infer Total Copy Number

This step takes the read-depth matrix generated from Step 1, and outputs a total copy-number matrix for downstream CNVeil analysis.

Script:

```bash
CNVeil_total_CN.py
```

Required conda environment:

```bash
conda activate CNVeil_py3
```


This script expects the Step 1 output directory structure:

```text
<work_dir>/Total_CN/
```

Required input file:

```text
<work_dir>/Total_CN/RDmatrix.csv
```



##### Arguments

| Argument | Short | Required | Description |
|---|---|---|---|
| `--work_dir` | `-w` | Yes | Working directory for CNVeil. The script reads input from `<work_dir>/Total_CN/RDmatrix.csv` and writes outputs to `<work_dir>/Total_CN/`. |
| `--reftype` | `-rt` | Yes | Reference genome version. Must be either `hg19` or `hg38`. |
| `--seq_type` | `-st` | Yes | Sequencing type. Must be either `paired-end` or `single-end`. This argument is currently parsed but not directly used in the script. |
| `--data_type` | `-dtype` | Yes | Data type. Must be either `10x` or `other`. This affects QC and low-coverage bin filtering behavior. |
| `--sample` | `-sp` | No | Sample name used in output file names. Default: `sample`. |
| `--num_threads` | `-t` | No | Number of threads. Default: `20`. This argument is currently parsed but not directly used in the main CN inference steps. |
| `--cell_node` | `-cn` | No | Optional cell-node mapping file. Mainly designed for simulation data. If provided, cell names are transformed to node labels in the final CNV output. |


##### Example usage

For real 10x data:

```bash
python CNVeil_total_CN.py \
  --work_dir /path/to/work_dir \
  --reftype hg38 \
  --seq_type <paired-end or single-end> \
  --data_type 10x \
  --sample SAMPLE_NAME \
  --num_threads 20
```

For other single-cell sequencing data (non-10x):

```bash
python CNVeil_total_CN.py \
  --work_dir /path/to/work_dir \
  --reftype hg38 \
  --seq_type <paired-end or single-end> \
  --data_type other \
  --sample SAMPLE_NAME
```

For simulation data with a cell-node mapping file:

```bash
python CNVeil_total_CN.py \
  --work_dir /path/to/work_dir \
  --reftype hg38 \
  --seq_type <paired-end or single-end> \
  --data_type 10x \
  --sample SAMPLE_NAME \
  --cell_node /path/to/cell_node.tsv
```



##### Main output files

All outputs are written under:

```text
<work_dir>/Total_CN/
```

Main total copy-number matrix:

```text
<work_dir>/Total_CN/CNV_<sample>.csv
```

This file has cells as rows and genomic bins as columns.

Transposed copy-number matrix:

```text
<work_dir>/Total_CN/CNV_<sample>.csv.T
```

This file has genomic bins as rows and cells as columns.


#### Step 3: Infer allele-specific copy number

This step uses the total copy-number results from Step 2 together with SNP allele-count information to infer allele-specific copy number (ASCN).


Script:

```bash
CNVeil_allele_CN.py
```

Required conda environment:

```bash
conda activate CNVeil_py3
```

##### Required input from previous steps

This step expects the CNVeil working directory to already contain the total CN result from Step 2:

```text
<work_dir>/Total_CN/CNV_<sample>.csv
```

It also requires a directory containing SNP/RDS-derived input files:

```text
<var_dir>/
```

The SNP loader expects this directory to contain:

```text
var_all.rds
ref_all.rds
alt_all.rds
```

These files are used to generate per-cell allele-count files for allele-specific CN inference. To get allele-count information, check [here](preprocess/allele_count/README.md).

#### Arguments

| Argument | Short | Required | Description |
|---|---|---|---|
| `--work_dir` | | Yes | Working directory for the CNVeil pipeline. The script reads total CN results from `<work_dir>/Total_CN/` and writes allele-specific CN results to `<work_dir>/Allele_CN/`. |
| `--var_dir` | | Yes | Directory containing SNP/RDS-derived files, including `var_all.rds`, `ref_all.rds`, and `alt_all.rds`. |
| `--ref_type` | | Yes | Reference genome type. Must be `hg19` or `hg38`. |
| `--snp_source` | | Yes | SNP source type. Must be `tumor` or `tumor-normal`.|
| `--convert` | | No | Optional cell name/barcode converter file. Use this if cell names in allele-count files and total CN files need to be matched. |
| `--snp_file` | | No | Optional SNP genotype file. Usually only needed when external genotype information is available. |
| `--num_threads` | `-t` | No | Number of parallel threads. Default: `22`. |

##### Example usage

Basic usage:

```bash
python CNVeil_allele_CN.py \
  --work_dir /path/to/work_dir \
  --var_dir /path/to/var_rds_dir \
  --ref_type hg38 \
  --snp_source tumor-normal \
  --num_threads 22
```

With a cell-name converter:

```bash
python CNVeil_allele_CN.py \
  --work_dir /path/to/work_dir \
  --var_dir /path/to/var_rds_dir \
  --ref_type hg38 \
  --convert /path/to/collected_barcodes_SAMPLE.txt \
  --snp_source tumor-normal \
  --num_threads 22
```

With an optional SNP genotype file:

```bash
python CNVeil_allele_CN.py \
  --work_dir /path/to/work_dir \
  --var_dir /path/to/var_rds_dir \
  --ref_type hg38 \
  --convert /path/to/collected_barcodes_SAMPLE.txt \
  --snp_source tumor \
  --snp_file /path/to/snp_genotype.tsv \
  --num_threads 22
```

##### Main output files

The script writes allele-specific CN results under:

```text
<work_dir>/Allele_CN/
```

Important outputs include:

```text
<work_dir>/Allele_CN/ascn_seg.csv
```

Allele-specific copy number inferred at the cell/segment level.


#### Step 4 Infer Haplotype-specific CN


## Benchmark
### Prepare scDNA-seq data

Here is an example of generating the bam files from fastq files downloaded from NCBI cell by cell. The patient is named as 'KTN126'.

```
patient="KTN126"
file=YOUR_PATH_FOR_SRR_LIST
srapath=YOUR_PATH_TO_STORE_SRA_FILE
fastqpath=YOUR_PATH_TO_STORE_FASTQ_FILE
mkdir ${fastqpath}

 while read line; do
     echo $line
     YOUR_PATH_TO_SRATOOLKIT/bin/fastq-dump -I --split-files $srapath$line"/"$line".sra" 
     -O $fastqpath$line
 done < $file
```
Generate the fastq files from the. sra file with sra toolkit.  
```
ml Intel/2017.4.196 BWA/0.7.17
ml GCC/11.3.0 SAMtools/1.18
sampath=${srapath}/sampath_SraAccList${patient}/
align_dir=${srapath}/sampath_SraAccList${patient}/
readgroup_path=${srapath}/bampath_SraAccList${patient}_readgroup/
dedup_path=${srapath}/bampath_SraAccList${patient}_dedup/

mkdir ${sampath} ${readgroup_path} ${dedup_path}
 
while read line; do
    echo $line
    cd $align_dir 
     bwa mem -M -t 8 YOUR_PATH_TO_REF_GENOME_FILE $fastqpath/$line/"$line"_1.fastq $fastqpath/$line/"$line"_2.fastq > $align_dir/"$line".sam
     samtools view -bS $align_dir/"$line".sam > $align_dir/"$line".bam 
     java -jar YOUR_PATH_TO_PICARD/picard.jar SortSam -I $align_dir/"$line".bam -SORT_ORDER coordinate -O $align_dir/"$line"_hg38.sorted.bam 
     java -jar YOUR_PATH_TO_PICARD/picard.jar AddOrReplaceReadGroups -I $align_dir/"$line"_hg38.sorted.bam -O $readgroup_path/"$line"_hg38.sorted.rg.bam -RGID $readgroup_path/"$line"_hg38/ -RGLB NAVIN_Et_Al -RGPL ILLUMINA -RGPU machine -RGSM $readgroup_path/"$line"_hg38/
     java -jar YOUR_PATH_TO_PICARD/picard.jar MarkDuplicates -REMOVE_DUPLICATES true -I $readgroup_path/"$line"_hg38.sorted.rg.bam -O $dedup_path/"$line"_hg38.sorted.rg.dedup.bam -METRICS_FILE $dedup_path/"$line"_hg38.sorted.rg.dedup.metrics.txt -PROGRAM_RECORD_ID MarkDuplicates -PROGRAM_GROUP_VERSION null -PROGRAM_GROUP_NAME MarkDuplicates
     java -jar YOUR_PATH_TO_PICARD/picard.jar BuildBamIndex -I $dedup_path/"$line"_hg38.sorted.rg.dedup.bam -O $dedup_path/"$line"_hg38.sorted.rg.dedup.bai
done < $file
```
Generate the bam files after alignment.

Next, we generate the bed files with hg19/hg38 reference.
```
ml GCC/8.2.0 BEDTools/2.28.0
bedpath_hg38=${srapath}/bedpath_SraAccList${patient}_hg38
bedpath_hg19=${srapath}/bedpath_SraAccList${patient}_hg19
mkdir ${bedpath_hg38} ${bedpath_hg19}
while read line; do
    echo $line
    bamToBed -i $dedup_path/$line'_hg38.sorted.rg.dedup.bam' > \
        ${bedpath_hg38}/$line'_grch38_rmdup.bed'
done < $file

while read line; do
    echo $line
    YOUR_PATH_TO_LIFTOVER/liftOver \
        ${bedpath_hg38}/$line'_grch38_rmdup.bed' \
        'YOUR_PATH_TO_LIFTOVER_CHAINFILE/liftOver_ChainFiles/hg38ToHg19.over.chain' \
        ${bedpath_hg19}/$line'_hg19_rmdup.bed' \
        ${bedpath_hg19}/$line'_hg19_unlifted_rmdup.bed'
done < $file
```
Zipped the bed file, which is recommended for the Ginkgo web platform to analyze.
```
cd ${bedpath_hg19}
```















