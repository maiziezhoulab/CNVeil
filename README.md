# CNVeil

## Table of content
- [CNVeil-Python](#CNVeil-Python)
  - [Install python version CNVeil](#install-python-version-cnveil)
  - [Run CNVeil-Python](#Run-CNVeil-Python)
- [CNVeil-R](#CNVeil-R)
  - [Install R version CNVeil](#install-R-version-cnveil)
  - [Run CNVeil-R](#Run-CNVeil-R)
- Benchmark
  - [Prepare scDNA-seq data](#prepare-scDNA-seq-data)
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
output_path=${srapath}/'bedpath_SraAccList'${patient}'_zipped'
mkdir ${output_path}
for file in *"_hg19_rmdup.bed"; do
    gzip -c "$file" > "${output_path}/${file%.bed}.bed.gz"
```
