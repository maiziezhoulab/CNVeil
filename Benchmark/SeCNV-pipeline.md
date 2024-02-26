# Orginal Vignettes
The package paths are:

```bash
https://github.com/raphael-group/chisel](https://github.com/deepomicslab/SeCNV)
```

# Environments
SeCNV is a Python based method. I installed `conda` and set up the package environment in my path. You can use that environment by adding the code chunk below in your bash before using any functions in SeCNV package:

```bash
ml Anaconda3/5.0.1
source activate /home/liy109/.conda/envs/SeCNV

ml GCC/10.2.0
ml SAMtools/1.12
```

# For our data
### Prepare the data for SeCNV
To run SeCNV, the bigwig files, hg19_mappability.bigWig and hg38_mappability.bigWig, should be downloaded from (https://drive.google.com/drive/folders/1XGuXa9muRtOiAjtz1J4gWO5Qk3c5-PwS) and put under the Script folder. And picard.jar file should also under the dataset folder.

Code chunk for generation of Input datasets (Using T10 dataset with hg19 reference file as an example):

```bash
module load Anaconda3/5.0.1
source activate /home/liy109/.conda/envs/SeCNV

ml Intel/2017.4.196 BWA/0.7.17
ml GCC/10.2.0 SAMtools/1.12

name="SraAccListT10"
file=/data/maiziezhou_lab/Datasets/singlecell_data/Navin/raw_data/SRRlist_DNA/"$name".txt
fastqpath=/data/maiziezhou_lab/Datasets/singlecell_data/Navin/raw_data/fastqpath_"$name"
sampath=/data/maiziezhou_lab/Datasets/singlecell_data/Navin/raw_data/samBWA_"$name"/
align_dir=/data/maiziezhou_lab/Datasets/singlecell_data/SeCNV/T10_hg19

while read line; do
    echo $line
    cd $align_dir 
     bwa mem -M -t 8 /data/maiziezhou_lab/yikang/SeCNV/Scripts/hg38.fa $fastqpath/$line/"$line"_1.fastq > $align_dir/T10/"$line".sam
     samtools view -bS $align_dir/T10/"$line".sam > $align_dir/T10/"$line".bam 
     java -jar picard.jar SortSam -I $align_dir/T10/"$line".bam -SORT_ORDER coordinate -O $align_dir/Input/"$line"_hg38.sorted.bam 
     java -jar picard.jar AddOrReplaceReadGroups -I $align_dir/Input/"$line"_hg38.sorted.bam -O $align_dir/Input/"$line"_hg38.sorted.rg.bam -RGID $align_dir/Input/"$line"_hg38/ -RGLB NAVIN_Et_Al -RGPL ILLUMINA -RGPU machine -RGSM $align_dir/Input/"$line"_hg38/
     java -jar picard.jar MarkDuplicates -REMOVE_DUPLICATES true -I $align_dir/Input/"$line"_hg38.sorted.rg.bam -O $align_dir/Input/"$line"_hg38.sorted.rg.dedup.bam -METRICS_FILE $align_dir/Input/"$line"_hg38.sorted.rg.dedup.metrics.txt -PROGRAM_RECORD_ID MarkDuplicates -PROGRAM_GROUP_VERSION null -PROGRAM_GROUP_NAME MarkDuplicates
     java -jar picard.jar BuildBamIndex -I $align_dir/Input/"$line"_hg38.sorted.rg.dedup.bam -O $align_dir/Input/"$line"_hg38.sorted.rg.dedup.bai
done < $file
```

For other datasets, e.g. T16 and TNBC datasets with reference hg38. We only have to change the `name` and `data path` of the previous code chunk.

### Use SeCNV to infer CNV
After preparation, we can use tool `SeCNV` to infer CNV:

```bash
module load Anaconda3/5.0.1
source /home/liy109/.conda/envs/SeCNV

ml GCC/10.2.0
ml SAMtools/1.12

cd Scripts
#python SeCNV.py input_folder output_folder reference_file
python SeCNV.py /data/maiziezhou_lab/Datasets/singlecell_data/SeCNV/T10_hg19/Input /data/maiziezhou_lab/Datasets/singlecell_data/SeCNV/T10_hg19/output hg19.fa
```

Note that, please make sure the genome reference is the same during the whole procedure. For different datasets, we only have to change the path of input folder, output folder and reference file.
