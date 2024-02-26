# Orginal Vignettes
The package paths are:

```bash
https://bioconda.github.io/recipes/chisel/README.html

https://github.com/raphael-group/chisel
```

# Environments
CHISEL is a Python based method. I installed `conda` and set up the package environment in my path. You can use that environment by adding the code chunk below in your bash before using any functions in CHISEL package:

```bash
source /data/maiziezhou_lab/Zihang/conda/bin/activate
conda activate chisel
```

# For scDNA-seq data
### Prepare the data for CHISEL
`chisel` function uses bam files containing CB:Z:` tags and 'Phased SNPs file' in tsv format as their input. So we need to prepare the data for that.


The barcode bam files should contain the `CB:Z:` tags. If the data do not contain `CB:Z:` tag, one option is to generate `RG:Z:` tag and then use `chisel_prep` to continue to generate `CB:Z:` tag.

Code chunk for generation of `RG:Z:` in tumor cells (normal cells are similar):

```bash
ml picard/2.18.27
patient='3020'
type='DNA'$patient
base='/data/maiziezhou_lab/Zihang/CNV/'
database=$base'dataset/TNBC/'
file=$database'DNA_bowtie2_rmdup_bcftools_phasing/KTN'$patient'.txt'
while read line; do
    echo $line
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
        I=$database'DNA_bowtie2_rmdup/'$type'_rmdup/'$line'_grch38_rmdup.bam'  \
        O=$database'DNA_bowtie2_rmdup/'$type'_rmdup_rg/'$line'_grch38_rmdup_rg.bam'  \
        RGID=4 \
        RGLB=lib1 \
        RGPL=illumina \
        RGPU=unit1 \
        RGSM=$line
done < $file
```

Code chunk for generation of `CB:Z:` in tumor cells (normal cells are similar):

```bash
patient='3020'
type='DNA'$patient
base='/data/maiziezhou_lab/Zihang/CNV/'
database=$base'dataset/TNBC/'
file=$database'DNA_bowtie2_rmdup_bcftools_phasing/KTN'$patient'.txt'
while read line; do
    echo $line
    chisel_prep -o $database'DNA_bowtie2_rmdup/'$type'_rmdup_rg_cb/'$line'_grch38_rmdup_rg_cb.bam' \
        -x $database'DNA_bowtie2_rmdup/'$type'_rmdup_rg_cb/' \
        $database'DNA_bowtie2_rmdup/'$type'_rmdup_rg/'$line'_grch38_rmdup_rg.bam'
done < $file
```
After this step, we use `samtool merge` to merge individual bam file into one bam file:

```bash
source /data/maiziezhou_lab/Zihang/CNV/CHISEL/setup
cd $database'DNA_bowtie2_rmdup/'$type'_rmdup_rg/
samtools merge -o all_tumor_cells.bam -@ 15 *.bam
```
Code chunk for generate the .bcf and .vcf file given the merged bam file:

```bash
source /data/maiziezhou_lab/Zihang/CNV/CHISEL/setup

bcftools mpileup -f /data/maiziezhou_lab/Zihang/CNV/index/human_ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb/all_tumor_cells.bam | \
    bcftools call -mv -Ob -o /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf/all_tumor_cells_grch38.bcf

bcftools view /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf/all_tumor_cells_grch38.bcf > \
    /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf/all_tumor_cells_grch38.vcf
```
The next move is to generate indexed bam file to run `whatshap`.

```bash
source /data/maiziezhou_lab/Zihang/CNV/CHISEL/setup
samtools index all_tumor_cells.bam 
```

The `chisel` function requests the user to merge all the single-cell bam files into one single bam file. So you can merge the tumor cells' bam file into one bam file and merge the normal cells' bam file into one bam file, too. Then I used `whatshap` to generate the 'Phased SNPs file'.

```bash
source /data/maiziezhou_lab/Zihang/conda/bin/activate
conda activate whatshap

whatshap phase -o /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf_whatshap/all_tumor_cells_grch38_phased.vcf \
    --reference=/data/maiziezhou_lab/Zihang/CNV/index/human_ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf/all_tumor_cells_grch38.vcf\
    /data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb/all_tumor_cells.bam \
    --ignore-read-groups
```

Then I used `Python` to convert the vcf file into tsv file:

```python
import os
import numpy as np
import csv

patient = '3020'
nucleic_type = 'DNA' + patient
base = '/data/maiziezhou_lab/Zihang/CNV/'
database = base + 'dataset/TNBC/'

file_vcf = database + 'DNA_bowtie2_rmdup_bcftools_phasing/' + nucleic_type + '_vcf_whatshap/all_tumor_cells_grch38_phased.vcf'
file_tsv = database + 'DNA_bowtie2_rmdup_bcftools_phasing/' + nucleic_type + '_vcf_whatshap_tsv/all_tumor_cells_grch38_phased.tsv'
out_array = []
# count=0
chr_list = []
for i in range(22):
    chr_list.append('chr'+str(i+1))

with open(file_vcf, "r") as fr:
    for line_vcf in fr:
        # count = count + 1
        line_vcf_all = line_vcf.split('\t')
        if line_vcf_all[0] not in chr_list:
            continue
        # print(line_vcf)
        line_vcf_1 = line_vcf_all[0]
        line_vcf_2 = line_vcf_all[1]
        line_vcf_3 = line_vcf_all[9].split(':')[0]
        if line_vcf_3 in ['1|0', '0|1']:
            out_array.append([line_vcf_1, line_vcf_2, line_vcf_3])
with open(file_tsv, "w") as fw:
    tsv_write = csv.writer(fw, delimiter='\t')
    tsv_write.writerows(out_array)
```


### Use CHISEL to infer CNV
After preparation, we can use function `chisel` to infer CNV:

```bash
patient='3020'
base='/data/maiziezhou_lab/Zihang/CNV/'
database=$base'dataset/TNBC/'
file=$database'DNA_bowtie2_rmdup_bcftools_phasing/KTN'$patient'.txt'
mkdir $base'CHISEL/CHISEL_bowtie2_rmdup_bcftools_whatshap/all_50kb'
# cd $base'CHISEL/CHISEL_bowtie2_rmdup_bcftools_whatshap/all_50kb'
chisel -x $base'CHISEL/TNBC_bowtie2_rmdup_bcftools_whatshap/all_50kb' \
    -t $database'DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb/all_tumor_cells.bam' \
    -n $database'DNA_bowtie2_rmdup/DNA3022_rmdup_rg_cb/all_tumor_cells.bam' \
    -r $base'index/human_ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna' \
    -l $database'DNA_bowtie2_rmdup_bcftools_phasing/DNA3020_vcf_whatshap_tsv/all_tumor_cells_grch38_phased.tsv' \
    -b 50kb --seed 1 -j 15
```

Note that, please make sure the genome reference is the same during the whole procedure.
