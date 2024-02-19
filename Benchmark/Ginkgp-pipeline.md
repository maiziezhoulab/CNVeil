# Statement
Ginkgo is a 'web-based' tool which you may only use online. The website is

```bash
http://qb.cshl.edu/ginkgo/
```

The code path is:

```bash
https://github.com/robertaboukhalil/ginkgo
```

Ginkgo only has a limited number of references which means that it only has hg19 other than hg38. Applying the alignment using grch38 as reference genome, I chose to convert the bam files to hg19-ref-based bam files and then performed Ginkgo. After that, I converted the Ginkgo result to grch38-ref-based bam files.

# For scDNA-seq data
### Prepare the data for Ginkgo (grch38 -> hg19)
Firstly, generate the bed files from bam files and use `liftOver` to convert the grch38-ref-based bed files to hg19-ref-based bed files. Example code chunck is provided here[Benckmark](https://github.com/maiziezhoulab/CNVeil/blob/main/README.md).

### Use Ginkgo to infer CNV
Upload the hg19-ref-based bed files to the Ginkgo website. Choose `hg19` as reference and choose bin-based before running.

Note that, each time you open the Ginkgo website, they will give a new website link for the new task. So, please keep the website link for further use, such as redirecting to the result page.

After running, you can copy the CNV results into a txt file. I named my file as `DNA3020.txt`. Then I extracted the first three lines and generated a new bed file named `hg19_bin500kb_cor.bed`. I also added the line number in to that file as the fourth column for 'alignment'. I used `liftOver` again similarly to convert the hg19-ref to grch38-ref and named that `hg38_bin500kb_cor.bed`:

```bash
/data/maiziezhou_lab/Zihang/CNV/Ginkgo/TNBC_bowtie2_rmdup/500kb/hg19:
 liftOver hg19_bin500kb_CNV.bed \
   $base'index/liftOver_ChainFiles/hg19ToHg38.over.chain' \
   hg38_bin500kb_cor.bed \
   hg38_unlifted_bin500kb_cor.bed
```

All the other tools we tried to compare with has the bins with beginning as x00001 and ending as (x+5)00000. So, the lifted grch38-file should be 500000-based smoothed. I gave an example in R and named the smoothed file `hg38_bin500kb_cor_shifted.bed`:

```r
hg38[,2] = round(hg38[,2]/500e3)*500e3 + 1
hg38[,3] = round(hg38[,3]/500e3)*500e3
```
