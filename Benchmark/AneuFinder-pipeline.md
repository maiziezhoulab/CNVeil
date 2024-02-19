# Ginkgo orginal Vignettes
The package paths are:

```bash
https://bioconductor.org/packages/release/bioc/html/AneuFinder.html

https://github.com/ataudt/aneufinder
```

The path of vignettes written by the authors is:

```bash
http://bioconductor.org/packages/release/bioc/vignettes/AneuFinder/inst/doc/AneuFinder.pdf
```

# For our data
### Mappability correction
The first step in the workflow is to produce a reference file for mappability correction.  It's also possible to proceed without mappability correction, in which case the algorithm defaults to using fixed-width bins. 

### Blacklisting
To further improve the quality of the results and remove artifacts caused by high mappability repeat regions, e.g. near centromers, a blacklist can be used with option `Aneufinder(..., blacklist)`. All reads falling into the regions specified by the blacklist will be discarded when importing the read files. You can either download a blacklist from the UCSC genome browser, e.g. the “DAC Blacklisted Regions from ENCODE/DAC(Kundaje)” mappability track, or make your own. For optimal results, we advice to make your own blacklist from a euploid reference. The code chunk below is an example of how I used to generate the blacklist file:

```r
library(AneuFinder)
bedfile <- "/data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_bowtie2_rmdup/DNA3020_rmdup_rg_cb/all_tumor_cells.bam"
bins <- binReads(bedfile, assembly='hg38', binsizes=500e3,
                 chromosomes=paste0("chr",c(1:22)))[[1]]
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.999)

blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist <- reduce(blacklist)
blacklist.file <- "/data/maiziezhou_lab/Zihang/CNV/AneuFinder/TNBC_bowtie2_rmdup/blacklist/bl3020"
exportGRanges(blacklist, filename=blacklist.file, header=FALSE,
              chromosome.format='UCSC')
```
Note that, in this data chunk:

- `bedfile` should be a bam file containing all the cells you try to infer CNV, you can use `samtools merge` if you have already had individual cells

- `assembly` should be changed to `hg19` or so if you generate your data using hg19 reference.

- `binsizes` should be set to the specific bin size you used.

-  Bins with read count above and below the 0.999 and 0.05 quantile are taken as blacklist

- `blacklist.file` is the path of the result. Note that when running the example, the name of the output file would be `bl3020.bed.gz`

- `chromosome.format` Please notice that this is quite important. If you used the ref from NCBI, you should replace `UCSC` by `NCBI` to avoid errors.

### GC-content correction
Just set the parameters in function `Aneufinder` like:

```r
correction.method = 'GC'
GC.BSgenome = BSgenome.Hsapiens.UCSC.hg38
```
Note that, in this data chunk:

- `GC.BSgenome` should be changed to `BSgenome.Hsapiens.UCSC.hg19` or so if you generate your data using hg19 reference.

### Run Aneufinder

```r
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg38)

basepath = '/data/maiziezhou_lab/Zihang/CNV/'
var.width.ref <- NULL
blacklist <- paste0(basepath, 'AneuFinder/TNBC_bowtie2_rmdup/blacklist/bl3020.bed.gz')
datafolder <- paste0(basepath, 'dataset/TNBC/DNA_bowtie2_rmdup/DNA3020_rmdup_bed')
outputfolder <- paste0(basepath, 'AneuFinder/TNBC_bowtie2_rmdup/500kb/nomapcor_paired')

## Library for GC correction

## Produce output files
Aneufinder(inputfolder = datafolder, outputfolder = outputfolder,
           pairedEndReads = TRUE, assembly = 'hg38', numCPU = 4,
           binsizes = 500e3, variable.width.reference = var.width.ref,
           chromosomes = paste0("chr",c(1:22)), blacklist = blacklist,
           correction.method = 'GC', GC.BSgenome = BSgenome.Hsapiens.UCSC.hg38,
           refine.breakpoints=FALSE,
           method = 'edivisive')
```
Note that, in this data chunk:

- `blacklist` is the path of the blacklist file generated before

- `datafolder` is a folder contains all the `bed` or `bam` files of the individual cells exclusively.

- `outputfolder` should be the place where you put your results

- `assembly` should be changed to `hg19` or so if you generate your data using hg19 reference.

- `binsizes` should be set to the specific bin size you used.

- `GC.BSgenome` should be changed to `BSgenome.Hsapiens.UCSC.hg19` or so if you generate your data using hg19 reference.






