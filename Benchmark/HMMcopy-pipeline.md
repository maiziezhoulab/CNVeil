# Orginal Vignettes
The package paths are:

```bash
https://bioconductor.org/packages/release/bioc/html/HMMcopy.html

https://github.com/jujubix/HMMcopy
```

The path of vignettes written by the authors is:

```bash
https://bioconductor.org/packages/devel/bioc/vignettes/HMMcopy/inst/doc/HMMcopy.pdf
```


# For scDNA-seq data
I used functions in SCOPE to generate the readdepth matrix `Y_raw` for HMMcopy and then ran HMMcopy (see details in SCOPE section for the readdepth matrix part):

```r
library(HMMcopy)
library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)

base_path = "/data/maiziezhou_lab/Zihang/"
dataset = "/TNBC/DNA_bowtie2_rmdup/"
dataset_path = paste0(base_path, "CNV/dataset", dataset)
result_path = paste0(base_path, "CNV/HMMcopy/result", dataset)

data.type.pre = "DNA3020_rmdup"   #tumor
bamfolder.pre = paste0(dataset_path,"/",data.type.pre)
bamFile.pre <- list.files(bamfolder.pre, pattern = '*_rmdup.bam$')
bamdir.pre <- file.path(bamfolder.pre, bamFile.pre)
sampname_raw.pre <- sapply(strsplit(bamFile.pre, "_", fixed = TRUE), "[", 1)

data.type.mid = "DNA3022_rmdup"   #control
bamfolder.mid = paste0(dataset_path,"/", data.type.mid)
bamFile.mid <- list.files(bamfolder.mid, pattern = '*_rmdup.bam$')
bamdir.mid <- file.path(bamfolder.mid, bamFile.mid)
sampname_raw.mid <- sapply(strsplit(bamFile.mid, "_", fixed = TRUE), "[", 1)  

bamdir = c(bamdir.pre, bamdir.mid)
sampname_raw = c(sampname_raw.pre, sampname_raw.mid)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                         hgref = "hg38", resolution = 500)
ref_raw <- bambedObj$ref

# 2.2 Getting GC content and mappability
mapp <- get_mapp(ref_raw, hgref = "hg38")
# head(mapp)

gc <- get_gc(ref_raw, hgref = "hg38")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))

ref_raw = data.frame(ref_raw)
names(ref_raw)[1] = "chr"
names(ref_raw)[7] = "map"
ref_raw$chr = factor(ref_raw$chr, levels = unique(ref_raw$chr))
ref_raw = ref_raw[,c(-4,-5)]
ref_raw$gc = ref_raw$gc/100

coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = 'paired-end', hgref = "hg38")
Y_raw <- coverageObj$Y

segments_list = list()
for(ii in seq(ncol(Y_raw))){
    copy_temp = data.frame(ref_raw[,1:3], Y_raw[,ii], ref_raw[,4:5])
    names(copy_temp)[4] = "reads"
    
    copy_temp = data.table(copy_temp)
    tumour_copy <- correctReadcount(copy_temp)
    segments_list[[ii]] = HMMsegment(tumour_copy)
}

saveRDS(segments_list, paste0(result_path,"segments_list_500k.rds"))
```
