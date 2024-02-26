# Orginal Vignettes
The package paths are:

```bash
https://www.bioconductor.org/packages/release/bioc/html/SCOPE.html

https://github.com/rujinwang/SCOPE
```

The path of vignettes written by the authors is:

```bash
http://bioconductor.org/packages/devel/bioc/vignettes/SCOPE/inst/doc/SCOPE_vignette.html
```
# Installation

You should load R/4.0.5 from accre, then install SCOPE according to the instruction in their github. The instructions are also copied here.

```r
install.packages('devtools')
devtools::install_github("rujinwang/WGSmapp")
devtools::install_github("rujinwang/SCOPE")
```

# For scDNA-seq data
### 1. Loading
Loading the bam files:

```r
library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(WGSmapp)

align.type = "bowtie2"
database_path = "/data/maiziezhou_lab/Zihang/CNV/dataset/TNBC/DNA_"

data.type.pre = "3020"   #tumor
bamfolder.pre = paste0(database_path,align.type,"_rmdup/DNA",data.type.pre,"_rmdup")
bamFile.pre <- list.files(bamfolder.pre, pattern = '*_rmdup.bam$')
bamdir.pre <- file.path(bamfolder.pre, bamFile.pre)
sampname_raw.pre <- sapply(strsplit(bamFile.pre, "_", fixed = TRUE), "[", 1)

data.type.mid = "3022"   #control
bamfolder.mid = paste0(database_path,align.type,"_rmdup/DNA",data.type.mid,"_rmdup")
bamFile.mid <- list.files(bamfolder.mid, pattern = '*_rmdup.bam$')
bamdir.mid <- file.path(bamfolder.mid, bamFile.mid)
sampname_raw.mid <- sapply(strsplit(bamFile.mid, "_", fixed = TRUE), "[", 1)

bamdir = c(bamdir.pre, bamdir.mid)
sampname_raw = c(sampname_raw.pre, sampname_raw.mid)
```

### 2. Pre-computation and Quality Control
#### 2.1 Pre-preparation

```r
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                         hgref = "hg38", resolution = 5000)
ref_raw <- bambedObj$ref
```

#### 2.2 Getting GC content and mappability

```r
mapp <- get_mapp(ref_raw, hgref = "hg38")

gc <- get_gc(ref_raw, hgref = "hg38")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
```

#### 2.3 Getting coverage

```r
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, 
                                  seq = 'paired-end', hgref = "hg38")
```
#### 2.4 Quality control

```r
QCmetric_raw <- get_samp_QC(bambedObj)
Y_raw <- coverageObj$Y

qcObj <- perform_qc(Y_raw = Y_raw, 
                    sampname_raw = sampname_raw, ref_raw = ref_raw, 
                    QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric
```

Note that:
- `hgref` is the genome reference.


### 3. Running SCOPE
If normal cells are not applicable, you can apply gini coefficient to generate the normal cells index:

```r
Gini <- get_gini(Y)
norm_index = which(Gini<=0.12)
```
otherwise, if normal cells are not applicable:

```r
norm_index = which(sampname %in% sampname_raw.mid)
```

SCOPE first used CODEX2 with no latent factors to gain initial parameters and then run SCOPE to estimate the Poisson Distribution parameters:

```r
# first-pass CODEX2 run with no latent factors
normObj <- normalize_codex2_ns_noK(Y_qc = Y,
                                   gc_qc = ref$gc,
                                   norm_index = norm_index)
# Ploidy initialization
ploidy <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)

normObj.scope <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref$gc,
                                         K = 1, ploidyInt = ploidy,
                                         norm_index = norm_index, T = 1:5,
                                         beta0 = normObj$beta.hat, nCores = 2)
```

If group information is available, you can choose to use the group-wise version:

```r
# first-pass CODEX2 run with no latent factors
normObj <- normalize_codex2_ns_noK(Y_qc = Y,
                                   gc_qc = ref$gc,
                                   norm_index = norm_index)
# Group-wise ploidy initialization
clones <- rep(0, length(sampname))
clones[which(sampname %in% sampname_raw.pre)] = "pre"
clones[which(sampname %in% sampname_raw.mid)] = "mid"

ploidy <- initialize_ploidy_group(Y = Y, Yhat = normObj$Yhat,
                                  ref = ref, groups = clones)

# Group-wise normalization
normObj.scope <- normalize_scope_group(Y_qc = Y,
                                       gc_qc = ref$gc,
                                       K = 1, ploidyInt = ploidy,
                                       norm_index = which(clones %in% "mid"),
                                       groups = clones,
                                       T = 1:5,
                                       beta0 = normObj$beta.hat)
```


You can then perform segmentation and infer CNVs using SCOPE:

```r
Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]

# Cross-sample segmentation by SCOPE
chrs <- unique(as.character(seqnames(ref)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs
for (chri in chrs) {
  message('\n', chri, '\n')
  segment_cs[[chri]] <- segment_CBScs(Y = Y,
                                      Yhat = Yhat,
                                      sampname = colnames(Y),
                                      ref = ref,
                                      chr = chri,
                                      mode = "integer", max.ns = 1)
}
iCN <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))

saveRDS(segment_cs, "segment_cs.rds")
saveRDS(qcObj, "qcObj.rds")
```
