.libPaths(c("~/R/scDNA-seq", .libPaths()))
library(FLCNA)
library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg19)
library(WGSmapp)

# Edit paths for the dataset; use BAMs actually aligned to hg19.
samplename <- "TN1"
bamfolder <- "/path/to/hg19_cell_bams"
bam_pattern <- "\\.bam$"
outputfolder <- "/path/to/FLCNA/TN1"
dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)

bamFile <- list.files(bamfolder, pattern = bam_pattern)
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, "_", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw,
                         hgref = "hg19", resolution = 500)
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw, hgref = "hg19")
gc <- get_gc(ref_raw, hgref = "hg19")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40,
                                  seq = "single-end", hgref = "hg19")
Y_raw <- coverageObj$Y
QCobject <- FLCNA_QC(Y_raw = Y_raw, ref_raw = ref_raw,
                     mapp_thresh = 0.9, gc_thresh = c(20, 80))
log2Rdata <- FLCNA_normalization(Y = QCobject$Y, gc = QCobject$ref$gc,
                                map = QCobject$ref$mapp)
output_FLCNA <- FLCNA(K = c(3, 4, 5), lambda = 3, Y = data.matrix(log2Rdata))
# Original output call retained, including ref_raw.
CNA.output <- CNA.out(mean.matrix = output_FLCNA$mu.hat.best, LRR = log2Rdata,
                      Clusters = output_FLCNA$s.hat.best, ref = ref_raw,
                      cutoff = 0.35, L = 100)
saveRDS(output_FLCNA, file.path(outputfolder, "output_FLCNA"))
saveRDS(CNA.output, file.path(outputfolder, "CNA.output"))
