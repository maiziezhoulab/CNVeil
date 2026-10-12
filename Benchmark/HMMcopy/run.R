.libPaths(c("~/R/HMMcopy", .libPaths()))
library(HMMcopy)
library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg19)
library(WGSmapp)
library(data.table)

bamfolder <- "/path/to/hg19_cell_bams"
bam_pattern <- "_grch19_rmdup_rg_cb\\.bam$"
result_path <- "/path/to/HMMcopy/TN1"
correction_script <- "/path/to/original/correction.R"
source(correction_script)
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

bamFile <- list.files(bamfolder, pattern = bam_pattern)
bamdir <- file.path(bamfolder, bamFile)
sampname_raw <- sapply(strsplit(bamFile, "_", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw,
                         hgref = "hg19", resolution = 500)
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw, hgref = "hg19")
gc <- get_gc(ref_raw, hgref = "hg19")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
ref_raw <- data.frame(ref_raw)
names(ref_raw)[1] <- "chr"
names(ref_raw)[7] <- "map"
ref_raw$chr <- factor(ref_raw$chr, levels = unique(ref_raw$chr))
ref_raw <- ref_raw[, c(-4, -5)]
ref_raw$gc <- ref_raw$gc / 100
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40,
                                  seq = "single-end", hgref = "hg19")
Y_raw <- coverageObj$Y
QCmetric_raw <- get_samp_QC(bambedObj)
# Same zero-coverage exclusion, retaining matrix dimensions for a single cell.
keep <- colSums(Y_raw) > 0
Y <- Y_raw[, keep, drop = FALSE]
cell_names <- matrix(colnames(Y_raw)[keep], nrow = 1)
saveRDS(Y_raw, file.path(result_path, "Y_raw_se"))
saveRDS(cell_names, file.path(result_path, "cell_names"))
saveRDS(Y, file.path(result_path, "Y"))
saveRDS(Y_raw, file.path(result_path, "Y_raw"))
saveRDS(ref_raw, file.path(result_path, "ref_raw"))
segments_list <- list()
for (ii in seq_len(ncol(Y))) {
    copy_temp <- data.frame(ref_raw[, 1:3], Y[, ii], ref_raw[, 4:5])
    names(copy_temp)[4] <- "reads"
    copy_temp <- data.table(copy_temp)
    print(ii)
    tumour_copy <- correctReadcount_github(copy_temp)
    segments_list[[ii]] <- HMMsegment(tumour_copy)
}
saveRDS(segments_list, file.path(result_path, "segments_list_500k.rds"))
