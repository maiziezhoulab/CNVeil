#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript run.R HG19_CELL_BAM_DIR CORRECTION_R NEW_OUTPUT_DIR")
input_dir <- normalizePath(args[1], mustWork = TRUE)
correction_file <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
if (file.exists(output_dir)) stop("Use a new output directory")
if (nzchar(Sys.getenv("HMMCOPY_R_LIB"))) .libPaths(c(Sys.getenv("HMMCOPY_R_LIB"), .libPaths()))
suppressPackageStartupMessages({
    library(HMMcopy); library(SCOPE); library(BSgenome.Hsapiens.UCSC.hg19)
    library(WGSmapp); library(data.table)
})
source(correction_file)
if (!exists("correctReadcount_github", mode = "function")) stop("Custom correction.R must define correctReadcount_github()")
bams <- sort(list.files(input_dir, pattern = "\\.bam$", full.names = TRUE))
if (!length(bams)) stop("No BAM files found")
cell_names <- sub("_.*$", "", basename(bams))
if (anyDuplicated(cell_names)) stop("Duplicate cell IDs before the first underscore")
dir.create(output_dir, recursive = TRUE)
setwd(normalizePath(output_dir))
bambedObj <- get_bam_bed(bamdir = bams, sampname = cell_names, hgref = "hg19", resolution = 500)
bins <- bambedObj$ref
gc <- get_gc(bins, hgref = "hg19") / 100
mapp <- get_mapp(bins, hgref = "hg19")
# Construct named columns; do not depend on GRanges data-frame column positions.
ref_raw <- data.frame(chr = factor(as.character(GenomicRanges::seqnames(bins)),
                                  levels = unique(as.character(GenomicRanges::seqnames(bins)))),
                      start = BiocGenerics::start(bins), end = BiocGenerics::end(bins),
                      gc = gc, map = mapp)
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = "single-end", hgref = "hg19")
Y_raw <- coverageObj$Y
if (nrow(Y_raw) != nrow(ref_raw)) stop("Coverage and reference bin counts differ")
if (any(!is.finite(Y_raw))) stop("Nonfinite coverage")
if (is.null(colnames(Y_raw))) {
    if (ncol(Y_raw) != length(cell_names)) stop("Cell counts differ")
    colnames(Y_raw) <- cell_names
}
keep <- colSums(Y_raw) > 0
if (!any(keep)) stop("All cells have zero coverage")
Y <- Y_raw[, keep, drop = FALSE]
cell_names <- colnames(Y)
segments_list <- setNames(vector("list", ncol(Y)), cell_names)
saveRDS(Y_raw, "Y_raw.rds")
saveRDS(Y, "Y.rds")
saveRDS(ref_raw, "ref_raw.rds")
saveRDS(cell_names, "cell_names.rds")
for (i in seq_len(ncol(Y))) {
    copy_temp <- data.table(chr = ref_raw$chr, start = ref_raw$start, end = ref_raw$end,
                            reads = Y[, i], gc = ref_raw$gc, map = ref_raw$map)
    message("Processing ", cell_names[i])
    tumour_copy <- correctReadcount_github(copy_temp)
    segments_list[[i]] <- HMMsegment(tumour_copy)
}
saveRDS(segments_list, "segments_list_500k.rds")
writeLines(bams, "input_files.txt")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
