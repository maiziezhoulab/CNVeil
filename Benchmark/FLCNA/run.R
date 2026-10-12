#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run.R HG19_CELL_BAM_DIR NEW_OUTPUT_DIR")
input_dir <- normalizePath(args[1], mustWork = TRUE)
output_dir <- args[2]
if (file.exists(output_dir)) stop("Use a new output directory")
if (nzchar(Sys.getenv("FLCNA_R_LIB"))) .libPaths(c(Sys.getenv("FLCNA_R_LIB"), .libPaths()))
suppressPackageStartupMessages({
    library(FLCNA); library(SCOPE); library(BSgenome.Hsapiens.UCSC.hg19); library(WGSmapp)
})
# Select actual hg19 BAMs; filenames alone do not establish the alignment build.
bams <- sort(list.files(input_dir, pattern = "\\.bam$", full.names = TRUE))
if (!length(bams)) stop("No BAM files found")
cell_names <- sub("_.*$", "", basename(bams))
if (anyDuplicated(cell_names)) stop("Duplicate cell IDs before the first underscore")
# Retain the supplied package API and matrix orientation; do not silently migrate.
if (!"Y" %in% names(formals(FLCNA::FLCNA)) || !"LRR" %in% names(formals(FLCNA::CNA.out))) {
    stop("Installed FLCNA API differs from the supplied workflow; use the original working version")
}
dir.create(output_dir, recursive = TRUE)
setwd(normalizePath(output_dir))
bambedObj <- get_bam_bed(bamdir = bams, sampname = cell_names, hgref = "hg19", resolution = 500)
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw, hgref = "hg19")
gc <- get_gc(ref_raw, hgref = "hg19")
S4Vectors::mcols(ref_raw) <- cbind(S4Vectors::mcols(ref_raw), S4Vectors::DataFrame(gc, mapp))
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = "single-end", hgref = "hg19")
QCobject <- FLCNA_QC(Y_raw = coverageObj$Y, ref_raw = ref_raw,
                     mapp_thresh = 0.9, gc_thresh = c(20, 80))
log2Rdata <- FLCNA_normalization(Y = QCobject$Y, gc = QCobject$ref$gc, map = QCobject$ref$mapp)
output_FLCNA <- FLCNA(K = c(3, 4, 5), lambda = 3, Y = data.matrix(log2Rdata))
# Use retained-bin coordinates after QC, instead of the original unfiltered ref.
CNA.output <- CNA.out(mean.matrix = output_FLCNA$mu.hat.best, LRR = log2Rdata,
                      Clusters = output_FLCNA$s.hat.best, ref = QCobject$ref,
                      cutoff = 0.35, L = 100)
saveRDS(output_FLCNA, "output_FLCNA.rds")
saveRDS(CNA.output, "CNA.output.rds")
saveRDS(QCobject$ref, "retained_bins.rds")
writeLines(bams, "input_files.txt")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
