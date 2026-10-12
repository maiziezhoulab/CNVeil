#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript run.R HG19_CELL_BAM_DIR NEW_OUTPUT_DIR")
input_dir <- normalizePath(args[1], mustWork = TRUE)
output_dir <- args[2]
if (file.exists(output_dir)) stop("Use a new output directory")
if (nzchar(Sys.getenv("SCOPE_R_LIB"))) .libPaths(c(Sys.getenv("SCOPE_R_LIB"), .libPaths()))
suppressPackageStartupMessages({
    library(SCOPE); library(BSgenome.Hsapiens.UCSC.hg19); library(WGSmapp)
})
bams <- sort(list.files(input_dir, pattern = "\\.bam$", full.names = TRUE))
if (!length(bams)) stop("No BAM files found")
cell_names <- sub("_.*$", "", basename(bams))
if (anyDuplicated(cell_names)) stop("Duplicate cell IDs before the first underscore")
dir.create(output_dir, recursive = TRUE)
setwd(normalizePath(output_dir))

# Rebuild hg19 preprocessing; do not reuse the original hg38 qcObj checkpoint.
bambedObj <- get_bam_bed(bamdir = bams, sampname = cell_names, hgref = "hg19", resolution = 500)
ref_raw <- bambedObj$ref
mapp <- get_mapp(ref_raw, hgref = "hg19")
gc <- get_gc(ref_raw, hgref = "hg19")
S4Vectors::mcols(ref_raw) <- cbind(S4Vectors::mcols(ref_raw), S4Vectors::DataFrame(gc, mapp))
coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40, seq = "paired-end", hgref = "hg19")
QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = coverageObj$Y, sampname_raw = cell_names,
                     ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)
saveRDS(qcObj, "qcObj.rds")
Y <- qcObj$Y
ref <- qcObj$ref
if (is.null(colnames(Y))) colnames(Y) <- qcObj$sampname
if (!nrow(Y) || !ncol(Y) || nrow(Y) != length(ref)) stop("Empty or misaligned QC data")
Gini <- get_gini(Y)
norm_index <- which(Gini <= 0.12)
if (!length(norm_index)) stop("No normal-reference cells selected at Gini <= 0.12")
normObj <- normalize_codex2_ns_noK(Y_qc = Y, gc_qc = ref$gc, norm_index = norm_index)
ploidy.sim <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)
normObj.scope <- normalize_scope(Y_qc = Y, gc_qc = ref$gc, K = 1,
                                 ploidyInt = ploidy.sim, norm_index = norm_index,
                                 T = 1:5, beta0 = normObj$beta.hat)
Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
chrs <- unique(as.character(GenomicRanges::seqnames(ref)))
segment_cs <- setNames(vector("list", length(chrs)), chrs)
for (chr in chrs) {
    message("Segmenting ", chr)
    segment_cs[[chr]] <- segment_CBScs(Y = Y, Yhat = Yhat, sampname = colnames(Y),
                                      ref = ref, chr = chr, mode = "integer", max.ns = 1)
}
iCN <- do.call(rbind, lapply(segment_cs, function(z) z[["iCN"]]))
saveRDS(segment_cs, "segment_cs.rds")
saveRDS(iCN, "iCN.rds")
writeLines(bams, "input_files.txt")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
