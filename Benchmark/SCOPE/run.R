.libPaths(c("~/R/scDNA-seq", .libPaths()))
library(SCOPE)
library(BSgenome.Hsapiens.UCSC.hg19)
library(WGSmapp)

# Original checkpoint-based workflow. This checkpoint must be hg19.
qc_file <- "/path/to/hg19/qcObj.rds"
outputfolder <- "/path/to/SCOPE/output"
dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)

# Optional original preprocessing, translated to hg19; leave commented when resuming.
# bamfolder <- "/path/to/hg19_cell_bams"
# bamFile <- list.files(bamfolder, pattern = "\\.bam$")
# bamdir <- file.path(bamfolder, bamFile)
# sampname_raw <- sapply(strsplit(bamFile, "_", fixed = TRUE), "[", 1)
# bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw,
#                          hgref = "hg19", resolution = 500)
# ref_raw <- bambedObj$ref
# mapp <- get_mapp(ref_raw, hgref = "hg19")
# gc <- get_gc(ref_raw, hgref = "hg19")
# values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
# coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40,
#                                   seq = "paired-end", hgref = "hg19")
# QCmetric_raw <- get_samp_QC(bambedObj)
# qcObj <- perform_qc(Y_raw = coverageObj$Y, sampname_raw = sampname_raw,
#                      ref_raw = ref_raw, QCmetric_raw = QCmetric_raw)
# saveRDS(qcObj, qc_file)

qcObj <- readRDS(qc_file)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric
Gini <- get_gini(Y)
norm_index <- which(Gini <= 0.12)
normObj <- normalize_codex2_ns_noK(Y_qc = Y, gc_qc = ref$gc,
                                  norm_index = norm_index)
ploidy.sim <- initialize_ploidy(Y = Y, Yhat = normObj$Yhat, ref = ref)
normObj.scope <- normalize_scope(Y_qc = Y, gc_qc = ref$gc, K = 1,
                                 ploidyInt = ploidy.sim, norm_index = norm_index,
                                 T = 1:5, beta0 = normObj$beta.hat)
Yhat <- normObj.scope$Yhat[[which.max(normObj.scope$BIC)]]
fGC.hat <- normObj.scope$fGC.hat[[which.max(normObj.scope$BIC)]]
chrs <- unique(as.character(seqnames(ref)))
segment_cs <- vector("list", length = length(chrs))
names(segment_cs) <- chrs
for (chri in chrs) {
    message("\n", chri, "\n")
    segment_cs[[chri]] <- segment_CBScs(Y = Y, Yhat = Yhat,
                                        sampname = colnames(Y), ref = ref,
                                        chr = chri, mode = "integer", max.ns = 1)
}
iCN <- do.call(rbind, lapply(segment_cs, function(z) z[["iCN"]]))
saveRDS(segment_cs, file.path(outputfolder, "segment_cs.rds"))
