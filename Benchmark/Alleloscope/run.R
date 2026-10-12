#!/usr/bin/env Rscript
# TN1 workflow. Paths are arguments; benchmark choices are explicit below.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
    stop("Usage: Rscript run.R SAMPLE VARTRIX_DIR COVERAGE_CSV NORMAL_RDS SIZES_TSV OUTPUT_DIR")
}
sample_name <- args[1]
input_dir <- args[2]
coverage_csv <- args[3]
normal_rds <- args[4]
sizes_tsv <- args[5]
output_dir <- args[6]
if (nzchar(Sys.getenv("ALLELOSCOPE_R_LIB"))) {
    .libPaths(c(Sys.getenv("ALLELOSCOPE_R_LIB"), .libPaths()))
}
suppressPackageStartupMessages({library(Alleloscope); library(Matrix)})

# Preserve TN1 settings. Do not automatically apply these to another dataset.
reference_region <- Sys.getenv("REFERENCE_REGION", "8:5e+05")
coverage_multiplier <- as.numeric(Sys.getenv("COVERAGE_MULTIPLIER", "2"))
coverage_columns <- seq(1L, 124L, by = 4L)
if (!is.finite(coverage_multiplier) || coverage_multiplier <= 0) stop("Invalid multiplier")
if (sample_name != "TN1") warning("Review the TN1 reference and coverage-column settings for this sample")

require_file <- function(path) {
    if (!file.exists(path)) stop("Missing file: ", path)
}
check_ids <- function(ids, label) {
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
        stop(label, " must contain unique, nonempty identifiers")
    }
}
# Strip known SRR filename suffixes for matching; retain original barcode labels.
cell_key <- function(ids) sub("^(SRR[0-9]+).*$", "\\1", ids)
bin_key <- function(ids) sub("^(chr[^.]+)\\.([0-9]+)\\.([0-9]+)$", "\\1-\\2-\\3", ids)

Cbn_matrix <- function(dir_path, samplename) {
    barcode_path <- file.path(dir_path, "barcodes.tsv")
    require_file(barcode_path)
    barcodes <- readLines(barcode_path)
    check_ids(barcodes, "Vartrix barcodes")
    alt_parts <- ref_parts <- var_parts <- vector("list", 22L)
    for (chr in seq_len(22L)) {
        alt_path <- file.path(dir_path, "chr", paste0("chr", chr, "_matrix_alt"))
        ref_path <- file.path(dir_path, "chr", paste0("chr", chr, "_matrix_ref"))
        vcf_path <- file.path(dir_path, "vcf", paste0(samplename, "_phased_chr", chr, "_chrprefix.vcf"))
        for (path in c(alt_path, ref_path, vcf_path)) require_file(path)
        alt <- Matrix::readMM(alt_path)
        ref <- Matrix::readMM(ref_path)
        variants <- read.table(vcf_path, header = FALSE, comment.char = "#",
                               sep = "\t", quote = "", stringsAsFactors = FALSE)
        if (!identical(dim(alt), dim(ref)) || nrow(alt) != nrow(variants) ||
            ncol(alt) != length(barcodes)) stop("Matrix/VCF/barcode dimensions disagree on chromosome ", chr)
        # Historical files named 'alt' contained reference counts, so the
        # original positive-'alt' filter actually selected positive REF coverage.
        keep <- Matrix::rowSums(ref) > 0
        alt_parts[[chr]] <- alt[keep, , drop = FALSE]
        ref_parts[[chr]] <- ref[keep, , drop = FALSE]
        var_parts[[chr]] <- variants[keep, , drop = FALSE]
        message("Combined chromosome ", chr)
    }
    alt_all <- do.call(rbind, alt_parts)
    ref_all <- do.call(rbind, ref_parts)
    var_all <- do.call(rbind, var_parts)
    if (!nrow(alt_all)) stop("No SNPs remain after filtering")
    colnames(alt_all) <- colnames(ref_all) <- barcodes
    list(alt_all = alt_all, ref_all = ref_all, var_all = var_all)
}

for (path in c(coverage_csv, normal_rds, sizes_tsv)) require_file(path)
if (dir.exists(output_dir)) stop("Use a new output directory to avoid overwriting previous results")
dir.create(output_dir, recursive = TRUE)
output_dir <- normalizePath(output_dir)
alleles <- Cbn_matrix(input_dir, sample_name)
for (name in names(alleles)) saveRDS(alleles[[name]], file.path(output_dir, paste0(name, ".rds")))

# The first CSV column contains cells; the remaining columns contain bins.
coverage <- read.csv(coverage_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(coverage) < 2L || nrow(coverage) < 1L) stop("Empty coverage CSV")
if (!all(vapply(coverage[-1L], is.numeric, logical(1)))) stop("Coverage bins must be numeric")
raw_counts <- t(as.matrix(coverage[-1L]))
rownames(raw_counts) <- bin_key(names(coverage)[-1L])
colnames(raw_counts) <- cell_key(as.character(coverage[[1L]]))
check_ids(rownames(raw_counts), "Coverage bins")
check_ids(colnames(raw_counts), "Coverage cells")
allele_keys <- cell_key(colnames(alleles$alt_all))
check_ids(allele_keys, "Normalized allele-matrix cell IDs")
cell_order <- match(colnames(raw_counts), allele_keys)
if (anyNA(cell_order)) stop("Coverage cells missing from allele matrices: ",
                           paste(colnames(raw_counts)[is.na(cell_order)], collapse = ", "))
alleles$alt_all <- alleles$alt_all[, cell_order, drop = FALSE]
alleles$ref_all <- alleles$ref_all[, cell_order, drop = FALSE]
colnames(alleles$alt_all) <- colnames(alleles$ref_all) <- colnames(raw_counts)
barcodes <- data.frame(V1 = colnames(raw_counts), stringsAsFactors = FALSE)

ref_counts <- as.matrix(readRDS(normal_rds))
if (is.null(rownames(ref_counts))) {
    if (Sys.getenv("TRUST_NORMAL_BIN_ORDER") != "1") {
        stop("Normal reference has no bin IDs. Verify its bin order, then set TRUST_NORMAL_BIN_ORDER=1 if appropriate")
    }
    if (nrow(ref_counts) != nrow(raw_counts)) stop("Normal bin count differs")
    rownames(ref_counts) <- rownames(raw_counts)
} else {
    rownames(ref_counts) <- bin_key(rownames(ref_counts))
    check_ids(rownames(ref_counts), "Normal-reference bins")
    normal_order <- match(rownames(raw_counts), rownames(ref_counts))
    if (anyNA(normal_order)) stop("Normal reference is missing coverage bins")
    ref_counts <- ref_counts[normal_order, , drop = FALSE]
}
for (m in list(raw_counts, ref_counts)) {
    if (any(!is.finite(m)) || any(m < 0)) stop("Coverage must be finite and nonnegative")
}
size <- read.table(sizes_tsv, stringsAsFactors = FALSE)

# All Alleloscope-generated files are written inside this run's directory.
setwd(output_dir)
checkpoint <- Sys.getenv("ALLELOSCOPE_CHECKPOINT")
if (nzchar(checkpoint)) {
    require_file(checkpoint)
    warning("Resuming an existing object: verify it was created using these same cells, bins and allele labels")
    Obj_filtered <- readRDS(checkpoint)
} else {
    Obj <- Createobj(alt_all = alleles$alt_all, ref_all = alleles$ref_all,
                     var_all = alleles$var_all, samplename = sample_name,
                     genome_assembly = "GRCh37", dir_path = "./", barcodes = barcodes,
                     size = size, assay = "scDNAseq")
    Obj_filtered <- Matrix_filter(Obj = Obj, cell_filter = 1000, SNP_filter = 20,
                                  min_vaf = 0.1, max_vaf = 0.9)
    Obj_filtered <- Segmentation(Obj_filtered = Obj_filtered, raw_counts = raw_counts,
                                 ref_counts = ref_counts, plot_seg = TRUE)
    Obj_filtered <- Segments_filter(Obj_filtered = Obj_filtered, nSNP = 2000)
    Obj_filtered <- Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000,
                                plot_stat = TRUE, cont = TRUE)
    saveRDS(Obj_filtered, "Obj_filtered_pre_normal.rds")
}
Obj_filtered <- Select_normal(Obj_filtered = Obj_filtered, raw_counts = raw_counts,
                              plot_theta = TRUE)
# TN1 override: automatic selection yielded predominantly diploid calls.
Obj_filtered$ref <- reference_region
Obj_filtered <- Genotype_value(Obj_filtered = Obj_filtered, type = "cellline",
                               raw_counts = raw_counts, ref_counts = ref_counts)
saveRDS(Obj_filtered, "Obj_filtered_before_scaling.rds")
if (coverage_multiplier != 1) {
    if (ncol(Obj_filtered$genotype_values) < max(coverage_columns)) {
        stop("Original TN1 scaling requires at least 121 genotype-value columns; review coverage_columns")
    }
    Obj_filtered$genotype_values[, coverage_columns] <-
        Obj_filtered$genotype_values[, coverage_columns] * coverage_multiplier
}
Obj_filtered <- Genotype(Obj_filtered = Obj_filtered)
saveRDS(Obj_filtered, "Obj_filtered.rds")
# Retained tool-native lineage output from the supplied workflow.
linplot <- Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000, nclust = 10)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
message("Saved results to ", output_dir)
