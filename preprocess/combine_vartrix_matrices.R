#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(optparse)
})

option_list <- list(
  make_option(c("--dir"), type = "character", default = "./",
              help = "Directory containing VarTrix output matrices"),
  make_option(c("--vcf-dir"), type = "character", default = "./vcf_sub/",
              help = "Directory containing per-chromosome VCF files"),
  make_option(c("--barcodes"), type = "character", default = "./barcodes.tsv",
              help = "Barcode file, one barcode per line"),
  make_option(c("--sample"), type = "character", default = "Sample",
              help = "Sample name used in VCF filenames"),
  make_option(c("--out-dir"), type = "character", default = "./",
              help = "Directory to save combined RDS files"),
  make_option(c("--plot-qc"), action = "store_true", default = FALSE,
              help = "Generate QC plots"),
  make_option(c("--chromosomes"), type = "character", default = "1:22",
              help = "Chromosomes to process, e.g. '1:22' or '1,2,3'")
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_chromosomes <- function(chr_string) {
  if (grepl(":", chr_string)) {
    parts <- as.integer(strsplit(chr_string, ":")[[1]])
    return(seq(parts[1], parts[2]))
  } else {
    return(as.integer(strsplit(chr_string, ",")[[1]]))
  }
}

combine_vartrix_matrices <- function(dir_path,
                                     vcf_dir,
                                     barcodes_path,
                                     sample_name,
                                     out_dir,
                                     plot_qc = FALSE,
                                     chromosomes = 1:22) {
  if (!dir.exists(dir_path)) {
    stop("Input matrix directory does not exist: ", dir_path)
  }

  if (!dir.exists(vcf_dir)) {
    stop("VCF directory does not exist: ", vcf_dir)
  }

  if (!file.exists(barcodes_path)) {
    stop("Barcode file does not exist: ", barcodes_path)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  barcodes <- read.table(
    barcodes_path,
    sep = "\t",
    stringsAsFactors = FALSE
  )[, 1]

  alt_all <- NULL
  ref_all <- NULL
  var_all <- NULL

  message("Combining VarTrix matrices...")

  for (chr in chromosomes) {
    alt_file <- file.path(dir_path, paste0("chr", chr, "_matrix_alt"))
    ref_file <- file.path(dir_path, paste0("chr", chr, "_matrix_ref"))

    vcf_file <- file.path(
      vcf_dir,
      paste0(sample_name, "_phased_chr", chr, "_chrprefix.vcf")
    )

    if (!file.exists(alt_file)) {
      warning("ALT matrix missing for chr", chr, "; skipping.")
      next
    }

    if (!file.exists(ref_file)) {
      warning("REF matrix missing for chr", chr, "; skipping.")
      next
    }

    if (!file.exists(vcf_file)) {
      warning("VCF file missing for chr", chr, "; skipping.")
      next
    }

    message("Processing chr", chr)

    alt <- readMM(alt_file)
    ref <- readMM(ref_file)
    var <- read.table(vcf_file, header = FALSE, stringsAsFactors = FALSE)

    if (nrow(alt) != nrow(var)) {
      stop(
        "Row mismatch for chr", chr,
        ": ALT matrix has ", nrow(alt),
        " rows but VCF has ", nrow(var), " rows."
      )
    }

    if (ncol(alt) != length(barcodes)) {
      stop(
        "Column mismatch for chr", chr,
        ": matrix has ", ncol(alt),
        " columns but barcode file has ", length(barcodes), " barcodes."
      )
    }

    keep <- rowSums(alt) > 0

    alt_all <- rbind(alt_all, alt[keep, ])
    ref_all <- rbind(ref_all, ref[keep, ])
    var_all <- rbind(var_all, var[keep, ])
  }

  if (is.null(alt_all) || nrow(alt_all) == 0) {
    stop("No valid SNPs were retained. Please check input matrices and VCF files.")
  }

  colnames(alt_all) <- barcodes
  colnames(ref_all) <- barcodes

  message("Total SNP number: ", nrow(alt_all))
  message("Total cell number: ", ncol(alt_all))

  saveRDS(alt_all, file.path(out_dir, "alt_all.rds"))
  saveRDS(ref_all, file.path(out_dir, "ref_all.rds"))
  saveRDS(var_all, file.path(out_dir, "var_all.rds"))

  message("Saved:")
  message("  ", file.path(out_dir, "alt_all.rds"))
  message("  ", file.path(out_dir, "ref_all.rds"))
  message("  ", file.path(out_dir, "var_all.rds"))

  total_all <- alt_all + ref_all
  af <- rowSums(alt_all) / rowSums(total_all)
  af[is.na(af)] <- 0

  if (plot_qc) {
    qc_pdf <- file.path(out_dir, paste0(sample_name, "_vartrix_qc.pdf"))

    pdf(qc_pdf, width = 8, height = 8)
    par(mfrow = c(3, 1))

    hist(
      colSums(total_all),
      main = paste0(sample_name, " cell-level coverage"),
      xlab = "Total coverage per cell",
      ylab = "Frequency",
      breaks = 100
    )

    hist(
      rowSums(total_all),
      main = paste0(sample_name, " SNP-level coverage"),
      xlab = "Total coverage per SNP",
      ylab = "Frequency",
      breaks = 100,
      xlim = c(0, min(100, max(rowSums(total_all))))
    )

    hist(
      af,
      breaks = 100,
      main = paste0(sample_name, " pooled VAF per SNP"),
      xlab = "Pooled VAF",
      ylab = "Frequency"
    )

    dev.off()

    message("QC plot saved to: ", qc_pdf)
  }

  invisible(list(
    alt = alt_all,
    ref = ref_all,
    var = var_all
  ))
}

chromosomes <- parse_chromosomes(opt$chromosomes)

combine_vartrix_matrices(
  dir_path = opt$dir,
  vcf_dir = opt$`vcf-dir`,
  barcodes_path = opt$barcodes,
  sample_name = opt$sample,
  out_dir = opt$`out-dir`,
  plot_qc = opt$`plot-qc`,
  chromosomes = chromosomes
)
