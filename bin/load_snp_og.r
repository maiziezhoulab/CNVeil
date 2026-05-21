# .libPaths("/home/luoc1/R/rlib-4.0.5")
## -----------------------------
## Set R library path from conda
## -----------------------------
conda_path <- Sys.getenv("CONDA_PREFIX")

if (conda_path == "") {
  conda_path <- dirname(dirname(Sys.which("R")))
}

lib_path <- file.path(conda_path, "lib", "R", "library")
.libPaths(c(lib_path, .libPaths()))

message("R library paths:")
print(.libPaths())


merge_snp_cell_counts <- function(var_all, ref_all, alt_all, out_tsv,
                                  min_total = 1L, integerize = TRUE) {
  # Needs: Matrix, data.table
  stopifnot(requireNamespace("Matrix", quietly = TRUE),
            requireNamespace("data.table", quietly = TRUE))

  # sanity checks
  if (nrow(var_all) != nrow(ref_all) ||
      nrow(var_all) != nrow(alt_all)) {
    stop("Row mismatch: nrow(var_all) must equal nrow(ref_all) and nrow(alt_all).")
  }
  if (ncol(ref_all) != ncol(alt_all)) {
    stop("Column mismatch: ncol(ref_all) must equal ncol(alt_all).")
  }

  # coerce to triplet form (dgTMatrix) and summarize nonzeros
  refT <- methods::as(ref_all, "dgTMatrix")
  altT <- methods::as(alt_all, "dgTMatrix")

  sr <- Matrix::summary(refT)  # columns: i, j, x (1-based)
  sa <- Matrix::summary(altT)

  # data.table for speed; merge on (i,j)
  dt_ref <- data.table::as.data.table(sr)
  data.table::setnames(dt_ref, c("i","j","x"), c("i","j","ref"))
  dt_alt <- data.table::as.data.table(sa)
  data.table::setnames(dt_alt, c("i","j","x"), c("i","j","alt"))

  dt <- dt_ref[dt_alt, on = c("i","j")]
  # fill missing with 0
  if (!"ref" %in% names(dt)) dt[, ref := NA_real_]
  if (!"alt" %in% names(dt)) dt[, alt := NA_real_]
  dt[is.na(ref), ref := 0]
  dt[is.na(alt), alt := 0]

  # optionally integerize counts (dgTMatrix@x can be numeric)
  if (integerize) {
    dt[, `:=`(ref = as.integer(round(ref)),
              alt = as.integer(round(alt)))]
  }

  # filter by total depth (keep where ref+alt >= min_total)
  dt[, total := ref + alt]
  if (min_total > 0L) dt <- dt[total >= as.integer(min_total)]
  dt[, total := NULL]

  # map matrix indices to metadata:
  # i = row (variant), j = column (cell). Matrix indices are 1-based.
  chrom_vec <- var_all[[1]]         # V1
  pos_vec   <- var_all[[2]]         # V2
  ref_alle  <- var_all[[4]]         # V4 (REF)
  alt_alle  <- var_all[[5]]         # V5 (ALT)
  cell_ids  <- colnames(ref_all)

  # attach columns
  dt[, chrom := chrom_vec[i]]
  dt[, pos   := pos_vec[i]]
  dt[, cell  := if (is.null(cell_ids)) as.character(j) else cell_ids[j]]
  # (optional: include REF/ALT alleles from var_all if you want them in the TSV)
  # dt[, ref_allele := ref_alle[i]]
  # dt[, alt_allele := alt_alle[i]]

  # final ordering and write
  out <- dt[, .(chrom, pos, cell, ref, alt)]
  data.table::fwrite(out, file = out_tsv, sep = "\t", quote = FALSE)

  invisible(out)
}



run_merge_snp_cell_counts <- function(vcf_path, out_dir, sample,
                                      min_total = 1L,
                                      integerize = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  message(vcf_path)

  var_file <- file.path(vcf_path, "var_all.rds")
  alt_file <- file.path(vcf_path, "alt_all.rds")
  ref_file <- file.path(vcf_path, "ref_all.rds")

  var_all <- readRDS(var_file)
  alt_all <- readRDS(alt_file)
  ref_all <- readRDS(ref_file)

  out_file <- file.path(out_dir, paste0(sample, "_baf.tsv"))

  res <- merge_snp_cell_counts(
    var_all,
    ref_all,
    alt_all,
    out_tsv = out_file,
    min_total = min_total,
    integerize = integerize
  )

  return(res)
}

split_baf_by_chrom <- function(outfile, outdir, sample, num_threads = 4) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  cmds <- vapply(1:22, function(chr) {
    out_chr <- file.path(outdir, paste0(sample, "_baf_chr", chr, ".tsv"))

    sprintf(
      "grep -w %s %s > %s",
      shQuote(paste0("chrom\\|chr", chr)),
      shQuote(outfile),
      shQuote(out_chr)
    )
  }, character(1))

  cmd <- sprintf(
    "parallel -j %d ::: %s",
    as.integer(num_threads),
    paste(shQuote(cmds), collapse = " ")
  )

  message(cmd)
  system(cmd)
}

load_snp <- function(vcf_path, out_dir, sample, num_threads = 4,
                 min_total = 1L, integerize = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  message("Running sample: ", sample)
  message("Input RDS folder: ", vcf_path)
  message("Output folder: ", out_dir)

  # 1. Merge ref/alt SNP counts across cells
  res <- run_merge_snp_cell_counts(
    vcf_path = vcf_path,
    out_dir = out_dir,
    sample = sample,
    min_total = min_total,
    integerize = integerize
  )

  # 2. Split merged BAF file by chromosome
  baf_file <- file.path(out_dir, paste0(sample, "_baf.tsv"))

  split_baf_by_chrom(
    outfile = baf_file,
    outdir = out_dir,
    sample = sample,
    num_threads = num_threads
  )

  message("Done.")
  message("Merged BAF file: ", baf_file)
  message("Per-chromosome BAF files written to: ", out_dir)

  invisible(res)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript load_snp.R <vcf_path> <work_dir> <sample> [num_threads]")
}

vcf_path <- args[1]
work_dir <- args[2]
out_dir <- file.path(work_dir, "Allele_CN", "allele_count_by_cell")
sample <- args[3]
num_threads <- ifelse(length(args) >= 4, as.integer(args[4]), 4)

load_snp(vcf_path, out_dir, sample, num_threads)