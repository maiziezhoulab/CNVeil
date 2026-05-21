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


## -----------------------------
## Infer normal cells from total CN
## -----------------------------
infer_normal_cells_from_total_cn <- function(work_dir, sample, tol = 0.1) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))

  cn_file_tsv <- file.path(work_dir, "Total_CN", paste0("CNV_", sample, ".tsv"))
  cn_file_csv <- file.path(work_dir, "Total_CN", paste0("CNV_", sample, ".csv"))

  if (file.exists(cn_file_tsv)) {
    cn_file <- cn_file_tsv
  } else if (file.exists(cn_file_csv)) {
    cn_file <- cn_file_csv
  } else {
    stop("Cannot find CNV file: ", cn_file_tsv, " or ", cn_file_csv)
  }

  message("Reading total CN file for normal-cell inference: ", cn_file)

  cn <- data.table::fread(cn_file)

  if (!"cell" %in% names(cn)) {
    stop("CNV file must have a 'cell' column.")
  }

  cell_ids <- cn[["cell"]]

  cn_cols <- setdiff(names(cn), "cell")
  cn_mat <- as.data.frame(cn[, cn_cols, with = FALSE])
  cn_mat[] <- lapply(cn_mat, as.numeric)

  avg_cn <- rowMeans(cn_mat, na.rm = TRUE)

  normal_cells <- cell_ids[abs(avg_cn - 2) <= tol]

  message("Total cells in CNV file: ", length(cell_ids))
  message("Normal cells inferred: ", length(normal_cells))
  message(
    "Normal-cell percentage: ",
    round(100 * length(normal_cells) / length(cell_ids), 2),
    "%"
  )

  return(normal_cells)
}


## -----------------------------
## Merge SNP-cell ref/alt counts
## -----------------------------
merge_snp_cell_counts <- function(var_all, ref_all, alt_all, out_tsv,
                                  min_total = 1L, integerize = TRUE) {
  stopifnot(requireNamespace("Matrix", quietly = TRUE),
            requireNamespace("data.table", quietly = TRUE))

  if (nrow(var_all) != nrow(ref_all) ||
      nrow(var_all) != nrow(alt_all)) {
    stop("Row mismatch: nrow(var_all) must equal nrow(ref_all) and nrow(alt_all).")
  }

  if (ncol(ref_all) != ncol(alt_all)) {
    stop("Column mismatch: ncol(ref_all) must equal ncol(alt_all).")
  }

  refT <- methods::as(ref_all, "dgTMatrix")
  altT <- methods::as(alt_all, "dgTMatrix")

  sr <- Matrix::summary(refT)
  sa <- Matrix::summary(altT)

  dt_ref <- data.table::as.data.table(sr)
  data.table::setnames(dt_ref, c("i", "j", "x"), c("i", "j", "ref"))

  dt_alt <- data.table::as.data.table(sa)
  data.table::setnames(dt_alt, c("i", "j", "x"), c("i", "j", "alt"))

  dt <- dt_ref[dt_alt, on = c("i", "j")]

  if (!"ref" %in% names(dt)) dt[, ref := NA_real_]
  if (!"alt" %in% names(dt)) dt[, alt := NA_real_]

  dt[is.na(ref), ref := 0]
  dt[is.na(alt), alt := 0]

  if (integerize) {
    dt[, `:=`(
      ref = as.integer(round(ref)),
      alt = as.integer(round(alt))
    )]
  }

  dt[, total := ref + alt]

  if (min_total > 0L) {
    dt <- dt[total >= as.integer(min_total)]
  }

  dt[, total := NULL]

  chrom_vec <- var_all[[1]]
  pos_vec   <- var_all[[2]]
  cell_ids  <- colnames(ref_all)

  dt[, chrom := chrom_vec[i]]
  dt[, pos   := pos_vec[i]]
  dt[, cell  := if (is.null(cell_ids)) as.character(j) else cell_ids[j]]

  out <- dt[, .(chrom, pos, cell, ref, alt)]

  data.table::fwrite(
    out,
    file = out_tsv,
    sep = "\t",
    quote = FALSE
  )

  invisible(out)
}


run_merge_snp_cell_counts <- function(vcf_path, out_dir, sample,
                                      min_total = 1L,
                                      integerize = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  message("Reading RDS files from: ", vcf_path)

  var_file <- file.path(vcf_path, "var_all.rds")
  alt_file <- file.path(vcf_path, "alt_all.rds")
  ref_file <- file.path(vcf_path, "ref_all.rds")

  if (!file.exists(var_file)) stop("Missing file: ", var_file)
  if (!file.exists(alt_file)) stop("Missing file: ", alt_file)
  if (!file.exists(ref_file)) stop("Missing file: ", ref_file)

  var_all <- readRDS(var_file)
  alt_all <- readRDS(alt_file)
  ref_all <- readRDS(ref_file)

  out_file <- file.path(out_dir, paste0(sample, "_baf.tsv"))

  res <- merge_snp_cell_counts(
    var_all = var_all,
    ref_all = ref_all,
    alt_all = alt_all,
    out_tsv = out_file,
    min_total = min_total,
    integerize = integerize
  )

  message("Merged SNP-cell BAF count file written: ", out_file)

  return(res)
}


## -----------------------------
## Split merged BAF file by chromosome
## -----------------------------
split_baf_by_chrom <- function(outfile, outdir, sample, num_threads = 4) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  if (!file.exists(outfile)) {
    stop("Merged BAF file does not exist: ", outfile)
  }

  message("Splitting merged BAF file by chromosome...")

  cmds <- vapply(1:22, function(chr) {
    out_chr <- file.path(outdir, paste0(sample, "_baf_chr", chr, ".tsv"))

    sprintf(
      "grep -w %s %s > %s",
      shQuote(paste0("chr", chr)),
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

  message("Finished splitting BAF file by chromosome.")
}


## -----------------------------
## Compute position-level BAF features
## Exclude inferred normal cells before aggregation
## -----------------------------
compute_pos_baf_one_chrom <- function(in_chr, out_chr_pos, normal_cells = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))

  if (!file.exists(in_chr)) {
    warning("Missing chromosome BAF file: ", in_chr)
    return(invisible(NULL))
  }

  if (file.info(in_chr)$size == 0) {
    warning("Empty chromosome BAF file: ", in_chr)
    return(invisible(NULL))
  }

  dt <- data.table::fread(
    in_chr,
    sep = "\t",
    header = FALSE,
    showProgress = FALSE
  )

  if (nrow(dt) == 0) {
    warning("No rows in chromosome BAF file: ", in_chr)
    return(invisible(NULL))
  }

  if (ncol(dt) < 5) {
    warning("File has fewer than 5 columns: ", in_chr)
    return(invisible(NULL))
  }

  data.table::setnames(dt, 1:5, c("chrom", "pos", "cell", "ref", "alt"))

  ## Important:
  ## Remove normal cells before computing position-level BAF.
  if (!is.null(normal_cells) && length(normal_cells) > 0) {
    n_before <- nrow(dt)

    dt <- dt[!(cell %in% normal_cells)]

    message(
      "Excluded normal cells from ", basename(in_chr),
      ": rows ", n_before, " -> ", nrow(dt)
    )

    if (nrow(dt) == 0) {
      warning("No rows left after excluding normal cells: ", in_chr)
      return(invisible(NULL))
    }
  }

  dt[, pos := as.integer(pos)]
  dt[, ref := as.numeric(ref)]
  dt[, alt := as.numeric(alt)]

  dt[is.na(ref), ref := 0]
  dt[is.na(alt), alt := 0]

  pos_dt <- dt[
    ,
    .(
      ref_sum = sum(ref, na.rm = TRUE),
      alt_sum = sum(alt, na.rm = TRUE),
      n_cells = data.table::uniqueN(cell),
      n_obs = .N
    ),
    by = .(chrom, pos)
  ]

  pos_dt[, depth := ref_sum + alt_sum]
  pos_dt <- pos_dt[depth > 0]

  if (nrow(pos_dt) == 0) {
    warning("No valid position-level BAF rows after depth filter: ", in_chr)
    return(invisible(NULL))
  }

  pos_dt[, baf := alt_sum / depth]
  pos_dt[, mbaf := pmin(baf, 1 - baf)]

  data.table::setorder(pos_dt, chrom, pos)

  data.table::fwrite(
    pos_dt,
    file = out_chr_pos,
    sep = "\t",
    quote = FALSE
  )

  message("Position-level BAF written: ", out_chr_pos)

  invisible(pos_dt)
}


compute_pos_baf_by_chrom <- function(outdir, sample, num_threads = 4,
                                     normal_cells = NULL) {
  stopifnot(requireNamespace("parallel", quietly = TRUE))

  message("Computing position-level BAF features by chromosome...")

  if (!is.null(normal_cells) && length(normal_cells) > 0) {
    message("Normal cells will be excluded during position-level BAF aggregation.")
    message("Number of normal cells to exclude: ", length(normal_cells))
  } else {
    message("No normal cells provided; using all cells for position-level BAF.")
  }

  jobs <- lapply(1:22, function(chr) {
    list(
      chr = chr,
      in_chr = file.path(outdir, paste0(sample, "_baf_chr", chr, ".tsv")),
      out_chr_pos = file.path(outdir, paste0(sample, "_pos_baf_chr", chr, ".tsv.gz"))
    )
  })

  parallel::mclapply(
    jobs,
    function(job) {
      message("Processing chr", job$chr)

      compute_pos_baf_one_chrom(
        in_chr = job$in_chr,
        out_chr_pos = job$out_chr_pos,
        normal_cells = normal_cells
      )
    },
    mc.cores = as.integer(num_threads)
  )

  message("Finished computing position-level BAF features.")
  message("Position-level BAF files written to: ", outdir)

  invisible(NULL)
}


## -----------------------------
## Main wrapper
## -----------------------------
load_snp <- function(vcf_path, work_dir, out_dir, sample, num_threads = 4,
                     min_total = 1L,
                     integerize = TRUE,
                     exclude_normal_cells = TRUE,
                     normal_tol = 0.1) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  message("Running sample: ", sample)
  message("Input RDS folder: ", vcf_path)
  message("Work directory: ", work_dir)
  message("Output folder: ", out_dir)
  message("Threads: ", num_threads)

  ## 1. Merge ref/alt SNP counts across cells
  res <- run_merge_snp_cell_counts(
    vcf_path = vcf_path,
    out_dir = out_dir,
    sample = sample,
    min_total = min_total,
    integerize = integerize
  )

  ## 2. Split merged cell-level BAF count file by chromosome
  baf_file <- file.path(out_dir, paste0(sample, "_baf.tsv"))

  split_baf_by_chrom(
    outfile = baf_file,
    outdir = out_dir,
    sample = sample,
    num_threads = num_threads
  )

  ## 3. Infer normal cells from total CN
  normal_cells <- NULL

  if (exclude_normal_cells) {
    normal_cells <- infer_normal_cells_from_total_cn(
      work_dir = work_dir,
      sample = sample,
      tol = normal_tol
    )
  }

  ## 4. Compute position-level BAF features by chromosome
  ##    using only non-normal cells.
  compute_pos_baf_by_chrom(
    outdir = out_dir,
    sample = sample,
    num_threads = num_threads,
    normal_cells = normal_cells
  )

  message("Done.")
  message("Merged cell-level BAF file: ", baf_file)
  message("Per-chromosome cell-level files: ", out_dir)
  message("Per-chromosome position-level BAF files: ", out_dir)

  invisible(res)
}


## -----------------------------
## CLI
## -----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript load_snp.R <vcf_path> <work_dir> <sample> [num_threads]")
}

vcf_path <- args[1]
work_dir <- args[2]
sample <- args[3]
num_threads <- ifelse(length(args) >= 4, as.integer(args[4]), 4)

out_dir <- file.path(work_dir, "Allele_CN", "allele_count_by_cell")

load_snp(
  vcf_path = vcf_path,
  work_dir = work_dir,
  out_dir = out_dir,
  sample = sample,
  num_threads = num_threads,
  exclude_normal_cells = TRUE,
  normal_tol = 0.1
)