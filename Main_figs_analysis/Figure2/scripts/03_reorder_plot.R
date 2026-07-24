.libPaths("~/R/scDNA-seq")

library(graphics)
library(grDevices)

ROOT <- "/data/maiziezhou_lab/Weiman/CNVeil/benchmark_sc/total/cnv_results"
 
TARGET_FOLDERS <- c("T16", "KTN302")
# or full paths:
# TARGET_FOLDERS <- c(paste0("TN", 1:8), )
# TARGET_FOLDERS <- c(paste0("S0_section_", c("A","B","C","D","E")))
# TARGET_FOLDERS <- paste0("S0_section_", c("A","B","C","D","E"))

RECURSIVE_SEARCH <- TRUE  
 
read_cn_table <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!("Cell" %in% colnames(df))) stop("Expected column 'Cell' in: ", path)
  cells <- df$Cell
  df$Cell <- NULL
  mat <- as.matrix(df)
  rownames(mat) <- cells
  storage.mode(mat) <- "numeric"
  mat
}

write_cn_table <- function(mat, path) {
  out <- data.frame(Cell = rownames(mat), mat, check.names = FALSE)
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
 
cnveil_self_order <- function(mat_cells_x_bins) {
  x <- mat_cells_x_bins
  x[x < 0] <- 0
  x[x > 7] <- 7
  hc <- hclust(dist(x))
  rownames(x)[hc$order]
}
 
reorder_by_cnveil <- function(tool_mat, cnveil_order) {
  tool_cells <- rownames(tool_mat)
  intersect_cells <- cnveil_order[cnveil_order %in% tool_cells]
  extra_cells <- tool_cells[!(tool_cells %in% cnveil_order)]
  new_order <- c(intersect_cells, extra_cells)
  tool_mat[new_order, , drop = FALSE]
}
 
list_cn_files <- function(dir) {
  list.files(dir, pattern = "\\.(tsv|txt|csv)$", full.names = TRUE, ignore.case = TRUE)
}

pick_cnveil_file <- function(files) {
  cand <- files[grepl("cnveil", basename(files), ignore.case = TRUE)]
  if (length(cand) == 0) return(NA_character_)
  cand[which.max(file.info(cand)$size)]
}
 
plot_fixed_order <- function(mat_cells_x_bins, outfile_png) {

  smart_image <- function(m, ...) {
    image(t(m[rev(seq_len(nrow(m))), , drop = FALSE]), useRaster = TRUE, ...)
  }

  hm_col <- c("#2166AC", "#92C5DE", "#FDFDFD", "#FDDBC7",
              "#F4A582", "#D6604D", "#B2182B", "#67001F")

  dat <- round(mat_cells_x_bins)
  dat[dat < 0] <- 0
  dat[dat > 7] <- 7
 
  bins <- colnames(dat)
  chr <- sub(":.*$", "", bins)
  chr_levels <- unique(chr)

  chr_len <- vapply(chr_levels, function(x) sum(chr == x), integer(1))
  chr_pos <- cumsum(chr_len)
 
  xpos <- (c(0, chr_pos[-length(chr_pos)]) + chr_pos) / 2 / ncol(dat)
  chr_print <- sub("^chr", "", chr_levels)
 
  band_is_gray <- (seq_along(chr_levels) %% 2 == 1)
 
  # gray band -> black text; black band -> light text (grey/white-ish)
  chr_text_col <- ifelse(band_is_gray, "black", "grey85")
 
  band_h <- 0.032

  png(outfile_png, width = 2500, height = 1500, pointsize = 25, type = "cairo")

  # top band, bottom heatmap
  layout(matrix(c(1, 2), nrow = 2), heights = c(band_h, 1 - band_h))
  par(mar = c(0, 0, 0, 0))

  # ---- 1) chromosome band (TOP)
  anno.chrom <- unlist(mapply(function(i, n) rep(ifelse(i %% 2 == 1, 1, 2), n),
                              seq_along(chr_levels), chr_len))
  anno.chrom <- matrix(anno.chrom, nrow = 1)

  image(t(anno.chrom), col = c("gray", "black"),
        xaxs = "i", yaxs = "i", axes = FALSE, useRaster = TRUE)

  # put text slightly LOWER than 0.5 to avoid "too up" effect in thin band
  # and color alternates with the band (fixes your even chr invisibility)
  text(xpos, 0.2, chr_print, cex = 1, col = chr_text_col)

  # ---- 2) heatmap (BOTTOM)
  par(mar = c(0, 0, 0, 0))

  vals <- sort(unique(as.numeric(dat)))
  vals <- vals[vals >= 0 & vals <= 7]
  if (length(vals) == 0) vals <- 0:7

  smart_image(dat, col = hm_col[vals + 1], xaxs = "i", yaxs = "i", axes = FALSE)

  abline(v = 0, lwd = 2)
  for (p in chr_pos) abline(v = p / ncol(dat), lwd = 2)

  dev.off()
}
 
resolve_targets <- function(root, targets, recursive_search = TRUE) {
  # convert relative -> full path; keep full paths as-is
  tgt_paths <- vapply(targets, function(t) {
    if (grepl("^/", t)) t else file.path(root, t)
  }, character(1))

  tgt_paths <- tgt_paths[file.exists(tgt_paths)]
  if (length(tgt_paths) == 0) return(character(0))

  if (!recursive_search) return(unique(tgt_paths))

  # If recursive: include subdirs too, but still only under user targets
  all <- unique(unlist(lapply(tgt_paths, function(p) c(p, list.dirs(p, recursive = TRUE, full.names = TRUE)))))
  unique(all)
}
 
candidate_dirs <- resolve_targets(ROOT, TARGET_FOLDERS, RECURSIVE_SEARCH)
 
case_dirs <- character(0)
for (d in candidate_dirs) {
  files <- list_cn_files(d)
  if (length(files) == 0) next
  cnv <- pick_cnveil_file(files)
  if (!is.na(cnv)) case_dirs <- c(case_dirs, d)
}
case_dirs <- unique(case_dirs)

cat("[INFO] Will process", length(case_dirs), "case directories.\n")
if (length(case_dirs) > 0) cat(paste0("  ", case_dirs, collapse = "\n"), "\n")

for (case in case_dirs) {

  cat("\n[INFO] Case:", case, "\n")

  files <- list_cn_files(case)
  cnveil_file <- pick_cnveil_file(files)
  if (is.na(cnveil_file)) {
    cat("[WARN] No CNVeil file in:", case, "\n")
    next
  }

  out_dir <- file.path(case, "ordered_and_plotted_by_CNVeil")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cnv <- tryCatch(read_cn_table(cnveil_file), error = function(e) e)
  if (inherits(cnv, "error")) {
    cat("[ERROR] CNVeil read failed:", cnveil_file, "\n  ", cnv$message, "\n")
    next
  }

  cnveil_order <- cnveil_self_order(cnv)
  saveRDS(cnveil_order, file.path(out_dir, "CNVeil.self_order.rds"))
 
  cnv_ord <- cnv[cnveil_order, , drop = FALSE]
  write_cn_table(cnv_ord, file.path(out_dir, "CNVeil.ordered.tsv"))
  plot_fixed_order(cnv_ord, file.path(out_dir, "CNVeil.ordered.png"))
 
  tool_files <- setdiff(files, cnveil_file)
  for (f in tool_files) {
    m <- tryCatch(read_cn_table(f), error = function(e) e)
    if (inherits(m, "error")) {
      cat("  [ERROR] Tool read failed:", f, "\n    ", m$message, "\n")
      next
    }

    m_ord <- reorder_by_cnveil(m, cnveil_order)
    base <- sub("\\.(tsv|txt|csv)$", "", basename(f), ignore.case = TRUE)

    write_cn_table(m_ord, file.path(out_dir, paste0(base, ".ordered.tsv")))
    plot_fixed_order(m_ord, file.path(out_dir, paste0(base, ".ordered.png")))
  }

  cat("[OK] Output:", out_dir, "\n")
}

cat("\n[ALL DONE]\n")
