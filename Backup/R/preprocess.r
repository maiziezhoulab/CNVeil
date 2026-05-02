library(data.table)
library(GenomicRanges)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SCOPE)
library(WGSmapp)
 
work_path <- commandArgs(trailingOnly = TRUE)[1]
ref <- commandArgs(trailingOnly = TRUE)[2]
bam_dir <- list.files(commandArgs(trailingOnly = TRUE)[3], commandArgs(trailingOnly = TRUE)[4])
bam_file <- list.files(commandArgs(trailingOnly = TRUE)[3], commandArgs(trailingOnly = TRUE)[4], full.names = TRUE)
bin_size <- as.integer(commandArgs(trailingOnly = TRUE)[5])

construct_bins <- function(ref, bin_size, work_path) {
  cat("Constructing genome-wide consecutive bins...\n")
  genome_size <- file.path(work_path, "genome_size.txt")
  genome_consecutive_bins <- file.path(work_path, "genome_consecutive_bins.bed")
  genome_consecutive_bins_add <- file.path(work_path, "genome_consecutive_bins_add.bed")
  
  if (ref == "hg19") {
    genome <- BSgenome.Hsapiens.UCSC.hg19
    genomesize <- as.character(seqlengths(genome)[1:22])
    writeLines(paste(chosen_chr, genomesize), con = genome_size)
    
  } else if (ref == "hg38") {
    genome <- BSgenome.Hsapiens.UCSC.hg38
    genomesize <- as.character(seqlengths(genome)[1:22])
    writeLines(paste(chosen_chr, genomesize), con = genome_size)
  }
  
  f_out <- file(genome_consecutive_bins, "w")
  for (chr_num in 1:22) {
    options(scipen = 999)
    length <- as.integer(genomesize[chr_num])
    start <- 1
    
    if (paste0("chr",chr_num) %in% chosen_chr) {
      while ((length - start) > (bin_size - 1)) {
        options(scipen = 999)
        writeLines(paste(chr_num, start, start + bin_size - 1), f_out)
        start <- start + bin_size
      }
      writeLines(paste(chr_num, start, length), f_out)
    }
  }
  close(f_out)
  
  f_out <- file(genome_consecutive_bins_add, "w")
  count <- 1
  
  for (chr_num in 1:22) {
    length <- as.integer(genomesize[chr_num])
    start <- 1
    
    if (paste0("chr",chr_num) %in% chosen_chr) {
      while ((length - start) > (bin_size - 1)) {
        writeLines(paste(chr_num, start, start + bin_size - 1, count), f_out)
        start <- start + bin_size
        count <- count + 1
      }
      writeLines(paste(chr_num, start, length, count), f_out)
      count <- count + 1
    }
  }
  close(f_out)
}

cal_gc <- function(bed_file) {
  cat("Calculating GC content...\n")
  genome_gc <- file.path(work_path, "genome_gc.bed")
  gc <- get_gc(ref_raw, hgref = ref)
  gc <- as.character(gc)
  pos <- read.table(paste0(work_path,"genome_consecutive_bins.bed"), header = 0)
  writeLines(paste(pos$V1, pos$V2, pos$V3, gc), con = genome_gc)
}

cal_map <- function(bed_file) {
  cat("Calculating mappability...\n")
  genome_map <- file.path(work_path, "genome_mappability.tab")
  mapp <- get_mapp(ref_raw, hgref = ref)
  mapp <- as.character(mapp)
  writeLines(mapp, con = genome_map)
}

perform_qc <- function(Y_raw, sampname_raw, ref_raw, QCmetric_raw,
                       cov_thresh = 0, minCountQC = 20, 
                       mapq20_thresh = 0.3, mapp_thresh = 0.9,
                       gc_thresh = c(20, 80), nMAD = 3) {
  if (length(ref_raw) != nrow(Y_raw)) {
    stop("Invalid inputs: length of ref and # of rows
            in read count matrix must be the same")
  }
  if (length(sampname_raw) != ncol(Y_raw)) {
    stop("Invalid inputs: length of sample names and # of cols in
            read count matrix must be the same")
  }
  if (nrow(QCmetric_raw) != ncol(Y_raw)) {
    stop("Invalid inputs: # of rows in QC metric and # of cols in
            read count matrix must be the same")
  }
  mapp <- ref_raw$mapp
  gc <- ref_raw$gc
  sampfilter1 <- (apply(Y_raw, 2, sum) <= cov_thresh)
  message("Removed ", sum(sampfilter1),
          " samples due to failed library preparation.")
  sampfilter2 <- (apply(Y_raw, 2, mean) <= minCountQC)
  message("Removed ", sum(sampfilter2),
          " samples due to failure to meet min coverage requirement.")    
  sampfilter3 <- (QCmetric_raw[, "mapped_prop"] < mapq20_thresh)
  message("Removed ", sum(sampfilter3),
          " samples due to low proportion of mapped reads.")
  if (sum(sampfilter1 | sampfilter2 | sampfilter3) != 0) {
    Y <- Y_raw[, !(sampfilter1 | sampfilter2 | sampfilter3)]
    sampname <- sampname_raw[!(sampfilter1 | sampfilter2 
                               | sampfilter3)]
    QCmetric <- QCmetric_raw[!(sampfilter1 | sampfilter2 
                               | sampfilter3), ]
  } else {
    Y <- Y_raw
    sampname <- sampname_raw
    QCmetric <- QCmetric_raw
  }
  binfilter1 <- (gc < gc_thresh[1] | gc > gc_thresh[2])
  message("Excluded ", sum(binfilter1),
          " bins due to extreme GC content.")
  binfilter2 <- (mapp < mapp_thresh)
  message("Excluded ", sum(binfilter2),
          " bins due to low mappability.")
  if (sum(binfilter1 | binfilter2) != 0) {
    ref <- ref_raw[!(binfilter1 | binfilter2)]
    Y <- Y[!(binfilter1 | binfilter2), ]
  } else {
    ref <- ref_raw
    Y <- Y
  }
  Y.nonzero <- Y[apply(Y, 1, function(x) {
    !any(x == 0)
  }), , drop = FALSE]
  if(dim(Y.nonzero)[1] <= 10){
    message("Adopt arithmetic mean instead of geometric mean")
    pseudo.sample <- apply(Y, 1, mean)
    N <- apply(apply(Y, 2, function(x) {
      x/pseudo.sample
    }), 2, median, na.rm = TRUE)
  } else{
    pseudo.sample <- apply(Y.nonzero, 1, function(x) {
      exp(sum(log(x))/length(x))
    })
    N <- apply(apply(Y.nonzero, 2, function(x) {
      x/pseudo.sample
    }), 2, median)
  }
  sampfilter3 <- (N == 0)
  message("Removed ", sum(sampfilter3),
          " samples due to excessive zero read counts in 
            library size calculation.")
  if (sum(sampfilter3) != 0) {
    Y <- Y[, !(sampfilter3)]
    sampname <- sampname[!(sampfilter3)]
    QCmetric <- QCmetric[!(sampfilter3), ]
    N <- N[!(sampfilter3)]
  }
  Nmat <- matrix(nrow = nrow(Y), ncol = ncol(Y), data = N,
                 byrow = TRUE)
  bin.sum <- apply(Y/Nmat, 1, sum)
  binfilter3 <- (bin.sum >= (median(bin.sum) -
                               nMAD * mad(bin.sum))) & (bin.sum <= (median(bin.sum) +
                                                                      nMAD * mad(bin.sum)))
  Y <- Y[binfilter3, ]
  ref <- ref[binfilter3]
  QCmetric <- as.data.frame(QCmetric)
  message("There are ", ncol(Y), " samples and ",
          nrow(Y), " bins after QC step. ")
  list(Y = Y, sampname = sampname, ref = ref, QCmetric = QCmetric)
}

cal_cov <- function(bed_file, bam_files) {
  read_num <- file.path(work_path, "read_num.txt")
  sample_list_f <- file.path(work_path, "sample_list.txt")
  
  if (file.exists(read_num)) {
    file.remove(read_num)
  }
  
  writeLines(all_files, con = sample_list_f)
  sample_list <- readLines(sample_list_f)
  read_num_list <- integer(length(sample_list))
  
  cat("Calculating coverage...\n") 
  coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 0, seq = 'single-end', hgref = "hg38")
  Y_raw <- coverageObj$Y
  separated_raw <- lapply(dimnames(Y_raw)[[1]], function(x) {
    parts <- unlist(strsplit(x, "[:-]"))
    data.frame(chr = parts[1], start = parts[2], end = parts[3], stringsAsFactors = FALSE)
  })
   
  result_raw <- do.call(rbind, separated_raw) 
  genome_cov_raw_file <- file.path(work_path, "genome_cov_raw.bed") 
  combined_data_raw <- cbind(result_raw, Y_raw) 
  write.table(combined_data_raw, file = genome_cov_raw_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  QCmetric_raw <- get_samp_QC(bambedObj)
  gc_data <- read.table(paste0(work_path,"genome_gc.bed"), header = 0)
  mapp_data <- read.table(paste0(work_path,"genome_mappability.tab"), header = 0)
  ref_raw$gc <- gc_data$V4
  ref_raw$mapp <- map_data$V1
  qcObj <- perform_qc(Y_raw = Y_raw,
                      sampname_raw = sampname_raw, ref_raw = ref_raw,
                      QCmetric_raw = QCmetric_raw, cov_thresh = 0, minCountQC = 0,
                      mapq20_thresh = 0.2, mapp_thresh = 0.9,
                      gc_thresh = c(20,80), nMAD = 3)
  Y <- qcObj$Y
  separated <- lapply(dimnames(Y)[[1]], function(x) {
    parts <- unlist(strsplit(x, "[:-]"))
    data.frame(chr = parts[1], start = parts[2], end = parts[3], stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, separated)
  
  genome_cov_file <- file.path(work_path, "genome_cov.bed")
  
  combined_data <- cbind(result, Y)
  
  write.table(combined_data, file = genome_cov_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}

main <- function() {
  sampname_raw <- sapply(strsplit(bam_files, ".", fixed = TRUE), "[", 1)
  bamdir <- file.path(bamfolder, bam_files)
  bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                         hgref = ref, resolution = 500)
  ref_raw <- bambedObj$ref
  construct_bins(ref, bin_size)
  cal_gc(file.path(work_path, "genome_consecutive_bins.bed"))
  cal_map(file.path(work_path, "genome_consecutive_bins_add.bed"))
  cal_cov(file.path(work_path, "ref.bed"), bam_files)
}

main()
