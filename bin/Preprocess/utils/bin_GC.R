bin_GC_minimal <- function(
  inputfolder,
  outputfolder,
  pairedEndReads = FALSE,
  assembly = "hg19",
  numCPU = 1,
  binsizes = 500e3,
  stepsizes = binsizes,
  chromosomes = paste0("chr", 1:22),
  blacklist = NULL,
  correction.method = "GC",
  GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19,
  remove.duplicate.reads = TRUE,
  min.mapq = 10,
  use.bamsignals = FALSE,
  one.file.only = TRUE
) {
  binpath.uncorrected <- file.path(outputfolder, "binned")
  binpath.corrected <- paste0(binpath.uncorrected, "-GC")

  dir.create(outputfolder, showWarnings = FALSE, recursive = TRUE)
  dir.create(binpath.uncorrected, showWarnings = FALSE, recursive = TRUE)

  files <- list.files(
    inputfolder,
    full.names = TRUE,
    pattern = "\\.bam$|\\.bed$|\\.bed\\.gz$"
  )

  if (length(files) == 0) {
    stop("No BAM/BED files found in: ", inputfolder)
  }

  files <- sort(files)

  if (one.file.only) {
    files <- files[1]
    message("Running only one file: ", files)
  }

  bamfile <- grep("\\.bam$", files, value = TRUE)[1]

  if (!is.na(bamfile)) {
    chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
  } else {
    chrom.info <- GenomeInfoDb::getChromInfoFromUCSC(assembly)
    chrom.lengths <- chrom.info$size
    names(chrom.lengths) <- chrom.info$chrom
  }

  chrom.lengths <- chrom.lengths[chromosomes]
  chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]

  chrom.lengths.df <- data.frame(
    chromosome = names(chrom.lengths),
    length = as.numeric(chrom.lengths)
  )

  write.table(
    chrom.lengths.df,
    file = file.path(outputfolder, "chrominfo.tsv"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  bins <- fixedWidthBins(
    chrom.lengths = chrom.lengths,
    chromosomes = chromosomes,
    binsizes = binsizes,
    stepsizes = stepsizes
  )

  pattern <- paste0(
    "binsize_",
    format(binsizes, scientific = TRUE, trim = TRUE),
    "_stepsize_",
    format(stepsizes, scientific = TRUE, trim = TRUE)
  )

  for (file in files) {
    message("Binning: ", file)

    binReads(
      file = file,
      assembly = chrom.lengths.df,
      pairedEndReads = pairedEndReads,
      binsizes = NULL,
      variable.width.reference = NULL,
      reads.per.bin = NULL,
      bins = bins[pattern],
      chromosomes = chromosomes,
      remove.duplicate.reads = remove.duplicate.reads,
      min.mapq = min.mapq,
      blacklist = blacklist,
      outputfolder.binned = binpath.uncorrected,
      save.as.RData = TRUE,
      reads.store = FALSE,
      use.bamsignals = use.bamsignals,
      calc.complexity = FALSE
    )
  }

  if (!is.null(correction.method) && "GC" %in% correction.method) {
    dir.create(binpath.corrected, showWarnings = FALSE, recursive = TRUE)

    binfiles <- list.files(
      binpath.uncorrected,
      pattern = "\\.RData$",
      full.names = TRUE
    )

    binfiles <- sort(binfiles)

    if (one.file.only) {
      binfiles <- binfiles[1]
      message("Correcting only one binned file: ", binfiles)
    }

    if (length(binfiles) == 0 ) {
      stop("No binned .RData files found for GC correction.")
    }

    message("Running GC correction...")

    binned.data.list <- suppressMessages(
      correctGC(
        binfiles,
        GC.BSgenome,
        same.binsize = TRUE
      )
    )

    for (i in seq_along(binned.data.list)) {
      binned.data <- binned.data.list[[i]]
      savename <- file.path(
        binpath.corrected,
        basename(names(binned.data.list)[i])
      )
      save(binned.data, file = savename)
    }
  }

  message("Done.")
  invisible(NULL)
}