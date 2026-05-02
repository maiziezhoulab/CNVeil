args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript run_aneufinder.R <merged_bam> <split_bam_folder> <refversion: hg19|hg38> <outputfolder> [numCPU]")
}

bamfile      <- args[1]
datafolder   <- args[2]
refversion   <- args[3]
outputfolder <- args[4]

numCPU <- if (length(args) >= 5) as.integer(args[5]) else 1

message("Using numCPU = ", numCPU)

if (!dir.exists(outputfolder)) {
  dir.create(outputfolder, recursive = TRUE)
}

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
## Load packages
## -----------------------------
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(Rsamtools)
  library(GenomicAlignments)
  library(rtracklayer)
  library(foreach)
  library(doParallel)
})

if (refversion == "hg19") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  GC.BSgenome <- BSgenome.Hsapiens.UCSC.hg19
} else if (refversion == "hg38") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  GC.BSgenome <- BSgenome.Hsapiens.UCSC.hg38
} else {
  stop("refversion must be either 'hg19' or 'hg38'")
}

## -----------------------------
## Source helper functions
## -----------------------------
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
utils_dir <- file.path(script_dir, "utils")

for (f in list.files(utils_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

## -----------------------------
## Settings
## -----------------------------
binsize <- 500e3
chromosomes <- paste0("chr", 1:22)

blacklist_prefix <- file.path(outputfolder, paste0("p3_blacklist_", refversion))
blacklist_bed <- paste0(blacklist_prefix, ".bed.gz")

## -----------------------------
## Step 1: Generate blacklist
## -----------------------------
bins <- binReads(
  file = bamfile,
  assembly = refversion,
  bamindex = paste0(bamfile, ".bai"),
  binsizes = binsize,
  chromosomes = chromosomes,
  calc.complexity = FALSE
)[[1]]

lcutoff <- quantile(bins$counts, 0.05, na.rm = TRUE)
ucutoff <- quantile(bins$counts, 0.999, na.rm = TRUE)

blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist <- reduce(blacklist)

exportGRanges(
  blacklist,
  filename = blacklist_prefix,
  header = FALSE,
  chromosome.format = "UCSC"
)

## -----------------------------
## Step 2: Run minimal AneuFinder
## -----------------------------
bin_GC_minimal(
  inputfolder = datafolder,
  outputfolder = outputfolder,
  pairedEndReads = FALSE,
  assembly = refversion,
  numCPU = numCPU,
  binsizes = binsize,
  chromosomes = chromosomes,
  blacklist = blacklist_bed,
  correction.method = "GC",
  GC.BSgenome = GC.BSgenome,
  one.file.only = FALSE
)