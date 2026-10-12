.libPaths(c("~/R/scDNA-seq", .libPaths()))
library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg19)

patient <- "TN1"
merged_bam <- "/path/to/all_tumor_cells.bam"
datafolder <- "/path/to/individual_cell_bams"
outputfolder <- "/path/to/AneuFinder/TN1"
dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)

# Blacklist from pooled coverage, using the original thresholds.
bins <- binReads(merged_bam, assembly = "hg19", binsizes = 500e3,
                  chromosomes = paste0("chr", 1:22))[[1]]
lcutoff <- quantile(bins$counts, 0.05)
ucutoff <- quantile(bins$counts, 0.999)
blacklist <- bins[bins$counts <= lcutoff | bins$counts >= ucutoff]
blacklist <- reduce(blacklist)
blacklist.file <- file.path(outputfolder, patient)
exportGRanges(blacklist, filename = blacklist.file, header = FALSE,
               chromosome.format = "UCSC")
Aneufinder(inputfolder = datafolder, outputfolder = outputfolder,
            pairedEndReads = FALSE, assembly = "hg19", numCPU = 1,
            binsizes = 500e3, variable.width.reference = NULL,
            chromosomes = paste0("chr", 1:22),
            blacklist = paste0(blacklist.file, ".bed.gz"),
            correction.method = "GC", GC.BSgenome = BSgenome.Hsapiens.UCSC.hg19,
            refine.breakpoints = FALSE, method = "edivisive")
