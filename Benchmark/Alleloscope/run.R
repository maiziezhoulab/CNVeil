.libPaths(c("~/R/Alleloscope", .libPaths()))
library(Alleloscope)
library(DescTools)
library(data.table)
library(Matrix)
source("Cbn_matrix.R")

samplename <- "TN1"
vcf_path <- "/path/to/merged_single"
coverage_csv <- "/path/to/secnv/TN1/cnv_matrix.csv"
normal_rds <- "/path/to/normalRD.rds"
barcodes_path <- "/path/to/single_cell_barcodes.tsv"
size_file <- "/path/to/sizes.cellranger-atac-hg19-1.2.0.txt"
checkpoint <- "/path/to/Obj_filtered"
outputfolder <- "/path/to/Alleloscope/TN1"
dir.create(outputfolder, recursive = TRUE, showWarnings = FALSE)

# Optional matrix combination, as in the original commented preprocessing.
# Cbn_matrix(dir_path = vcf_path, vcf_path = vcf_path,
#            barcodes_path = barcodes_path, samplename = samplename)
alt_all <- readRDS(file.path(vcf_path, "alt_all.rds"))
ref_all <- readRDS(file.path(vcf_path, "ref_all.rds"))
var_all <- readRDS(file.path(vcf_path, "var_all.rds"))
barcodes <- read.table(barcodes_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
test <- read.csv(coverage_csv)
cell_names <- test[, 1]
bin_names <- colnames(test)[-1]
# Equivalent numeric bin-by-cell matrix; avoid coercing the cell-name column.
matrix_data <- t(as.matrix(test[, -1, drop = FALSE]))
cell_names <- gsub("^(SRR[0-9]+).*", "\\1", cell_names)
bin_names <- gsub("^(chr[^.]+)\\.([0-9]+)\\.([0-9]+)$", "\\1-\\2-\\3", bin_names)
rownames(matrix_data) <- bin_names
colnames(matrix_data) <- cell_names
ref_counts <- as.matrix(readRDS(normal_rds))
rownames(ref_counts) <- rownames(matrix_data)
genome_assembly <- "GRCh37"
size <- read.table(size_file, stringsAsFactors = FALSE)

# Original workflow resumes after Est_regions().
# For a fresh object, prepare consistently aligned cells in all input matrices,
# then run the original object-building steps below instead of readRDS().
# Obj <- Createobj(alt_all = alt_all, ref_all = ref_all, var_all = var_all,
#                  samplename = samplename, genome_assembly = "GRCh37",
#                  dir_path = outputfolder, barcodes = barcodes, size = size,
#                  assay = "scDNAseq")
# Obj_filtered <- Matrix_filter(Obj = Obj, cell_filter = 1000, SNP_filter = 20,
#                               min_vaf = 0.1, max_vaf = 0.9)
# Obj_filtered <- Segmentation(Obj_filtered = Obj_filtered, raw_counts = matrix_data,
#                              ref_counts = ref_counts, plot_seg = TRUE)
# Obj_filtered <- Segments_filter(Obj_filtered = Obj_filtered, nSNP = 2000)
# Obj_filtered <- Est_regions(Obj_filtered = Obj_filtered, max_nSNP = 30000,
#                             plot_stat = TRUE, cont = TRUE)
Obj_filtered <- readRDS(checkpoint)
setwd(outputfolder)
Obj_filtered <- Select_normal(Obj_filtered = Obj_filtered, raw_counts = matrix_data,
                              plot_theta = TRUE)
# Manual TN1 reference: automatic selection yielded mostly diploid calls.
Obj_filtered$ref <- "8:5e+05"
Obj_filtered <- Genotype_value(Obj_filtered = Obj_filtered, type = "cellline",
                               raw_counts = matrix_data, ref_counts = ref_counts)
saveRDS(Obj_filtered, "Obj_filtered_1")
Obj_filtered$genotype_values[, seq(1, 124, by = 4)] <-
    Obj_filtered$genotype_values[, seq(1, 124, by = 4)] * 2
Obj_filtered <- Genotype(Obj_filtered = Obj_filtered)
linplot <- Lineage_plot(Obj_filtered = Obj_filtered, nSNP = 2000, nclust = 10)
saveRDS(Obj_filtered, "Obj_filtered.rds")
