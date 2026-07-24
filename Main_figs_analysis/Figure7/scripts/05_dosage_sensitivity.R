#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(DESeq2)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

DNA_CN_FILE <- "example_data/figure7/transcriptomics/dna_cnv_matched.csv"
RNA_COUNT_FILE <- "example_data/figure7/transcriptomics/gene_count_matrix.tsv"
CELL_MAP_FILE <- "example_data/figure7/transcriptomics/cell_id_map.csv"
SUBCLONE_FILE <- "example_data/figure7/copy_number/haplotype_cluster_assignments.tsv"
GENE_ANNOTATION_FILE <- "example_data/figure7/gencode_v38_gene_pos.symbol.txt"
OUTPUT_DIR <- "results/figure7"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS <- 20
SENSITIVE_GENES <- c("MTDH", "IGF1R", "CDH1")
INSENSITIVE_GENES <- c("MDM2", "MAP3K1", "SHC1")

dna_cn <- fread(DNA_CN_FILE, data.table = FALSE)
rownames(dna_cn) <- dna_cn[[1]]
dna_cn <- dna_cn[, -1, drop = FALSE]

counts <- fread(RNA_COUNT_FILE, data.table = FALSE)
rownames(counts) <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])

cell_map <- fread(CELL_MAP_FILE, data.table = FALSE)
subclone_df <- fread(SUBCLONE_FILE, data.table = FALSE)
cell_col <- intersect(c("Cell_ID", "Cell", "cell", "cell_id"), colnames(subclone_df))[1]
cluster_col <- intersect(c("Cluster", "cluster", "cluster_id"), colnames(subclone_df))[1]
subclone_labels <- setNames(as.character(subclone_df[[cluster_col]]), as.character(subclone_df[[cell_col]]))

shared_barcodes <- intersect(colnames(counts), cell_map$rna_barcode)
counts <- counts[, shared_barcodes, drop = FALSE]
barcode_map <- setNames(cell_map$srr_id, cell_map$rna_barcode)
colnames(counts) <- barcode_map[colnames(counts)]

shared_cells <- Reduce(intersect, list(colnames(counts), rownames(dna_cn), names(subclone_labels)))
counts <- counts[, shared_cells, drop = FALSE]
dna_cn <- dna_cn[shared_cells, , drop = FALSE]
subclone_labels <- subclone_labels[shared_cells]

cluster_sizes <- table(subclone_labels)
valid_clusters <- names(cluster_sizes[cluster_sizes >= MIN_CELLS])
keep_cells <- names(subclone_labels)[subclone_labels %in% valid_clusters]
counts <- counts[, keep_cells, drop = FALSE]
dna_cn <- dna_cn[keep_cells, , drop = FALSE]
subclone_labels <- subclone_labels[keep_cells]

seurat_obj <- CreateSeuratObject(counts = counts)
seurat_obj$subclone <- subclone_labels[colnames(seurat_obj)]
Idents(seurat_obj) <- "subclone"
pseudo_counts <- AggregateExpression(seurat_obj, assays = "RNA", group.by = "subclone", return.seurat = FALSE)$RNA
pseudo_counts <- as.matrix(pseudo_counts)
storage.mode(pseudo_counts) <- "integer"
pseudo_counts <- pseudo_counts[rowSums(pseudo_counts) > 0, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = pseudo_counts,
  colData = data.frame(subclone = colnames(pseudo_counts), row.names = colnames(pseudo_counts)),
  design = ~1
)
dds <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))
colnames(vst_mat) <- sub("^g", "", colnames(vst_mat))

cn_consensus <- do.call(cbind, lapply(valid_clusters, function(cl) {
  cells <- names(subclone_labels)[subclone_labels == cl]
  round(colMeans(dna_cn[cells, , drop = FALSE]))
}))
colnames(cn_consensus) <- valid_clusters
rownames(cn_consensus) <- colnames(dna_cn)

parse_bins <- function(bin_names) {
  parts <- strsplit(bin_names, ":|-")
  do.call(rbind, lapply(parts, function(x) {
    data.frame(bin = paste(x, collapse = ":"), chrom = x[1], start = as.integer(x[2]), end = as.integer(x[3]))
  }))
}

bin_df <- parse_bins(rownames(cn_consensus))
bin_df$bin <- rownames(cn_consensus)
genes <- fread(GENE_ANNOTATION_FILE, header = FALSE, data.table = FALSE)
colnames(genes) <- c("gene", "chrom", "start", "end")
if (grepl("^chr", bin_df$chrom[1]) && !grepl("^chr", genes$chrom[1])) genes$chrom <- paste0("chr", genes$chrom)
if (!grepl("^chr", bin_df$chrom[1]) && grepl("^chr", genes$chrom[1])) genes$chrom <- sub("^chr", "", genes$chrom)

gr_bins <- GRanges(seqnames = bin_df$chrom, ranges = IRanges(bin_df$start, bin_df$end), bin = bin_df$bin)
gr_genes <- GRanges(seqnames = genes$chrom, ranges = IRanges(genes$start, genes$end), gene = genes$gene)
hits <- findOverlaps(gr_genes, gr_bins)
map <- data.frame(gene = gr_genes$gene[queryHits(hits)], bin = gr_bins$bin[subjectHits(hits)])
map <- map[!duplicated(map$gene), ]
cn_genes <- cn_consensus[map$bin, , drop = FALSE]
rownames(cn_genes) <- map$gene

common_clusters <- intersect(colnames(vst_mat), colnames(cn_genes))

plot_gene <- function(gene) {
  expression <- as.numeric(vst_mat[gene, common_clusters])
  copy_number <- as.numeric(cn_genes[gene, common_clusters])
  test <- cor.test(copy_number, expression, method = "pearson")
  df <- data.frame(subclone = common_clusters, copy_number = copy_number, expression = expression)
  ggplot(df, aes(x = factor(copy_number), y = expression)) +
    stat_summary(fun = mean, geom = "col", alpha = 0.55) +
    geom_point(aes(color = factor(subclone)), position = position_jitter(width = 0.08), size = 2) +
    geom_line(
      data = df %>% group_by(copy_number) %>% summarize(expression = mean(expression), .groups = "drop"),
      aes(x = factor(copy_number), y = expression, group = 1), color = "black"
    ) +
    labs(
      title = sprintf("%s: r=%.2f, p=%.2g", gene, unname(test$estimate), test$p.value),
      x = "Integer copy number",
      y = "Mean subclonal gene expression"
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
}

sensitive <- Filter(Negate(is.null), lapply(SENSITIVE_GENES, function(g) if (g %in% rownames(vst_mat) && g %in% rownames(cn_genes)) plot_gene(g)))
insensitive <- Filter(Negate(is.null), lapply(INSENSITIVE_GENES, function(g) if (g %in% rownames(vst_mat) && g %in% rownames(cn_genes)) plot_gene(g)))

figure_e <- (wrap_plots(sensitive, nrow = 1) + plot_annotation(title = "Dosage-sensitive genes")) /
  (wrap_plots(insensitive, nrow = 1) + plot_annotation(title = "Dosage-insensitive genes"))

ggsave(file.path(OUTPUT_DIR, "Figure7e_dosage_sensitivity.pdf"), figure_e, width = 13, height = 6)
ggsave(file.path(OUTPUT_DIR, "Figure7e_dosage_sensitivity.png"), figure_e, width = 13, height = 6, dpi = 300)
