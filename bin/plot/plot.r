.libPaths(c("/home/yuanw2/R/scDNA-seq","/home/luoc1/R/rlib-4.4.2"))
library(optparse)
library(readxl) 
library(graphics)
library(cowplot)
library(gplots)
library(RColorBrewer)
library(grDevices)
library(ggplot2)
get_script_dir <- function() {
  # Case 1: Rscript execution
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  
  # Case 2: Sourced script
  if (!is.null(parent.frame(2)$ofile)) {
    return(dirname(normalizePath(parent.frame(2)$ofile)))
  }
  
  # Case 3: RStudio interactive
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getActiveDocumentContext()$path)))
  }
  
  # Fallback: current working directory
  return(normalizePath(getwd()))
}

script_dir <- get_script_dir()
print(script_dir)


source(file.path(script_dir, "plot_order.r"))

#!/usr/bin/env Rscript



# Define arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input dir", metavar = "FILE"),
  make_option(c("-o", "--output"), type = "character", default = "output.txt",
              help = "Output dir [default = %default]", metavar = "FILE")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Interactive fallback if arguments missing
if (is.null(opt$input)) {
  cat("No input provided. Please enter interactively:\n")
  opt$input <- readline(prompt = "Enter input file path: ")
}

cat("Using input:", opt$input, "\n")
cat("Output will be saved to:", opt$output, "\n")


if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}


# # Main logic
# if (opt$verbose) {
#   cat("Verbose mode is ON\n")
# }





samples = c(paste0("p", c(1, 3, 4, 5)), paste0("S0_section_", c("A", "B", "C", "D", "E")), paste0("TN", c(1:8)), "KTN302", "T10", "T16")
paths =  paste0(opt$input, samples, "/CNV_", samples, ".csv.T")


Cluster_number=c(3, 4, 5, 4, 3, 4, 3, 2)

reorder_cells <- function(cell_vector, reference, CNV) { 
  match_indices <- match(cell_vector, reference$cell_number) 
  match_indices <- na.omit(match_indices) 
  ordered_CNV <- CNV[order(match_indices),]
  return(ordered_CNV)
}
reorder_srr_ids <- function(srr_ids, other_vector) {
  other_positions <- setNames(seq_along(other_vector), other_vector)
  new_order <- srr_ids    
  filled_positions <- rep(FALSE, length(srr_ids))
  
  for (id in other_vector) {
    if (id %in% srr_ids) {
      position_in_srr_ids <- match(id, srr_ids)
      position_in_other_vector <- other_positions[[id]]
      new_order[position_in_other_vector] <- id
      filled_positions[position_in_srr_ids] <- TRUE
    }
  }
  remaining_ids <- srr_ids[!filled_positions]
  new_order[!filled_positions] <- remaining_ids
  
  return(new_order)
}
reorder_indices <- function(srr_ids, other_vector) {
  other_positions <- setNames(seq_along(other_vector), other_vector)
  order_indices <- seq_along(srr_ids)
  reassigned <- rep(FALSE, length(srr_ids))
  for (id in other_vector) {
    if (id %in% srr_ids) {
      current_index <- match(id, srr_ids)
      target_index <- other_positions[[id]]
      order_indices[target_index] <- current_index
      reassigned[current_index] <- TRUE
    }
  }
  
  remaining_indices <- order_indices[!reassigned]
  order_indices[!reassigned] <- remaining_indices
  
  return(order_indices)
}

for (i in 1:4){
    
    test=read.csv(paths[i])
    bin=test[,1] 
    cnv=as.matrix(test[,-1])
    ref=sub(":.*", "",  test$bin)
    cell_names = colnames(cnv)

  data_tumor <- read.table(paste0("/data/maiziezhou_lab/Weiman/Simulator/Cell_node/cell_node", c(1, 3, 4, 5)[i]), header = FALSE, stringsAsFactors = FALSE)
  data_normal <-read.table("/data/maiziezhou_lab/Weiman/Simulator/Cell_node/cell_node2", header = FALSE, stringsAsFactors = FALSE)
  data_tumor$cell_number <- as.numeric(sub("cell", "", data_tumor$V1))
  data_tumor$node_number <- as.numeric(sub("node", "", data_tumor$V2)) 
  data_tumor <-  data_tumor[order(data_tumor$node_number, data_tumor$cell_number), ]
  data_normal$cell_number <- as.numeric(sub("normalcell", "", data_normal$V1))
  data_normal$node_number <- rep(1, nrow(data_normal)) 
  data_normal <-  data_normal[order(data_normal$node_number, data_normal$cell_number), ]
  normal_cells <- grepl("^normal", cell_names)
  normal_cell_names <- cell_names[normal_cells]
  other_cell_names <- cell_names[!normal_cells]
  normal_cell_numbers <- as.numeric(gsub("normalcell(\\d+)", "\\1", normal_cell_names))
  other_cell_numbers <- as.integer(regmatches(other_cell_names, regexpr("(?<=_cell)\\d+", other_cell_names, perl = TRUE)))

  normal_CNV <- cnv[, normal_cells]
  tumor_CNV <- cnv[,!normal_cells]

  normal_CNV <- reorder_cells(normal_cell_numbers, data_normal, t(normal_CNV))
  tumor_CNV <- reorder_cells(other_cell_numbers, data_tumor, t(tumor_CNV))
  cn_matrix = rbind(tumor_CNV, normal_CNV)
  Gini <- runif(nrow(cn_matrix))
 
  plot.test <- plotiCN(t(cn_matrix), ref, Gini, paste0(opt$output,"/CNVeil_", samples[i], ".png"), annotation = NULL,
                  plot.dendrogram = FALSE, show.names = FALSE)  
    
} ##simulation, gold standard

for (i in 5:17){
    
    test=read.csv(paths[i])
    bin=test[,1] 
    cnv=as.matrix(test[,-1])
    ref=sub(":.*", "",  test$bin)
    cell_names = colnames(cnv)
    cell_names = sub("\\.1", "", cell_names)
    Gini <- runif(ncol(cnv))

    gold_cell=readRDS(paste0("/data/maiziezhou_lab/Weiman/CNVeil/aneufinder/", samples[i], "/cell_names.rds")) 
    gold_cell= sub("\\-1.bam", "", gold_cell)
    gold_order=readRDS(paste0("/data/maiziezhou_lab/Weiman/CNVeil/aneufinder/AF_", samples[i], ".pngrclust_order.rds")) 

    
    ordered_cell <- cell_names[match(cell_names, gold_cell)][gold_order]
    matched_indices <- match(gold_cell[gold_order], cell_names) 
    valid_indices <- matched_indices[!is.na(matched_indices)] 
     
    cnv <- cnv[, valid_indices] 
    plot.test<- plotiCN(cnv, ref, Gini, paste0(opt$output,"/CNVeil_", samples[i], ".png"), annotation = NULL,
                  plot.dendrogram = FALSE, show.names = FALSE)  
    
} #10X & ACT, by AneuFinder




for (i in 18){
    
    test=read.csv(paths[i])
    bin=test[,1] 
    cnv=as.matrix(test[,-1])
    ref=sub(":.*", "",  test$bin)
    cell_names = colnames(cnv)
    cell_names = sub("\\.1", "", cell_names)
    Gini <- runif(ncol(cnv))

    pre_names = read.csv("/data/maiziezhou_lab/Alumni/Zihang/CNV/code/yunfei_python/Data/TNBC_DNA_bowtie2_resolution_500k/treatment_pre.csv", 
                     row.names = 1, header = 1)
    index = cell_names %in% pre_names$x
    pre = cnv[, index]
    mid = cnv[, !index]
    cnv = cbind(pre, mid) 
        plot.test<- plotiCN(cnv, ref, Gini, paste0(opt$output,"/CNVeil_", samples[i], ".png"), annotation = NULL,
                      plot.dendrogram = FALSE, show.names = FALSE)  
        
} #KTN302

for (i in 19){
    test=read.csv(paths[i])
    bin=test[,1] 
    cnv=as.matrix(test[,-1])
    ref=sub(":.*", "",  test$bin)
    cell_names = colnames(cnv)
    cell_names = sub("\\.1", "", cell_names) 

    seg = readRDS("/data/maiziezhou_lab/Weiman/Newdataset/SCOPE/draft/getfgc_group/segment_cs.rds") 
    SCOPE = seg[[1]]$iCN
    for (j in seq(2,22)){
        SCOPE_j = seg[[j]]$iCN
        SCOPE= rbind(SCOPE, SCOPE_j)
        }
    char_vector <- cell_names
    other_vector <- dimnames(SCOPE)[[2]]
    pattern <- "SRR\\d+"
    srr_ids <- regmatches(char_vector, gregexpr(pattern, char_vector))
    srr_ids <- unlist(srr_ids)

    reordered_srr_ids <- reorder_srr_ids(srr_ids, other_vector)
    reordered_indices <- reorder_indices(srr_ids, other_vector)
    # AneuFinder.CNV <- AneuFinder.CNV[, reordered_indices]
    gold = read.delim("/data/maiziezhou_lab/Weiman/CNVeil/cell2sector", header = FALSE, stringsAsFactors = FALSE) 
    split_data <- strsplit(gold$V1, "\\s+")
    first_elements <- sapply(split_data, `[`, 1)
    cnv <- cnv[, match(first_elements, srr_ids)]
    Gini <- runif(ncol(cnv))
    plot.test<- plotiCN(cnv, ref, Gini, paste0(opt$output,"/CNVeil_", samples[i], ".png"), annotation = NULL,
                  plot.dendrogram = FALSE, show.names = FALSE)  
    
} #T10


for (i in 20){
    
    test=read.csv(paths[i])
    bin=test[,1] 
    cnv=as.matrix(test[,-1])
    ref=sub(":.*", "",  test$bin)
    cell_names = colnames(cnv)
    cell_names = sub("\\.1", "", cell_names)
    Gini <- runif(ncol(cnv))

    # gold_cell=readRDS(paste0("/data/maiziezhou_lab/Weiman/CNVeil/aneufinder/", samples[i], "/cell_names.rds")) 
    # gold_cell= sub("\\-1.bam", "", gold_cell)
    # gold_order=readRDS(paste0("/data/maiziezhou_lab/Weiman/CNVeil/aneufinder/AF_", samples[i], ".pngrclust_order.rds")) 

     qcObj = readRDS("/data/maiziezhou_lab/Weiman/Newdataset/SCOPE/T16/qcObj.rds")
     gold_cell = qcObj$sampname
     gold_order= readRDS("/data/maiziezhou_lab/Weiman/Newdataset/SCOPE/T16/r_order.rds")

    ordered_cell <- cell_names[match(cell_names, gold_cell)][gold_order]
    matched_indices <- match(gold_cell[gold_order], cell_names) 
    valid_indices <- matched_indices[!is.na(matched_indices)] 

    cnv <- cnv[, valid_indices] 
    plot.test<- plotiCN(cnv, ref, Gini, paste0(opt$output, "/CNVeil_", samples[i], ".png"), annotation = NULL,
                  plot.dendrogram = FALSE, show.names = FALSE)  
    
} #T16