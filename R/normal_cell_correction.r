library(data.table)
library(dplyr)
library(matrixStats)
library(kernlab)
library(pracma)
library(gtools)
library(SCOPE)
library(DescTools)

work_path <- commandArgs(trailingOnly = TRUE)[1]
normal_cell_file <- commandArgs(trailingOnly = TRUE)[2]
saved_path <- commandArgs(trailingOnly = TRUE)[3]

read_matrix <- function(file_name) {
  f <- file(file_name, "r")
  matrix <- list()
  chr_name <- list()
  bin_list <- list()
  count <- 0
  
  while (length(line <- readLines(f, n = 1)) > 0) {
    count <- count + 1
    line <- strsplit(line, "\t")[[1]]
    
    if (count > 1) {
      matrix[[count - 1]] <- as.integer(line[4:length(line)])
      chr_name[[count - 1]] <- line[1]
      bin_list[[count - 1]] <- paste(line[1], line[2], line[3], sep = ":")
    } else {
      sample_list <- line[4:length(line)]
    }
  }
  
  matrix <- do.call(rbind, matrix)
  chr_name <- unlist(chr_name)
  bin_list <- unlist(bin_list)
  sample_list <- sapply(strsplit(sample_list, "/"), tail, 1)
  chr_name <- factor(chr_name)
  close(f)
  
  return(list(matrix = matrix, chr_name = chr_name, bin_list = bin_list, sample_list = sample_list))
}


# get_gini <- function(Y) {
#   Gini <- rep(NA, ncol(Y))
#   for (i in seq_len(ncol(Y))) {
#     y <- sort(Y[, i])
#     x <- c(0, seq_len(length(y))/length(y))
#     y <- c(0, cumsum(y)/sum(y))
#     Gini[i] <- 2 * round(0.5 - AUC(x, y), 4)
#   }
#   return(Gini)
# }

get_norm_cell <- function(Y, normal_cell_file = "None", quan_vec = 0.4, gini_threshold = 0.12) {
  # pca_result <- prcomp(t(Y), scale. = TRUE) 
  pca_result <- prcomp(t(Y), center = TRUE, scale. = TRUE) 
  sorted_bins <- sort(abs(pca_result$rotation[, 1]), decreasing = TRUE, index.return = TRUE)
  bin_index <- sorted_bins$ix[sorted_bins$x >= quantile(sorted_bins$x, quan_vec)]
  
  Y_selected <- Y[bin_index,] 
  gini_coefficients <- get_gini(Y_selected) 
  
  if (normal_cell_file == "None") { 
    norm_index <- which(gini_coefficients < gini_threshold)
    abnorm_index <- which(gini_coefficients >= gini_threshold)
    print(paste(names(data_set), ": num of normal cell", length(norm_index)))
  } else { 
    print(paste("Using predefined normal cells from", normal_cell_file))
    predefined_normal_cells <- read.table(normal_cell_file, stringsAsFactors = FALSE)[, 1]
    norm_index <- which(read_matrix(data_set)$sample %in% predefined_normal_cells)
    abnorm_index <- which(!(read_matrix(data_set)$sample %in% predefined_normal_cells))
  }
  
  return(list(norm_index = norm_index, abnorm_index = abnorm_index))
}

Bias_norm <- function(Y, norm_cell_index) {
  Y <- t(Y)
  
  if (length(norm_cell_index) > 0) {
    norm_cell_Y <- Y[norm_cell_index, , drop = FALSE]
    bias_matrix <- matrix(0, nrow = nrow(norm_cell_Y), ncol = ncol(norm_cell_Y))
    
    for (i in seq_len(nrow(norm_cell_Y))) {
      cell <- norm_cell_Y[i, ]
      median <- median(cell)
      bias_list <- cell / median
      bias_matrix[i, ] <- bias_list
    }
    
    ave_bias <- colMeans(bias_matrix)
    ave_bias[ave_bias == 0] <- 1
  } else {
    ave_bias <- scan("bias.txt")
  }
  
  gc_nor_Y <- Y / matrix(ave_bias, nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  
  return(t(gc_nor_Y))
}

main <- function() {
 cov_bed_path <- file.path(work_path, "genome_cov.bed")
 result_path <- saved_path
 read_data <- read_matrix(cov_bed_path)
 Y <- read_data$matrix
 chr_name <- read_data$chr_name
 bin_list <- read_data$bin_list
 sample_list <- read_data$sample_list

 var_norm_cell <- get_norm_cell(Y, normal_cell_file = "None", quan_vec = 0.4, gini_threshold = 0.12)
 cov_matrix <- Bias_norm(Y, var_norm_cell$norm_index)
 saveRDS(cov_matrix, paste0(result_path, "RDmatrix.rds"))
}
 
main()
