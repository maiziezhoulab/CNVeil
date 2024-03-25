library(optparse)
library(data.table)
library(ggplot2)
library(cluster)
library(stats)

option_list <- list(
  make_option(c('--input_dir', '-i'), type = "character", help = "Input directory should include genome_cov.bed and RDmatrix.rds"),
  make_option(c('--output_dir', '-o'), type = "character"),
  make_option(c('--prefix', '-x'), default = "sample", type = "character"),
  make_option(c('--cell_node', '-c'), type = "character", help = "cell node file; optional, if not given, will not transform cell to node in the end"),
  make_option(c('--change_rate_mode', '-r'), type = "character", default = 'h', help = "q -> 0.84 quantile; h -> 0.1 hardcode"),
  make_option(c('--ploidy_mode', '-p'), type = "character", default = 'cl', help = "global or by cluster")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Validate 'change_rate_mode' and 'ploidy_mode'
valid_change_rate_modes <- c('q', 'h')
valid_ploidy_modes <- c('gl', 'cl')

if (!args$change_rate_mode %in% valid_change_rate_modes) {
  stop("Invalid value for --change_rate_mode. Choose either 'q' or 'h'.")
}

if (!args$ploidy_mode %in% valid_ploidy_modes) {
  stop("Invalid value for --ploidy_mode. Choose either 'gl' or 'cl'.")
}

input_dir <- args$input_dir
output_dir <- args$output_dir
cell_node <- args$cell_node
crm <- args$change_rate_mode
pm <- args$ploidy_mode
prefix <- args$prefix

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE)
}

find_local_minima <- function(numbers) {
  local_minima <- numeric(0)
  n <- length(numbers)
  
  if (n < 3) {
    # If the list has less than 3 elements, it cannot have local minima
    return(local_minima)
  }
  
  idxl <- numeric(0)
  
  for (i in 2:(n - 1)) {
    if (numbers[i] < numbers[i - 1] && numbers[i] < numbers[i + 1]) {
      local_minima <- c(local_minima, numbers[i])
      idxl <- c(idxl, i)
    }
  }
  
  return(list(minima = local_minima, indices = idxl))
}

# Function to get ploidy from one cell CNV
get_ploidy_from_one_cell_cnv <- function(cnv_oc) {
  cand_p <- seq(1.15, 5.50, by = 0.05)
  rs_list <- numeric(length(cand_p))
  
  for (i in seq_along(cand_p)) {
    p <- cand_p[i]
    rs <- sum(abs(cnv_oc * p - round(cnv_oc * p)))
    rs_list[i] <- rs
  }
  min_idx <- which.min(rs_list)
  best_p <- cand_p[min_idx]
  idxl <- find_local_minima(rs_list)$indices
  lm_p <- cand_p[idxl]
  lm_rs <- rs_list[idxl] 
  best_p <- lm_p[which.min(lm_rs)]
  
  return(best_p)
}

# Function to get ploidy for multiple cells
get_ploidy_mt_cell <- function(dat) {
    n <- length(dat)
    m_list <- numeric()
    
    for (i in 1:n) {
        m <- get_ploidy_from_one_cell_cnv(dat[i])
        m_list <- c(m_list, m)
    }
    
    cat(max(m_list) - min(m_list), "\n")
    cat(m_list, "\n")
    
    p_med <- quantile(m_list, 0.5)
    p_mean <- get_ploidy_from_one_cell_cnv(colMeans(matrix(dat, nrow = 1)))
    
    return(p_mean)
}

# Function to calculate BND
BND <- function(cnv_list) {
    which(diff(cnv_list) != 0) + 0.5
}

# Function to calculate nbdiff
nbdiff <- function(cnv_list, flanking) {
  diffl <- numeric(flanking)
  
  for (i in (flanking + 1):(length(cnv_list) - flanking)) {
    lmed <- median(cnv_list[(i - flanking):(i - 1)])
    rmed <- median(cnv_list[i:(i + flanking - 1)])
    diff_val <- abs(lmed - rmed)
    diffl <- c(diffl, diff_val)
  }
  
  diffl <- c(diffl, numeric(flanking))
  return(diffl)
}

# Function to estimate eBND
eBND <- function(cnv_list, flanking, thresh) {
  dc <- nbdiff(cnv_list, flanking)
  bnd <- which(dc >= thresh)
  return(bnd)
}

# Function for hierarchical clustering
cluster <- function(y4, t) {
  cat("---------------------------extract fingerprint bins\n")
  diff <- apply(y4, 2, function(x) max(x) - min(x))
  thresh <- quantile(diff, 0.9)
  
  cat('proportion of finger print bins:', mean(diff >= thresh), "\n")

  cat("---------------------------hierachical clustering\n")
  keep <- which(diff >= thresh)
  data <- y4[, keep, drop = FALSE]
  
  # Using hclust for hierarchical clustering
  hc <- hclust(dist(data), method = "ward.D2")
  labels <- cutree(hc, k = 2)  # Assuming two clusters
  
  y4_by_cluster <- split(y4, labels)
  num_cell_clusters <- length(unique(labels))
  
  cat("num cell clusters: ", num_cell_clusters, "\n")
  cat(table(labels), "\n")
  
  return(as.numeric(labels))
}

# Function to normalize and find tumor/normal cells
get_norm_tm_cell <- function(dat_tnbc, t = 20) {
  cl_tnbc <- cluster(t(dat_tnbc), t)
  normal_cl_idx <- integer()
  tm_cl_idx <- integer()
  
  for (i in unique(cl_tnbc)) {
    p <- get_ploidy_from_one_cell_cnv(rowMeans(dat_tnbc[, cl_tnbc == i, drop = FALSE]))
    cat(p, "\n")
    if (1.7 < p & p < 2.2) {
      normal_cl_idx <- c(normal_cl_idx, which(cl_tnbc == i))
    } else {
      tm_cl_idx <- c(tm_cl_idx, which(cl_tnbc == i))
    }
  }
  
  cat("num normal", length(normal_cl_idx), "\n")
  cat("num tm", length(tm_cl_idx), "\n")
  
  return(list(normal_cl_idx, tm_cl_idx))
}

get_by_cell_cnv <- function(bnd_list, dat, p) {
  last_idx <- 1
  med_dat <- numeric()
  
  for (cidx in bnd_list) {
    med <- median(dat[last_idx:cidx])
    med_dat <- c(med_dat, rep(med, cidx - last_idx + 1))
    last_idx <- cidx + 1
  }

  est_cnv <- as.integer(round(p * med_dat))
  est_bnd <- BND(est_cnv)  # Assuming BND function is defined
  return(est_cnv)
}

smooth_block <- function(arr) {
  med <- median(arr)
  med_arr <- matrix(med, nrow = nrow(arr), ncol = ncol(arr))
  y <- (arr - med_arr) / 2 + med_arr
  return(y)
}

smooth_cnv <- function(bnd_list, dat) {
  last_idx <- 1
  result <- matrix(nrow = nrow(dat), ncol = 0)
  
  for (cidx in bnd_list) {
    x <- dat[, last_idx:cidx, drop = FALSE]
    y <- smooth_block(x)
    result <- cbind(result, y)
    last_idx <- cidx + 1
  }
  
  return(result)
}

# Define parameters
pdir <- paste0(input_dir, "/")
rdmatrix <- paste0(pdir, "RDmatrix.rds")
genome_cov <- paste0(pdir, "genome_cov.bed")

cat("Using input files:\n")
cat(rdmatrix, "\n")
cat(genome_cov, "\n")

if (!is.null(cell_node)) {
    cat(cell_node, "\n")
}

# Load data
result <- readRDS(rdmatrix)
df <- as.data.frame(result) 
df1 <- read.csv(genome_cov, sep = "\t")
colnames(df) <- colnames(df1)[4:ncol(df1)]

# Get reference bins
ref_bins <- apply(df1, 1, function(x) paste0(x["chr"], ":", x["start"], "-", x["end"]))
dc_ref <- setNames(seq_len(nrow(df1)), paste0(df1$chr, df1$end))

Y <- as.matrix(df)
Y_plus_1 <- Y + 1
column_means <- colMeans(Y_plus_1)

# Normalize each column of Y_plus_1 by dividing it by its respective column mean
Y_norm <- sweep(Y_plus_1, 2, column_means, FUN = "/")

# Calculate the column std of normal
cell_std <- apply(Y_norm, 2, sd)
std_thresh <- 0.25
cat(sprintf("Use standard deviation threshold %f to separate normal and cancer cells.\n", std_thresh))

# Create histogram
hist(cell_std, breaks = 50, main = "", xlab = "Standard Deviation")
ggsave(filename = paste0(output_dir, "/cell_std_", prefix, ".pdf"))

all_cells <- colnames(df)
n_cell_idx <- which(cell_std <= std_thresh)
c_cell_idx <- which(cell_std > std_thresh)

if (length(n_cell_idx) < 20) {
    idx_list <- get_norm_tm_cell(Y_norm, t = 50)
    n_cell_idx <- idx_list[[1]]
    c_cell_idx <- idx_list[[2]]
}

Y_norm_n_cell <- Y_norm[, n_cell_idx, drop = FALSE]
Y_norm_c_cell <- Y_norm[, c_cell_idx, drop = FALSE]
normal_cells <- all_cells[n_cell_idx]
cancer_cells <- all_cells[c_cell_idx]

cat(sprintf("Got %d normal cells, %d cancer cells\n", length(n_cell_idx), length(c_cell_idx)))
cat("Normal cells:\n")
print(normal_cells)
cat("Cancer cells:\n")
print(cancer_cells)

cat("Processing cancer cells...", "\n")
df4 <- as.data.frame(Y_norm_c_cell)
colnames(df4) <- cancer_cells

y4 <- t(Y_norm_c_cell)

# Extract fingerprint features
diff <- matrix(apply(y4, 2, max) - apply(y4, 2, min), ncol = 1)
temp <- t(diff)
thresh <- quantile(diff, 0.9)
cat("bin variance thresh:", thresh, "\n")
pos <- temp[1,] < thresh
keep <- temp[1,] >= thresh
cat(sum(keep), "out of", length(diff), "bins are extracted as fingerprint", "\n")

data <- y4[, keep, drop = FALSE]
data_norm <- sweep(data, 1, apply(data, 1, mean), "/")

# Hierarchical clustering for ploidy
cell_cl_dist <- 30
cat(paste("Use threshold cell distance threshold", cell_cl_dist, "to perform cell clustering"), "\n")
hc <- hclust(dist(data_norm), method = "ward.D2")
labels <- cutree(hc, h = cell_cl_dist)
cat("clustering result:\n", labels, "\n")
num_cell_clusters <- length(unique(labels))
cat("num cell clusters: ", num_cell_clusters, "\n")

y4_by_cluster <- list()
for (i in unique(labels)) {
  y4_by_cluster[[i]] <- y4[which(labels == i), , drop = FALSE]
}

labels <- c(labels, rep(length(unique(labels)) + 1, length(n_cell_idx)))
y4_by_cluster[[length(unique(labels))]] <- t(Y_norm_n_cell)
num_cell_clusters <- length(unique(labels))
cat("(adding normal cell) num cell clusters: ", num_cell_clusters, "\n")

# Get global ploidy
cat("Estimate cancer cell global ploidy\n") 
global_p <- get_ploidy_from_one_cell_cnv(colMeans(y4))

# CNV Estimation
CNV <- list()
new_cell_list <- list()
all_cells <- c(cancer_cells, normal_cells)
final_p_list <- rep(-1.1, length(labels))

for (j in unique(labels)) {
  n <- dim(y4_by_cluster[[j]])[1]
  
  new_cells <- all_cells[labels == j]
  new_cell_list <- c(new_cell_list, new_cells)
  cat("-------------Process cluster ", j + 1, paste0("(", n, " cells)"), "\n")
  
  # Get bin consensus
  Y_one_cl <- y4_by_cluster[[j]]
  Y_mean <- matrix(colMeans(Y_one_cl), nrow = 1)
  
  dat <- Y_mean[1, ]
  
  # Get changing rate
  f <- 6
  diffl <- diff(dat, lag = f)
  
  # Change rate cut-off
  if (crm == 'q') {
    t <- round(quantile(diffl, 0.84), 2)
  } else if (crm == 'h') {
    t <- 0.1
  } else {
    stop("change_rate_mode can only be 'q' or 'h'")
  }
  
  cat("changing rate cut off ", round(t, 2), "\n")
  
  # Estimate ploidy
  if (j == max(unique(labels))) {
    cat(paste0("Estimate normal cell cluster ", j, " ploidy"), "\n")
    p <- 2.0
    cat(p, "\n")
  } else if (n >= 5) {
    cat(paste0("Estimate cancer cell cluster ", j, " ploidy"), "\n")
    # In Python: p = get_ploidy_mt_cell(Y_one_cl)
    p <- get_ploidy_mt_cell(Y_mean)
    cat(p, "\n")
  } else {
    cat(paste0("small cancer cell cluster ", j, " (", n, " cells), use global ploidy ", global_p), "\n")
    p <- global_p
  }
  
  final_p_list[labels == j] <- p
}

cat("The ploidy list\n")
print(final_p_list)

cat("-----------------------------------------------------\n\n")
cat("                   Second round clustering           \n\n")
cat("-----------------------------------------------------\n\n")

# Hierarchical clustering
cell_cl_dist <- 14

hc <- hclust(dist(data_norm), method="ward.D2")
labels <- cutree(hc, h=cell_cl_dist)  

cat(sprintf("Use cell distance threshold %d to perform cell clustering\n", cell_cl_dist))
cat(sprintf("clustering result:\n"))
print(labels)

num_cell_clusters <- length(unique(labels))
cat(sprintf("num cell clusters: %d\n", num_cell_clusters))

y4_by_cluster <- vector("list", length(unique(labels)))

for (i in unique(labels)) {
  y4_by_cluster[[i]] <- y4[labels == i, , drop = FALSE]
}

labels <- c(labels, rep(max(labels) + 1, length(n_cell_idx)))
y4_by_cluster[[length(y4_by_cluster) + 1]] <- t(Y_norm_n_cell)

# Update the number of clusters
num_cell_clusters <- length(unique(labels))
cat(sprintf("(adding normal cell) num cell clusters: %d\n", num_cell_clusters))

number_of_bins <- ncol(Y_one_cl)
CNV <- matrix(nrow = 0, ncol = number_of_bins)
new_cell_list <- list()
all_cells <- c(cancer_cells, normal_cells)

for (j in seq_along(unique(labels))) {
  
  Y_one_cl <- y4_by_cluster[[j]]
  n <- nrow(Y_one_cl)
  
  new_cells <- all_cells[labels == unique(labels)[j]]
  new_cell_list <- c(new_cell_list, new_cells)
  cat(sprintf("-------------Process cluster %d (%d cells)\n", j, n))
  
  # ----------------- get bin consensus
  Y_mean <- matrix(colMeans(Y_one_cl), nrow = 1)
  dat <- Y_mean[1, ]
  
  # ------------- get changing rate
  f <- 6
  diffl <- nbdiff(dat, f)    
  #--------change rate cut-off
  if (crm == 'q') {
    t <- round(quantile(diffl, probs = 0.84), 2)
  } else if (crm == 'h') {
    t <- 0.1
  } else {
    stop("change_rate_mode can only be 'q' or 'h'")
  }
  
  cat(sprintf("changing rate cut off %.2f\n", t))
  
  #---------------get segmentation
  apx_ids <- which(diffl >= t)
  cl_list <- list(apx_ids[1])
  
  for (i in 2:length(apx_ids)) {
    cid <- apx_ids[i]
    lid <- tail(cl_list[[length(cl_list)]], 1)
    if (cid - lid <= 1) {
      cl_list[[length(cl_list)]] <- c(cl_list[[length(cl_list)]], cid)
    } else {
      cl_list <- c(cl_list, list(cid))
    }
  }
  clrg_list <- lapply(cl_list, function(cl) {
    if (length(cl) >= 4) {
      return(c(min(cl), max(cl)))
    }
  })
  clrg_list <- Filter(Negate(is.null), clrg_list) 
  # print(clrg_list)
  #---------transform bnd
  a2 <- sapply(clrg_list, function(cl) round(mean(cl) - 0.5))
  bnd_list <- c(as.integer(a2), length(dat))
  
  # ------------------ estimate CNV
  idx_list <- which(labels == unique(labels)[j]) 
  Y_one_cl_sm <- smooth_cnv(bnd_list, Y_one_cl)   
  for (m in 1:nrow(Y_one_cl)) {
    p_use <- final_p_list[idx_list[m]]
    cat("p_use:", p_use, "\n")
    est_cnv <- get_by_cell_cnv(bnd_list, Y_one_cl_sm[m, ], p_use)   
    CNV <- rbind(CNV, est_cnv)
  }
}
nbin <- length(CNV[[1]])

if (!is.null(cell_node) && nzchar(cell_node)) {
  df3 <- read.table(cell_node, sep = '\t', header = FALSE)
  dc_cn <- setNames(as.character(df3[,2]), df3[,1])
  cells <- ifelse(new_cell_list %in% names(dc_cn), paste(dc_cn[new_cell_list], new_cell_list, sep="_"), new_cell_list)
} else {
  cells <- new_cell_list
}

# Write to CSV 
bins <- apply(df1, 1, function(x) paste(x["chr"], paste0(x["start"], "-", x["end"]), sep=":"))
dfcnv <- data.frame(CNV)
colnames(dfcnv) <- bins
dfcnv$cell <- cells
dfcnv <- data.frame(dfcnv, row.names = 'cell')
output_path <- paste0(output_dir, "/CNV_", prefix, ".csv")
write.csv(dfcnv, file = output_path, row.names = TRUE)

# dfcnv1 <- as.data.frame(t(CNV_matrix))
# colnames(dfcnv1) <- cells
# dfcnv1$bin <- bins
# dfcnv1 <- data.frame(dfcnv1, row.names = 'bin')
# write.csv(dfcnv1, file = paste0(output_path, ".T"), row.names = TRUE)

