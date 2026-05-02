.libPaths(c("/home/luoc1/R/rlib-4.0.5/","/home/yuanw2/R/rlib-4.0.5","/data/maiziezhou_lab/Softwares/R_4.0_package_path/")) 
#.libPaths("~/R/rlib-4.0.5")
# Load the argparse package
library(argparse)

# Create an ArgumentParser object
parser <- argparse::ArgumentParser(description = "Description of your script")

# Add input arguments
parser$add_argument("--input_dir", '-i', dest = "input_dir", help = "Path to the input file")
parser$add_argument("--prefix", '-px', dest = "prefix", help = "prefix,default = sample", default = "sample")

# Parse the command-line arguments
args <- parser$parse_args()

input_dir <- args$input_dir
prefix <- args$prefix
# RDmatrix<-args$RDmatrix
# output_file <- args$output_file
# option1 <- args$option1
# importFrom DescTools AUC


library(RColorBrewer)
# library(SCOPE)
library(stringr)






# plotiCN <- function(iCNmat, ref, Gini, outfile, annotation = NULL,
#                      plot.dendrogram = TRUE, show.names = FALSE) {
plotiCN <- function(iCNmat, ref,  outfile, annotation = NULL,
                     plot.dendrogram = TRUE, show.names = FALSE) {
  smart_image <- function(mat, ...) {
    image(t(mat[rev(seq(nrow(mat))), ]),useRaster = TRUE, ...)
  }
  hm_col <- c("#2166AC", "#92C5DE", "#FDFDFD", "#FDDBC7",
              "#F4A582", "#D6604D", "#B2182B", "#67001F")
  
  if (!is.matrix(iCNmat)) {
    stop("Invalid plot object: must be an integer matrix. \n")
  }
  if (length(ref) != nrow(iCNmat)) {
    stop("Invalid GRanges object: length of ref and # of
            rows in iCNmat must be the same")
  }
  if (!is.null(annotation)) {
    if (!is.null(dim(annotation))) {
      stop("Invalid annotation object: has to be a vector or
                factor with the same # of cells as that of iCNmat")
    }
    if (length(annotation) != ncol(iCNmat)) {
      stop("Invalid annotation object: length of annotation and
                # of cells in iCNmat must be the same")
    }
  }
  if (show.names){
    if (is.null(colnames(iCNmat))){
      stop("Invalid plot object: cell names cannot be NULL")
    }
  }
  # if (length(Gini) != ncol(iCNmat)) {
  #   stop("Invalid Gini object: length of Gini coefficient and
  #           # of cells in iCNmat must be the same")
  # }
  
  # page setup
  if (is.null(annotation)) {
    if (plot.dendrogram) {
      if (show.names) {
        mm <- matrix(c(0, 0, 4, 0, 0,
                       2, 3, 1, 7, 5,
                       2, 3, 1, 7, 0,
                       2, 3, 1, 7, 6), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 5, 0.2, 0.5)
        mw <- mw/sum(mw) 
      } else {
        mm <- matrix(c(0, 0, 4, 0,
                       2, 3, 1, 5,
                       2, 3, 1, 0,
                       2, 3, 1, 6), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 5, 0.5)
        mw <- mw/sum(mw) 
      }
    } else {
      if (show.names) {
        mm <- matrix(c(0, 0, 3, 0, 0,
                       0, 2, 1, 6, 4,
                       0, 2, 1, 6, 0,
                       0, 2, 1, 6, 5), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 5, 0.2, 0.5)
        mw <- mw/sum(mw)
      } else {
        # print("here I am")
        # mm <- matrix(c(0, 0, 3, 0,
        #                0, 2, 1, 4,
        #                0, 2, 1, 0,
        #                0, 2, 1, 5), nrow = 4, byrow = TRUE)




        # mh <- c(2, 20, 20, 20)
        # mh <- mh/sum(mh)
        # mw <- c(0.25, 0.1, 5, 0.5)
        # mw <- mw/sum(mw)

        mm <- matrix(c(2, 
                        1, 
                        1, 
                        1 ), nrow = 4, byrow = TRUE)

        mh <- c(2,20,20,20)
        mh <- mh/sum(mh)
        mw <- c(1)
      }
    }
  } else {
    if (plot.dendrogram) {
      if (show.names) {
        mm <- matrix(c(0, 0, 0, 5, 0, 0,
                       2, 3, 4, 1, 9, 6,
                       2, 3, 4, 1, 9, 7,
                       2, 3, 4, 1, 9, 8), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 0.1, 5, 0.2, 0.5)
        mw <- mw/sum(mw)
      } else {
        mm <- matrix(c(0, 0, 0, 5, 0,
                       2, 3, 4, 1, 6,
                       2, 3, 4, 1, 7,
                       2, 3, 4, 1, 8), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 0.1, 5, 0.5)
        mw <- mw/sum(mw)
      }
    } else {
      if (show.names) {
        mm <- matrix(c(0, 0, 0, 4, 0, 0,
                       0, 2, 3, 1, 8, 5,
                       0, 2, 3, 1, 8, 6,
                       0, 2, 3, 1, 8, 7), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 0.1, 5, 0.2, 0.5)
        mw <- mw/sum(mw)
      } else {
        mm <- matrix(c(0, 0, 0, 4, 0,
                       0, 2, 3, 1, 5,
                       0, 2, 3, 1, 6,
                       0, 2, 3, 1, 7), nrow = 4, byrow = TRUE)
        mh <- c(2, 20, 20, 20)
        mh <- mh/sum(mh)
        mw <- c(0.25, 0.1, 0.1, 5, 0.5)
       
        mw <- mw/sum(mw)
      }
    }
  }
  png(outfile, width = 2500,
      height = 1500, pointsize = 25, type="cairo")
  
  # pdf(outfile, width = 25,
  #     height = 11,  )

  # print(mm)
  # print(mw)
  # print(mh)
  layout(mm, widths = mw, heights = mh)
  par(mar = rep(0, 4))
  
  iCNmat <- round(iCNmat)
  if (!is.null(annotation)) {
    annotation <- as.factor(annotation)
  }
  
  chr.pos <- rep(NA, length(unique(ref)))
  for (chri in seq_len(length(chr.pos))) {
    chr.pos[chri] <- length(ref[which(as.character(
      ref) == as.character(unique(
        ref))[chri])])
    
  }
  chr.pos <- cumsum(chr.pos)
  xpos <- round(c(0, chr.pos[seq_len(length(chr.pos)-1)]) +
                  (chr.pos - c(0, chr.pos[seq_len(length(chr.pos)-1)]))/2)
  
  # 1) iCN heatmap
  dat <- t(iCNmat)
  #dat <- t(CNV_total)
  dat[dat >= 7] <- 7
  dat[dat <= 0] <- 0
  iCNtab <- as.numeric(names(table(dat)))
  rclust <- hclust(dist(dat))
  # rclust = readRDS("/data/maiziezhou_lab/Weiman/Simulator/Tool/Weiman/dbscan+secnv/ploidy5/rclust")
  # r_order = readRDS("/data/maiziezhou_lab/Weiman/Newdataset/SCOPE/draft/getfgc_group/r_order.rds")
  dat <- dat[rclust$order, ]
  # print(dat[1,])
  
  smart_image(dat, col = hm_col[iCNtab + 1], xaxs = "i",
              yaxs = "i", axes = FALSE)
  abline(v = 0, lwd = 2)
  for (i in seq_len(length(chr.pos))) {
    abline(v = chr.pos[i]/length(ref), lwd = 2)
  }
  
  # 2) hclust
  if (plot.dendrogram) {
    plot(rev(as.dendrogram(rclust)), leaflab = "none",
         horiz = TRUE, axes = FALSE, yaxs = "i")
  }
  
  # 3) Gini annotation
  # anno.Gini <- matrix(Gini[rclust$order], nrow = nrow(dat),
  #                     ncol = 1)
  # col.Gini <- gplots::colorpanel(50, "#F7FBFF", "#084594")
  # smart_image(anno.Gini, col = col.Gini, xaxs = "i",
  #             yaxs = "i", axes = FALSE)
  
  if (!is.null(annotation)) {
    # 4) Customized annotation
    anno.level <- levels(annotation)
    anno.mat <- matrix(match(annotation[rclust$order], anno.level),
                       nrow = nrow(dat), ncol = 1)
    col.anno <- brewer.pal(n = 12, name = "Set3")[
      sort(unique(match(annotation[rclust$order], anno.level)))]
    smart_image(anno.mat, col = col.anno, xaxs = "i",
                yaxs = "i", axes = FALSE)
  }
  
  # 5) chromosome
  anno.chrom <- NULL
  for (i in seq_len(length(chr.pos))) {
    if (i%%2 == 1) {
      temp <- matrix(rep(1, length(which(as.character(
        ref) == as.character(unique(
          ref))[i]))), nrow = 1)
      anno.chrom <- cbind(anno.chrom, temp)
    } else {
      temp <- matrix(rep(2, length(which(as.character(
        ref) == as.character(unique(
          ref))[i]))), nrow = 1)
      anno.chrom <- cbind(anno.chrom, temp)
    }
  }

  # image(t(anno.chrom), col = c("white", "white"),
  #       xaxs = "i", yaxs = "i", axes = FALSE)
  # image(t(anno.chrom), col = c("gray", "black"),
  #       xaxs = "i", yaxs = "i", axes = FALSE)
  # pos.text <- xpos/length(ref)
  # chr.noprint <- as.character(unique(ref))
  # chr.print <- substr(chr.noprint, 4, nchar(chr.noprint))

  # text(pos.text, 0.2, chr.print, col = c("black", "grey"), 
  #      cex = 1.5)
  
  # 6) Gini legend
  # par(mar = c(2, 2, 2, 4))
  # image(1, seq_len(length(brewer.pal(n = 8, name = "Blues"))),
  #       t(as.matrix(seq_len(length(brewer.pal(n = 8, name = "Blues"))))),
  #       col = brewer.pal(n = 8, name = "Blues"), xlab = "", ylab = "",
  #       xaxt = "n", yaxt = "n", bty = "n")
  # axis(4, at = c(1, length(brewer.pal(n = 8, name = "Blues"))),
  #      labels = round(c(min(Gini), max(Gini)), 2), col.ticks = "white",
  #      col = NA, lwd.ticks = 0, cex.axis = 1.5, las = 2, font = 2)
  # title("Gini", cex.main = 1.5)
  
  if (!is.null(annotation)) {
    plot(0, 0, type = "n", axes = FALSE)
    legend("center", legend = sort(unique(annotation)),
           col = col.anno, pch = 15, bty = "n", cex = 1.5)
  }
  
  # 8) iCN legend
  # par(mar = c(2, 2, 2, 4))
  # image(1, seq_len(length(hm_col[iCNtab + 1])),
  #       t(as.matrix(seq_len(length(hm_col[iCNtab + 1])))),
  #       col = hm_col[iCNtab + 1], xlab = "", ylab = "",
  #       xaxt = "n", yaxt = "n", bty = "n")
  # axis(4, at = seq_len(length(hm_col[iCNtab + 1])),
  #      labels = c(0:7)[iCNtab + 1], col.ticks = "white", col = NA,
  #      lwd.ticks = 0, cex.axis = 1.5, las = 2, font = 2)
  # title("integer CN", cex.main = 1.5)
  
  # 9) cell names
  if (show.names) {
    par(mar = rep(0, 4))
    image(1, seq_len(ncol(iCNmat)), 
          t(as.matrix(seq_len(ncol(iCNmat)))),
          col = "white", xlab = "", ylab = "",
          xaxt = "n", yaxt = "n", bty = "n")
    text(1, rev(seq_len(ncol(iCNmat))), 
         labels = colnames(iCNmat)[rclust$order], las = 2, cex = 1.5, font = 2)
  }


  # image(t(anno.chrom), col = c("white", "white"),
  #       xaxs = "i", yaxs = "i", axes = FALSE)
  image(t(anno.chrom), col = c("gray", "black"),
        xaxs = "i", yaxs = "i", axes = FALSE,useRaster = TRUE)
  pos.text <- xpos/length(ref)
  chr.noprint <- as.character(unique(ref))
  chr.print <- substr(chr.noprint, 4, nchar(chr.noprint))
  # print(pos.text)
  # print(chr.print)
  text(pos.text, 0.2, chr.print, col = c("black", "grey"), 
       cex = 1.5)

  dev.off()
}










# Y = readRDS("/data/maiziezhou_lab/Weiman/Old/TNBC_cell_cluster/SeCNV_correction/ploidy5/RDmatrix.rds")
# Y = readRDS(RDmatrix)
# print("-------------------")
# print(input_dir)
RDmatrix<- paste0(input_dir,'/',"RDmatrix_qc.csv")
input_file <- paste0(input_dir,'/',sprintf('CNV_%s.csv.T', prefix))
# input_file <- paste0(input_dir,'/',sprintf('CNV_%s_ordered_by_gold.csv', prefix))
filename <- sprintf("heatmap_%s.png", prefix)
output_file <- paste0(input_dir,'/',filename)
Y <- as.matrix(read.csv(RDmatrix ))



#---------- test only
# Y <- Y[, -ncol(Y)]  # Removes the last column
# CNV <- readRDS("CNV.rds")
my_df <- read.csv(input_file, row.names = 1)
CNV<-as.matrix(my_df)

# Check if the number of columns match
if (ncol(CNV) != ncol(Y)) {
  # Remove the last column of Y if column numbers do not match
  Y <- Y[, -ncol(Y)]
}

# CNV <- readRDS(input_file)
# Gini <- get_gini(CNV)
# print(length(Gini))
bin_list <- rownames(CNV)
ref <- str_extract(bin_list, "chr\\d+")

# # Initialize an empty vector to store the annotation values
# annotation <- c()

# # Create a for loop to generate annotation values
# for (i in 1:cell_clusters) {
#   annotation <- c(annotation, rep(i, sum(clusters_result_cell == i)))
# }

# print(dim(CNV))
# print(length(annotation))

# plot.test<- plotiCN(CNV, ref, Gini, output_file,annotation = NULL,
                  # plot.dendrogram = FALSE, show.names = FALSE)
plot.test<- plotiCN(CNV, ref,  output_file,annotation = NULL,
                  plot.dendrogram = FALSE, show.names = FALSE)
