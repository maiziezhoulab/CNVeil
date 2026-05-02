#!/usr/bin/env Rscript
# .libPaths(c("/home/luoc1/R/rlib-4.0.5/","/data/maiziezhou_lab/Softwares/R_4.0_package_path/","/home/yuanw2/R/AneuFinder")) 

library("optparse")
# library(AneuFinder)

read_tab_file <- function(file_path) {
  # Read in the file as a data frame
  df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # Convert the data frame to a named list
  my_list <- as.list(df[[2]])
  names(my_list) <- df[[1]]
#   print(length(my_list))

  return(my_list)
}

filter_by_named_list <- function(named_list, df) {
  # Extract the row indices from the row names of the data frame
  row_indices <-  row.names(df)
#   print(row_indices)
  
  # Split the row indices by '-' and extract the first element
  first_elements <- sapply(strsplit(row_indices, ".bam"), function(x) x[1])
#   print(first_elements)
#   print(named_list)
#   print(first_elements)
  # Filter the data frame by the names in the named list
  filtered_df <- df[first_elements %in% names(named_list), ]
#   print(dim(filtered_df))
  return(filtered_df)
}

modify_row_names <- function(df, named_list) {
  # Extract the row indices from the row names of the data frame
  row_indices <-  row.names(df)
  
  # Split the row indices by '-' and extract the first element
  first_elements <- sapply(strsplit(row_indices, ".bam"), function(x) x[1])
#   print(row_indices)
#   print(first_elements)
  # Replace the row names of the data frame with the modified row names
  new_row_names <- paste0(named_list[first_elements], "_", first_elements)
  row.names(df) <- new_row_names
  
  return(df)
}



combine_coords <- function(chr_list, start_list, end_list) {
  if (length(chr_list) != length(start_list) || length(start_list) != length(end_list)) {
    stop("Error: input lists must have the same length.")
  }
  
  # Create a character vector to store the combined coordinates
  coords <- character(length(chr_list))
  
  # Loop over the input lists and combine the coordinates
  for (i in seq_along(chr_list)) {
    coords[i] <- paste0(chr_list[i], ":", start_list[i], "-", end_list[i])
  }
  
  return(coords)
}


# transpose_with_ref_bin <- function(df) {
#   # Transpose the data frame using the "ref_bin" column as column names
#   df_t <- tidyr::pivot_wider(df, names_from = ref_bin_list, values_from = names(df)[1:dim(df)[2]-1])
  
#   # Return the transposed data frame
#   return(df_t)
# }

loadFromFiles <- function(files, check.class=c('GRanges', 'GRangesList', 'aneuHMM', 'aneuBiHMM')) {

    # ptm <- startTimedMessage("Loading data from files ...")
    if (is.null(files)) {
        # stopTimedMessage(ptm)
        return(files)
    }
    if (any(! check.class %in% c('GRanges', 'GRangesList', "aneuHMM", "aneuBiHMM"))) {
        stop("Argument 'check.class' must contain any combination of c('", paste0(c('GRanges', 'GRangesList', "aneuHMM", "aneuBiHMM"), collapse="', '"), "').")
    }
    .check_class <- function(x)
        any(vapply(check.class, is, logical(1), object=x))
    modellist <- list()
    if (is.character(files)) {
        for (file in files) {
            # print(file)
            temp.env <- new.env()
            model <- get(load(file, envir=temp.env), envir=temp.env)
            # if (!.check_class(model)) {
            #     stop("File '", file, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            # }
            modellist[[file]] <- model
        }
    } else if (.check_class(files)) {
        modellist[[1]] <- files
    } else if (is.list(files)) {
        for (file in files) {
            model <- file
            if (!.check_class(model)) {
                stop("List entry '", length(modellist)+1, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[length(modellist)+1]] <- model
        }
        names(modellist) <- names(files)
    } else if (!.check_class(files)) {
        stop("Input does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
    }
    # stopTimedMessage(ptm)
    return(modellist)

}

load_afdata <- function(input_dir) {
    AneuFinder_results_path = input_dir
    # print(AneuFinder_results_path)
    AneuFinder_results_file = list.files(AneuFinder_results_path, full.names=TRUE)
    # Exclude paths containing "rg_cb"
    # AneuFinder_results_file <- AneuFinder_results_file[!grepl("rg_cb", AneuFinder_results_file)]
    AneuFinder_results <- loadFromFiles(AneuFinder_results_file, check.class=c('GRanges', 'GRangesList'))
    print("================haha")
    AneuFinder.CNV = c()
    cell_names = c()
    # print(AneuFinder_results)
    # print(length(AneuFinder_results))
    for(i_file in 1:length(AneuFinder_results)){
        # print(i)
    i_result = data.frame(AneuFinder_results[[i_file]]$bins)
    granges_data <- AneuFinder_results[[i_file]]@unlistData
    metadata <- granges_data@elementMetadata@listData$counts



    # print(metadata)

    # str(AneuFinder_results[[i_file]])
    # print(colnames(i_result))
    i_AneuFinder.CNV = matrix(metadata, ncol = 1)
    i_cellnames = unlist(strsplit(basename(AneuFinder_results_file[i_file]), split = "_"))[[1]]

    i_cellnames <- sub("\\.bam$", "", i_cellnames)

    if (!i_cellnames %in% cell_names) {
        cell_names <- c(cell_names, i_cellnames)
        AneuFinder.CNV <- cbind(AneuFinder.CNV, i_AneuFinder.CNV)
    }


    }

    # Example: Assume `compressed_grl` is your CompressedGRangesList object
    # Extract the unlisted GRanges object
    granges <- granges_data

    # Extract components
    chromosomes <- as.character(seqnames(granges)) # Extract chromosome names
    start_positions <- start(granges)             # Extract start positions
    end_positions <- end(granges)                 # Extract end positions

    # Combine into bins
    bins <- paste(chromosomes, paste(start_positions, end_positions, sep = "-"), sep = ":")

    # View the resulting bins
    # print(bins)


    # AneuFinder.ref = i_result[,1:3]
    # names(AneuFinder.ref) = c("chrom", "start", "end")

    # chrom=AneuFinder.ref$chrom
    # start=AneuFinder.ref$start
    # end=AneuFinder.ref$end

    # ref_bin_list = combine_coords(chrom, start, end)


    # CNV = t(AneuFinder.CNV)
    # colnames(CNV) =  ref_bin_list
    rownames(AneuFinder.CNV) = bins
    colnames(AneuFinder.CNV) = cell_names
    # rownames(CNV) = cell_names
    df = as.data.frame(AneuFinder.CNV)

    return(df)
}




option_list = list(
    make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="aneufinder working dir", metavar="character"),
    make_option(c("-t", "--cell_node_file"), type="character", default="", 
              help="cell to node file[default= %default]", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default="", 
              help="out file [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_dir = opt$input_dir
cell_node_file = opt$cell_node_file
# Assuming 'opt' is a dictionary or an object with a property 'cell_node_file'
# cell_node_file = opt.get('cell_node_file', None)  # Use get() to safely retrieve the file path

outfile = opt$outfile 

AneuFinder_results_path = input_dir

if (!dir.exists(AneuFinder_results_path)) {
    message(paste0(AneuFinder_results_path, " not exists"))
    q()
    stop()
}


df = load_afdata(input_dir)

# Only run the code if cell_node_file is provided
if (!is.null(cell_node_file) && cell_node_file != "") {
    # Read the cell node file
    cell_node_list = read_tab_file(cell_node_file)
    
    # Filter the dataframe by the cell node list
    df = filter_by_named_list(cell_node_list, df)
    
    # Modify row names in the dataframe
    df = modify_row_names(df, cell_node_list)
} else {
  cat("No cell_node_file provided. Skipping the processing.\n")
}

# Write data frame to CSV with index
write.csv(df, file = outfile , row.names = TRUE)


