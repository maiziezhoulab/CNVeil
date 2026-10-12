# Combine chromosome matrices using the original filenames and SNP filter.
Cbn_matrix <- function(dir_path = "./", vcf_path = "./vcf_sub/",
                       barcodes_path = "./barcodes.tsv", samplename = "Sample",
                       plot_stat = TRUE) {
    barcode <- read.table(barcodes_path, sep = "\t", stringsAsFactors = FALSE)
    alt_all <- ref_all <- var_all <- NULL
    for (ii in 1:22) {
        alt <- Matrix::readMM(file.path(dir_path, "chr", paste0("chr", ii, "_matrix_alt")))
        ref <- Matrix::readMM(file.path(dir_path, "chr", paste0("chr", ii, "_matrix_ref")))
        var <- read.table(file.path(vcf_path, paste0(samplename, "_phased_chr", ii,
                                                   "_chrprefix.vcf")),
                          header = FALSE, stringsAsFactors = FALSE)
        ind0 <- Matrix::rowSums(alt) > 0
        var_all <- rbind(var_all, var[ind0, , drop = FALSE])
        alt_all <- rbind(alt_all, alt[ind0, , drop = FALSE])
        ref_all <- rbind(ref_all, ref[ind0, , drop = FALSE])
        cat(paste0("chr", ii, " "))
    }
    colnames(alt_all) <- barcode[, 1]
    colnames(ref_all) <- barcode[, 1]
    saveRDS(alt_all, file.path(dir_path, "alt_all.rds"))
    saveRDS(ref_all, file.path(dir_path, "ref_all.rds"))
    saveRDS(var_all, file.path(dir_path, "var_all.rds"))
    total_all <- alt_all + ref_all
    af <- Matrix::rowSums(alt_all) / Matrix::rowSums(total_all)
    af[is.na(af)] <- 0
    if (plot_stat) {
        pdf(file.path(dir_path, "statistics__raw.pdf"))
        par(mfrow = c(2, 1))
        hist(Matrix::colSums(total_all), main = paste(samplename, "cell coverage"), breaks = 300)
        hist(Matrix::rowSums(total_all), main = paste(samplename, "SNP coverage"),
             xlim = c(0, 100), breaks = 1000000)
        hist(af, 100, main = paste(samplename, "VAF for each SNP"))
        dev.off()
    }
}
