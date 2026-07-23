.libPaths("/home/yuanw2/R/scDNA-seq")
library(ASCAT)
library(dplyr)

# args <- commandArgs(trailingOnly = TRUE)

sample_id = args[1]#"TN1"
normal_srr_id = args[2]
tumor_srr_id = args[3]
Tumour.bam = args[4]#path to bam
Normal.bam = args[5]
Tumour_name = paste0(sample_id, "Tumour")
Normal_name = paste0(sample_id, "Normal")
Gender = args[6]# e.g. 'XX'
Worksheet = args[7]#see details on ASCAT github

ascat.prepareTargetedSeq(
  '/data/maiziezhou_lab/Weiman/single_cell_project/Validation/exome_analysis/ASCAT/Worksheet_1.tsv', 
  '/data/maiziezhou_lab/Weiman/single_cell_project/Validation/exome_analysis/', 
  '/data/maiziezhou_lab/Weiman/reference/1000Genome/WES/hg19/G1000_allelesAll_hg19/G1000_alleles_hg19_chr', 
  '/data/maiziezhou_lab/Weiman/single_cell_project/Validation/exome_analysis/SeqCap_bed_ref/modified_overlapping_snps.bed', 
  '/data/maiziezhou_lab/Weiman/software/alleleCount/bin/alleleCounter', 'hg19', nthreads=4,
    minCounts=20, is_chr_based=FALSE, chrom_names=c(1:22, "X"), 
    min_base_qual=20, min_map_qual=35, ref.fasta=NA, plotQC=TRUE)

info = fread("./TN_sample_info.tsv", header = TRUE)
for (i in c(1:8)){
   ascat.prepareHTS(
  tumourseqfile = info$Tumor_file[i],
  normalseqfile = info$Normal_file[i],
  tumourname = info$Tumor_ID[i],
  normalname = info$Normal_ID[i],
  allelecounter_exe = "/data/maiziezhou_lab/Weiman/software/alleleCount/bin/alleleCounter",
  alleles.prefix = '/data/maiziezhou_lab/Weiman/reference/1000Genome/WES/hg19/G1000_allelesAll_hg19/G1000_alleles_hg19_chr',
  loci.prefix = '/data/maiziezhou_lab/Weiman/reference/1000Genome/WES/hg19/G1000_lociAll_hg19/G1000_loci_hg19_chr',
  gender = 'XX',
  genomeVersion = "hg19",
  nthreads = 8,
  tumourLogR_file = paste0("TN", i, "Tumor_LogR.txt"),
  tumourBAF_file = paste0("TN", i, "Tumor_BAF.txt"),  
  normalLogR_file = paste0("TN", i, "Germline_LogR.txt"), 
  normalBAF_file = paste0("TN", i, "Germline_BAF.txt")) 
  
ascat.bc = ascat.loadData(Tumor_LogR_file = paste0("TN", i, "Tumor_LogR.txt"), Tumor_BAF_file = paste0("TN", i, "Tumor_BAF.txt"), Germline_LogR_file = paste0("TN", i, "Germline_LogR.txt"), Germline_BAF_file = paste0("TN", i, "Germline_BAF.txt"), gender = 'XX', genomeVersion = "hg19", isTargetedSeq=T)
ascat.plotRawData(ascat.bc, img.prefix = paste0("Before_correction_", "TN", i))
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "./GC_G1000_hg19.txt", replictimingfile = "./RT_G1000_hg19.txt") ## download?
ascat.plotRawData(ascat.bc, img.prefix = paste0("After_correction_", "TN", i))
ascat.bc = ascat.aspcf(ascat.bc, penalty=25)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = paste0("TN", i, '_ASCAT_objects.Rdata'))
}
