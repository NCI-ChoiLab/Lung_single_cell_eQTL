# Merge cis nominal QTL results from TensorQTL

# Thong Luong
# Jan 3rd 2025

library(dplyr)
library(ggplot2)
library(future)
library(arrow)
library(Matrix)

# read list of cell types 
celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

chr_num <- c(1:22,"X")
chrs <- paste0("chr",chr_num)
tensor_path <- "/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/"
dir.create(file.path(tensor_path,"output_rsid"), showWarnings = FALSE)
full_sig_list <- readRDS("/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/full_sig_list.rds") # a list of results from each cell types only including eGenes
for (celltype in celltypes) {
  rst.list <- list()
  for (chr in chrs) {
    rst <- read_parquet(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/', celltype, '/', celltype, '.cis_qtl_pairs.', chr,'.parquet'))
    rst.list[[chr]] <- rst
  }
  qtl_rst_all <- Reduce(rbind, rst.list)
  qtl_rst_all$celltype <- celltype
  dir.create(file.path(tensor_path,"output_rsid",celltype), showWarnings = FALSE)
  output_path <- file.path(tensor_path,"output_rsid",celltype)
  write.table(qtl_rst_all, file = paste0(output_path,"/qtl_results_all.txt"), sep = "\t",
              quote = FALSE, row.names = FALSE)
  sig_eGenes <- full_sig_list[[celltype]]
  rownames(sig_eGenes) <- sig_eGenes$phenotype_id
  tmp <- qtl_rst_all %>% filter(phenotype_id %in% sig_eGenes$phenotype_id) %>% 
    group_by(phenotype_id) %>% 
    mutate(significance = pval_nominal < sig_eGenes[phenotype_id,"pval_nominal_threshold"])
  
  sig_qtl <- tmp[which(tmp$significance),]
  write.table(sig_qtl, file = paste0(output_path,"/top_qtl_results_all_FDR0.05.txt"), sep = "\t",
              quote = FALSE, row.names = FALSE)
}
