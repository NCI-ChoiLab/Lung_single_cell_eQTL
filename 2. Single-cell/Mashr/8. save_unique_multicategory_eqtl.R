#Author: Thong Luong
#Date: May 29th 2025

library(QTLExperiment)
library(dplyr)
library(lattice)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(stringr)
library(tidyverse)

#saving mashr defined unique and multicategory eQTL

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Overlap_eQTL_cCRE/Shared_Unique/Mashr/')

tensor_sig_qtl = data.frame()
top_sig = data.frame()

celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

epithelial = c("Alv_trans",'AT1','AT2','Basal','Club','Goblet','Multiciliated','Sec_trans')
immune = c("Alv_mph","Bcells","CD4","CD8","Cla_mono",'DC1',"DC2","Int_mph_peri","Mast","Mig_DC",
           "Mono_mph","NK","Noncla_mono")
endothelial = c("EC_aero_cap","EC_art","EC_gen_cap","EC_ven_pul","EC_ven_sys","Lym_EC_mat","Lym_EC_pro")
stromal = c("Adv_fib","Alv_fib","Peri_fib","SM","Sub_fib")

for (c in celltypes){
  ts = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  tensor_sig_qtl = rbind(tensor_sig_qtl, ts)
  
  ts %>% dplyr::group_by(phenotype_id) %>% dplyr::slice(which.min(pval_nominal)) -> top_eqtl
  top_eqtl = as.data.frame(top_eqtl)
  top_sig = rbind(top_sig, top_eqtl)
}

tensor_sig_qtl$qtl = paste(tensor_sig_qtl$phenotype_id,tensor_sig_qtl$variant_id,sep='|')
sig_qtl_list = unique(tensor_sig_qtl$qtl)

top_sig$qtl = paste(top_sig$phenotype_id,top_sig$variant_id,sep='|')
top_sig_list = unique(top_sig$qtl)

mash = readRDS('/data/Choi_lung/TTL/tensor/mashr/h5_output_Sum_Final_chr_pos/FDR/output/mashr_10K_qtl.rds')
lfsr = as.data.frame(lfsrs(mash))
rm(mash)
lfsr = subset(lfsr, rownames(lfsr) %in% top_sig_list)
lfsr$eqtl = rownames(lfsr)

get_eqtl_list = function(l){
  lfsr_list = list()
  for (c in l){
    cell = lfsr[,c(c,'eqtl')]
    colnames(cell)[1] = 'lfsr'
    cell = subset(cell, lfsr < 0.05)
    lfsr_list[[c]] = rownames(cell)
  }
  return(lfsr_list)
}

lfsr_list = get_eqtl_list(celltypes)

unique_list = list()
for (c in celltypes){
  rest_of_group = lfsr_list
  rest_of_group[[c]] = NULL
  cell = lfsr_list[[c]]
  unique_list[[c]] = setdiff(unlist(cell),unique(unlist(rest_of_group)))
}

unique = data.frame()
for (c in celltypes){
  cell = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  cell$qtl = paste(cell$phenotype_id,cell$variant_id,sep='|')
  cell = subset(cell, qtl %in% unlist(unique_list[[c]]))
  unique = rbind(unique, cell)
}


epi_list = get_eqtl_list(epithelial)
endo_list = get_eqtl_list(endothelial)
imm_list = get_eqtl_list(immune)
strom_list = get_eqtl_list(stromal)

Cell_group = list()
Cell_group[['Epithelial']] = unique(unlist(epi_list))
Cell_group[['Endothelial']] = unique(unlist(endo_list))
Cell_group[['Immune']] = unique(unlist(imm_list))
Cell_group[['Stromal']] = unique(unlist(strom_list))

uniq <- sort(unique(unlist(Cell_group)))
Multi_category = list()
Multi_category[['2 or more']] = uniq[rowSums(sapply(Cell_group, `%in%`, x = uniq)) > 1]
Multi_category[['3 or more']] = uniq[rowSums(sapply(Cell_group, `%in%`, x = uniq)) > 2]
Multi_category[['All']] = uniq[rowSums(sapply(Cell_group, `%in%`, x = uniq)) > 3]

multi2 = data.frame()
multi3 = data.frame()
multi4 = data.frame()
for (c in celltypes){
  cell = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  cell$qtl = paste(cell$phenotype_id,cell$variant_id,sep='|')
  cell2 = subset(cell, qtl %in% unlist(Multi_category[['2 or more']]))
  cell3 = subset(cell, qtl %in% unlist(Multi_category[['3 or more']]))
  cell4 = subset(cell, qtl %in% unlist(Multi_category[['All']]))
  multi2 = rbind(multi2, cell2)
  multi3 = rbind(multi3, cell3)
  multi4 = rbind(multi4, cell4)
}

saveRDS(unique,'unique.rds')
#saveRDS(multi2,'multi2.rds')
#saveRDS(multi3,'multi3.rds')
saveRDS(multi4,'multi4.rds')


