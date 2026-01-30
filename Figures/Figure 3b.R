#Author: Thong Luong
#Date: April 25th 2025

library(QTLExperiment)
library(dplyr)
library(lattice)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(ComplexUpset)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Mashr_heatmap/')

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
lfsr$eqtl = c(1:nrow(lfsr))

get_eqtl_list = function(l){
  lfsr_list = list()
  for (c in l){
    cell = lfsr[,c(c,'eqtl')]
    colnames(cell)[1] = 'lfsr'
    cell = subset(cell, lfsr < 0.05)
    lfsr_list[[c]] = cell$eqtl
  }
  return(lfsr_list)
}


lfsr_list = get_eqtl_list(celltypes)
uc = upset(fromList(lfsr_list), celltypes, max_degree = 1, width_ratio = .2,
           matrix=(
             intersection_matrix(
               geom=geom_point(shape=16, size=1) # Set dot size here
             )
           ))

#library(ggplot2)
plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/'
ggsave(filename = file.path(plot_save_dir, 'supp_upset.pdf'), plot = uc, width = 180, height = 150, units = "mm", dpi = 450, device = 'pdf',)



extract_interactions = function(l){
  df2 = data.frame(eqtl = unique(unlist(l)))
  df1 = lapply(l, function(x){
    data.frame(eqtl = x)
  }) %>% bind_rows(.id = 'path')
  
  df_int <- lapply(df2$eqtl,function(x){
    # pull the name of the intersections
    intersection <- df1 %>% 
      dplyr::filter(eqtl==x) %>% 
      arrange(path) %>% 
      pull("path") %>% 
      paste0(collapse = "|")
    
    # build the dataframe
    data.frame(eqtl = x,int = intersection)
  }) %>% 
    bind_rows()
  return(df_int)
}

interaction_list = extract_interactions(lfsr_list)

epi_list = get_eqtl_list(epithelial)
endo_list = get_eqtl_list(endothelial)
imm_list = get_eqtl_list(immune)
strom_list = get_eqtl_list(stromal)

Cell_group = list()

uniq <- sort(unique(unlist(epi_list)))
Cell_group[['Epithelial']] = uniq[rowSums(sapply(epi_list, `%in%`, x = uniq)) > 1]
epi_unique = subset(interaction_list, int %in% epithelial)
Cell_group[['Epithelial']] = c(unlist(Cell_group[['Epithelial']]),epi_unique$eqtl)

uniq <- sort(unique(unlist(endo_list)))
Cell_group[['Endothelial']] = uniq[rowSums(sapply(endo_list, `%in%`, x = uniq)) > 1]
endo_unique = subset(interaction_list, int %in% endothelial)
Cell_group[['Endothelial']] = c(unlist(Cell_group[['Endothelial']]),endo_unique$eqtl)

uniq <- sort(unique(unlist(imm_list)))
Cell_group[['Immune']] = uniq[rowSums(sapply(imm_list, `%in%`, x = uniq)) > 1]
imm_unique = subset(interaction_list, int %in% immune)
Cell_group[['Immune']] = c(unlist(Cell_group[['Immune']]),imm_unique$eqtl)

uniq <- sort(unique(unlist(strom_list)))
Cell_group[['Stromal']] = uniq[rowSums(sapply(strom_list, `%in%`, x = uniq)) > 1]
strom_unique = subset(interaction_list, int %in% stromal)
Cell_group[['Stromal']] = c(unlist(Cell_group[['Stromal']]),strom_unique$eqtl)
rm(uniq)


u = upset(fromList(Cell_group),c('Epithelial','Endothelial','Immune','Stromal'),name = 'Cell category', width_ratio = .1, max_degree = 1)

