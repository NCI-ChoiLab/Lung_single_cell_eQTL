#Author: Thong Luong
#Date: April 21st 2025

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

eqtl_heatmap = function(m,l){
  
  lfsr.all = lfsrs(m)
  pm.mash = betas(m)
  standard.error = errors(m)
  
  pm.mash.beta <- pm.mash*standard.error
  pm.mash.beta=pm.mash.beta[l,]
  lfsr.mash=lfsr.all[l,]
  
  thresh=0.05
  shared.fold.size=matrix(NA,nrow = ncol(lfsr.mash),ncol=ncol(lfsr.mash))
  colnames(shared.fold.size)=rownames(shared.fold.size)=colnames(standard.error)
  for(i in 1:ncol(lfsr.mash)){
    for(j in 1:ncol(lfsr.mash)){
      sig.row=which(lfsr.mash[,i]<thresh)
      sig.col=which(lfsr.mash[,j]<thresh)
      a=(union(sig.row,sig.col))
      quotient=(pm.mash.beta[a,i]/pm.mash.beta[a,j])
      shared.fold.size[i,j]=mean(quotient > 0.5 & quotient < 2)
    }
  }
  shared.fold.size = shared.fold.size[c(epithelial, immune, endothelial, stromal),
                                      c(epithelial, immune, endothelial, stromal)]
  
  return(shared.fold.size)
}

clrs <- colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF",
                               "#E0F3F8","#91BFDB","#4575B4")))(64)


pushViewport(viewport(gp = gpar(fontfamily = "Arial")))

hm = Heatmap(top_sig_mashr, col = clrs, name = 'Proportion Shared', 
        column_names_rot = 45, show_row_dend = F, row_names_side = 'left', width = unit(70,'mm'), height = unit(65,'mm'),
        row_names_gp = gpar(fontsize = 6, col = c(rep('green',8),rep('blue',13),rep('purple',7),rep('orange',5))),
        column_names_gp = gpar(fontsize = 6, col = c(rep('green',8),rep('blue',13),rep('purple',7),rep('orange',5))),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize=6),
                                    grid_width = unit(2, "mm"),
                                    grid_height = unit(4, "mm")))
ht = draw(hm, newpage = FALSE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)



png('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/Figure_3e.png', w = w, h = h, units = 'in', res = 450)
Heatmap(top_sig_mashr, col = clrs, name = 'Proportion Shared', 
        column_names_rot = 45, show_row_dend = F, row_names_side = 'left', width = unit(70,'mm'), height = unit(65,'mm'),
        row_names_gp = gpar(fontsize = 6, col = c(rep('green',8),rep('blue',13),rep('purple',7),rep('orange',5))),
        column_names_gp = gpar(fontsize = 6, col = c(rep('green',8),rep('blue',13),rep('purple',7),rep('orange',5))),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 6),title_gp = gpar(fontsize=6),
                                    grid_width = unit(2, "mm"),
                                    grid_height = unit(4, "mm")))
dev.off()



