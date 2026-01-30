#Author: Thong Luong
#Date: April 25th 2025

library(QTLExperiment)
library(dplyr)
library(lattice)
library(ComplexHeatmap)
library(circlize)
library(UpSetR)
library(ComplexUpset)
library(stringr)
library(Seurat)
library(ggrepel)

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
lfsr$Genes = str_split_fixed(rownames(lfsr),'\\|',2)[,1]

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
upset(fromList(lfsr_list), celltypes, max_degree = 1)

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
cell_type_specific = subset(interaction_list, int %in% celltypes)

celltype_specific = data.frame(matrix(vector(), length(unique(cell_type_specific$int)),2, 
                         dimnames = list(unique(cell_type_specific$int), c('Cell_type_specific', 'Sig_eQTL'))),
                  stringsAsFactors = F)

for (c in rownames(celltype_specific)){
  cell = subset(cell_type_specific, int == c)
  celltype_specific[c,'Cell_type_specific'] = nrow(cell)
  celltype_specific[c,"Sig_eQTL"] = length(lfsr_list[[c]])
}

celltype_specific$log_sig = log10(celltype_specific$Sig_eQTL)

print(setdiff(rownames(data),rownames(celltype_specific)))

col_specific <- list('Multiciliated' = 'darkseagreen1' ,'AT2'='green4', 'AT1'='green1','Club' = 'seagreen4','Alv_trans' = 'olivedrab', 'Sec_trans' = 'palegreen2', 'Goblet' = 'limegreen','Basal' = 'darkgreen',
             'Alv_mph' = 'royalblue1', 'NK' = 'dodgerblue4', 'Cla_mono'='steelblue4','Noncla_mono'='blue','CD4'='aquamarine4','Mono_mph'='lightblue', 'CD8'='darkturquoise', 
             'DC2'='cornflowerblue', 'Int_mph_peri'='dodgerblue', 'DC1'='skyblue', 'Mig_DC'='cadetblue',
             'Lym_EC_mat'='purple1','EC_art'='orchid1','EC_ven_pul'='mediumpurple3','Lym_EC_pro'='darkviolet','EC_gen_cap'='magenta',
              'Adv_fib'='lightsalmon','Alv_fib'='sienna')
celltype_specific = celltype_specific[c('Multiciliated','AT2','AT1','Club','Alv_trans','Sec_trans','Goblet','Basal',
              'Alv_mph','NK','Cla_mono','Noncla_mono','CD4','Mono_mph','CD8','DC2','Int_mph_peri','DC1','Mig_DC',
              'Lym_EC_mat','EC_art','EC_ven_pul','Lym_EC_pro','EC_gen_cap',
              'Adv_fib','Alv_fib'),]

t = ggplot(celltype_specific,aes(x=Sig_eQTL,y=Cell_type_specific)) + geom_point(aes(color = col_specific)) + 
  geom_label_repel(aes(label = rownames(celltype_specific)),
                   box.padding   = 0.3, 
                   point.padding = 0.2,
                   label.padding = .1,
                   segment.color = 'grey50', size = 1.75) +
  xlab('# of Significant eQTL (lfsr < 0.05)') +
  ylab('# Cell Type Specific eQTL') + 
  theme(text = element_text(family = 'Arial'),
        axis.text=element_text(size=5),
        axis.text.x = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color= cols),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1)) 

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'sup_fig8b.png'), plot = t, width = 85, height = 80, units = "mm", dpi = 450, device = 'png',)


lfsr$Genes = str_split_fixed(rownames(lfsr),'\\|',2)[,1]

get_gene_list = function(l){
  lfsr_list = list()
  for (c in l){
    cell = lfsr[,c(c,'Genes')]
    colnames(cell)[1] = 'lfsr'
    cell = subset(cell, lfsr < 0.05)
    lfsr_list[[c]] = unique(cell$Genes)
  }
  return(lfsr_list)
}

eGene_list = get_gene_list(celltypes)

data = data.frame(matrix(vector(), length(celltypes),3, 
                         dimnames = list(celltypes, c('Cell_number', 'eGenes','Mashr_eGenes'))),
                  stringsAsFactors = F)

sc = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/total_with_preinfo.rds')
numb = as.data.frame(table(sc$cell_types))
rownames(numb) = numb$Var1
numb$Var1 = NULL
rownames(numb) = c('Adv_fib','Alv_fib','Alv_mph','Alv_trans','AT1','AT2','Bcells','Basal',
                   'CD4','CD8','Cla_mono','Club','DC1','DC2','EC_aero_cap','EC_art',
                   'EC_gen_cap','EC_ven_pul','EC_ven_sys','Goblet','Int_mph_peri','Lym_EC_mat','Lym_EC_pro','Mast',
                   "Mesothelium",'Mig_DC','Mono_mph','Multiciliated','Myo_fib','Neuroendocrine','NK','Noncla_mono',
                   'Peri_fib','Plasma','Pla_DC','Pro_mph','Pro_NK','Pro_T','Sec_trans','SM','Sub_fib')



for (c in celltypes){
  data[c,'Cell_number'] = numb[c,'Freq']
  egenes = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  egenes = length(unique(egenes$phenotype_id))
  data[c,'eGenes'] = egenes
  data[c,'Mashr_eGenes'] = length(eGene_list[[c]])
}

cols <- list('Multiciliated' = 'darkseagreen1' ,'AT2'='green4', 'AT1'='green1','Club' = 'seagreen4','Alv_trans' = 'olivedrab', 'Sec_trans' = 'palegreen2', 'Goblet' = 'limegreen','Basal' = 'darkgreen',
          'Alv_mph' = 'royalblue1', 'NK' = 'dodgerblue4', 'Cla_mono'='steelblue4','Noncla_mono'='blue','CD4'='aquamarine4','Mono_mph'='lightblue', 'CD8'='darkturquoise', 
          'DC2'='cornflowerblue', 'Int_mph_peri'='dodgerblue', 'DC1'='skyblue', 'Mig_DC'='cadetblue', 'Bcells' = 'midnightblue','Mast'='steelblue',
          'Lym_EC_mat'='purple1','EC_art'='orchid1','EC_ven_pul'='mediumpurple3','EC_ven_sys'='mediumorchid2','Lym_EC_pro'='darkviolet','EC_gen_cap'='magenta','EC_aero_cap'='orchid4', 
          'SM'='tan', 'Adv_fib'='lightsalmon','Alv_fib'='sienna', 'Peri_fib'='gold', 'Sub_fib'='goldenrod3')
data = data[c('Multiciliated','AT2','AT1','Club','Alv_trans','Sec_trans','Goblet','Basal',
                'Alv_mph','NK','Cla_mono','Noncla_mono','CD4','Mono_mph','CD8','DC2','Int_mph_peri','DC1','Mig_DC','Bcells','Mast',
                'Lym_EC_mat','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_pro','EC_gen_cap','EC_aero_cap',
                'SM','Adv_fib','Alv_fib','Peri_fib','Sub_fib'),]


tensor_scatterplot = ggplot(data,aes(x=Cell_number,y=eGenes)) + geom_point(aes(color = cols)) + 
  geom_label_repel(aes(label = rownames(data)),
                   box.padding   = 0.3, 
                   point.padding = 0.2,
                   label.padding = .1,
                   segment.color = 'grey50', size = 1.75) +
  xlab('Number of Cells') +
  ylab('# Significant eGenes (qvalue < 0.05)') + 
  theme(axis.text.x = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color= cols),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1)) 



plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'sup_fig7a.png'), plot = tensor_scatterplot, width = 85, height = 85, units = "mm", dpi = 450, device = 'png',)



mashr_scatterplot = ggplot(data,aes(x=Cell_number,y=Mashr_eGenes)) + geom_point(aes(color = cols)) + 
  geom_label_repel(aes(label = rownames(data)),
                   box.padding   = 0.3, 
                   point.padding = 0.2,
                   label.padding = .1,
                   segment.color = 'grey50', size = 1.75) +
  xlab('Number of Cells') +
  ylab('# Significant eGenes (lfsr < 0.05)') + 
  theme(axis.text.x = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color= cols),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1)) 



plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'sup_fig7b.png'), plot = mashr_scatterplot, width = 85, height = 85, units = "mm", dpi = 450, device = 'png',)


