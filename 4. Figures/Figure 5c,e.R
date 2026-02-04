#Author: Thong Luong
#Date: August 19th 2025

library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(reshape2)
library(Matrix)
library(Hmisc)
library(wesanderson)
library(magrittr)
library(rtracklayer)
library(geomtextpath)
library(snpStats)
library(Seurat)
library(scales)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/boxplot_fig/')

fam <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.fam"
bim <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bim"
bed <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bed"
sample <- read.plink(bed, bim, fam)
genotype_mtx <- sample$genotypes
genotypes<- genotype_mtx@.Data

snp_info = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/snps.rds')
rownames(snp_info) = paste(snp_info$V1,snp_info$V4,snp_info$V6,snp_info$V5,sep='_')

time_finite = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/pseudotime_finite.rds')

sc = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/total_with_preinfo.rds')
sc[['RNA']] = JoinLayers(sc[['RNA']])
Idents(sc) = 'cell_types'
sc = subset(sc, cell_types %in% c('AT2','AT1','Alveolar Transitional Cells'))

cellInfo = sc@meta.data
cellInfo = cellInfo[rownames(time_finite),]

counts <- sc@assays$RNA$counts

#only keeping genes that is expressed in at least 10% of all cells
counts = counts[rowSums(counts >0) >= (ncol(counts)/10),]
#keeping genes with mean counts of greater than 0.1
counts = counts[rowMeans(counts) > 0.1,]

norm_counts = NormalizeData(counts)

#new_dataFrame <- apply(norm_counts, 1, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
#new_dataFrame = as.data.frame(new_dataFrame)
#new_dataFrame_sub = new_dataFrame[rownames(time_finite),]

norm_counts_sub = as.data.frame(t(norm_counts))
norm_counts_sub = norm_counts_sub[rownames(time_finite),]

genotypes_sub = as.data.frame(genotypes)
genotypes_sub$Sample = rownames(genotypes_sub)

int_effect = readRDS('../Q_pseudotime_output/int_eqtl_effect.rds')
get_sig_time = function(d){
  tested_df = d
  tested_df = as.data.frame(t(tested_df))
  tested_df$Pseudotime = c(1,2,3,4,5,6)
  df = as.data.frame(matrix(nrow = ncol(tested_df)-1, ncol = 2))
  rownames(df) = colnames(tested_df)[1:(ncol(tested_df)-1)]
  colnames(df) = c('r.squared','pval')
  df_quad = df
  
  for (i in 1:(ncol(tested_df)-1)){
    linear = lm(formula = tested_df[,i] ~ Pseudotime, data = tested_df)
    df[colnames(tested_df)[i],'r.squared'] = summary(linear)$r.squared
    df[colnames(tested_df)[i],'pval'] = summary(linear)$coefficients[,4][['Pseudotime']]
  }
  
  sig_ln_df = subset(df, pval < 0.05)
  print(nrow(sig_ln_df))
  
  for (i in 1:(ncol(tested_df)-1)){
    quad = lm(formula = tested_df[,i] ~ Pseudotime + I(Pseudotime^2), data = tested_df)
    df_quad[colnames(tested_df)[i],'r.squared'] = summary(quad)$r.squared
    df_quad[colnames(tested_df)[i],'pval'] = summary(quad)$coefficients[,4][['Pseudotime']]
  }
  
  sig_quad_df = subset(df_quad, pval < 0.05)
  print(nrow(sig_quad_df))
  
  both = intersect(rownames(sig_quad_df),rownames(sig_ln_df))
  linear_name = setdiff(rownames(sig_ln_df),both)
  quad_name = setdiff(rownames(sig_quad_df),both)
  signif_name = unique(c(rownames(sig_quad_df),rownames(sig_ln_df)))
  
  sig_ln = as.data.frame(subset(d, rownames(d) %in% linear_name))
  sig_quad = as.data.frame(subset(d, rownames(d) %in% quad_name))
  both_df = as.data.frame(subset(d, rownames(d) %in% both))
  non_sig = as.data.frame(subset(d, !(rownames(d) %in% signif_name)))
  
  sig_ln$Model = 'Linear'
  sig_quad$Model = 'Quadratic'
  both_df$Model = 'Both'
  non_sig$Model = 'Neither'
  print(nrow(both_df))
  print(nrow(non_sig))
  
  tested_df = rbind(sig_ln,both_df)
  tested_df = rbind(tested_df,sig_quad)
  tested_df = rbind(tested_df,non_sig)
  return(tested_df)
}

int_effect = get_sig_time(int_effect)
int_effect = as.data.frame(int_effect)
int_effect$Gene = str_split_fixed(rownames(int_effect),pattern = '\\|',2)[,1]
int_effect$variant_id = str_split_fixed(rownames(int_effect),pattern = '\\|',2)[,2]
int_effect = subset(int_effect, Model == 'Linear')
int_effect = subset(int_effect, Gene %in% colnames(norm_counts_sub))


variant_loc = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/snps.rds')
colnames(variant_loc)[2] = 'rsid'
variant_loc$variant_id = paste(variant_loc$V1,variant_loc$V4,variant_loc$V6,variant_loc$V5,sep = '_')
variant_loc = variant_loc[,c('variant_id','rsid')]
int_effect = merge(int_effect,variant_loc,by = 'variant_id')
int_effect = int_effect[order(int_effect$Q6, decreasing = T),]

generate_graph = function(goi, soi){
  df = cbind(time_finite, norm_counts_sub[,goi])
  colnames(df)[4] = goi
  df = cbind(df, cellInfo[,'Sample'])
  colnames(df)[5] = 'Sample'
  
  chr_pos = rownames(snp_info[which(snp_info$V2 == soi),])
  
  geno_sub = genotypes_sub[,c('Sample',chr_pos)]
  rownames(geno_sub) = NULL
  df = merge(df, geno_sub, by = 'Sample')
  colnames(df)[5] = 'Gene'
  colnames(df)[6] = 'snp'
  
  df$snp_ra <- "0|0"
  df$snp_ra[which(df$snp == "01")] <- "1|1"
  df$snp_ra[which(df$snp == "02")] <- "0|1"
  
  df$snp_order = 1
  df$snp_order[which(df$snp_ra == '0|1')] = 2
  df$snp_order[which(df$snp_ra == '1|1')] = 3
  
  df_final = df
  
  df_final$Gene = scale(df_final$Gene, center = T, scale = T)
  
  df_final$q = factor(df_final$q, levels = c('Q1','Q2','Q3','Q4','Q5','Q6'))
  
  df_final$q_snp = paste(df_final$q,df_final$snp_ra,sep = '_')
  df_final$q_snp = factor(df_final$q_snp, levels = c('Q1_0|0','Q1_0|1','Q1_1|1',
                                         'Q2_0|0','Q2_0|1','Q2_1|1',
                                         'Q3_0|0','Q3_0|1','Q3_1|1',
                                         'Q4_0|0','Q4_0|1','Q4_1|1',
                                         'Q5_0|0','Q5_0|1','Q5_1|1',
                                         'Q6_0|0','Q6_0|1','Q6_1|1'))
  
  g1 = ggplot(df_final, aes(fill=snp_ra, y=Gene, x=q_snp)) + 
    geom_violin(position = position_dodge(.5), trim=FALSE,linewidth = 0.25) + 
    geom_boxplot(position = position_dodge(.5), aes(middle = median(Gene)),width=0.1, linewidth = 0.25, outlier.size = .1) + 
    geom_smooth( mapping = aes(x = q_snp, y = Gene, group = q),formula = y~x, color = "gray", method='lm', size = .5, se =T,fill = alpha("gray", .5)) + 
    scale_x_discrete(labels = c('','Q1','','','Q2','','','Q3','','','Q4','','','Q5','','','Q6','' )) +
    labs(title="",x='', y = paste0("Scaled expression of ", goi))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(goi, ": ",soi, ": ", paste(chr_pos, "b38", sep = "_")))+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 6, color = 'black', face = 'plain'),
          axis.ticks.length.x = unit(rep(c(0,.1,0),6),"cm"),
          axis.text.x = element_text(color = "black", size = 5,face = "plain"),
          axis.text.y = element_text(color = "black", size = 5,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 5, face = "plain"),
          axis.title.y = element_text(color = "black", size = 5, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
  
          legend.position="none")
  return(g1) 
  
}

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_5/'

scgb3a2 = generate_graph('SCGB3A2','rs138215727')
ggsave(filename = file.path(plot_save_dir, 'scgb3a2_s.png'), plot = scgb3a2, width = 75, height = 50, units = "mm", dpi = 450, device = 'png')
akr1c1 = generate_graph('AKR1C1','rs11252856')
ggsave(filename = file.path(plot_save_dir, 'akr1c1_s.png'), plot = akr1c1, width = 75, height = 50, units = "mm", dpi = 450, device = 'png')
fam13a = generate_graph('FAM13A','rs1344312507')
ggsave(filename = file.path(plot_save_dir, 'fam13a_s.png'), plot = fam13a, width = 75, height = 50, units = "mm", dpi = 450, device = 'png')

#Figure 5e
generate_graph_sf = function(goi, soi){
  df = cbind(time_finite, norm_counts_sub[,goi])
  colnames(df)[4] = goi
  df = cbind(df, cellInfo[,'Sample'])
  colnames(df)[5] = 'Sample'
  
  chr_pos = rownames(snp_info[which(snp_info$V2 == soi),])
  
  geno_sub = genotypes_sub[,c('Sample',chr_pos)]
  rownames(geno_sub) = NULL
  df = merge(df, geno_sub, by = 'Sample')
  colnames(df)[5] = 'Gene'
  colnames(df)[6] = 'snp'
  
  df$snp_ra <- "0|0"
  df$snp_ra[which(df$snp == "01")] <- "1|1"
  df$snp_ra[which(df$snp == "02")] <- "0|1"
  
  df$snp_order = 1
  df$snp_order[which(df$snp_ra == '0|1')] = 2
  df$snp_order[which(df$snp_ra == '1|1')] = 3
  
  df_final = df
  
  df_final$Gene = scale(df_final$Gene, center = T, scale = T)
  
  df_final$q = factor(df_final$q, levels = c('Q1','Q2','Q3','Q4','Q5','Q6'))
  
  df_final$q_snp = paste(df_final$q,df_final$snp_ra,sep = '_')
  df_final$q_snp = factor(df_final$q_snp, levels = c('Q1_0|0','Q1_0|1','Q1_1|1',
                                                     'Q2_0|0','Q2_0|1','Q2_1|1',
                                                     'Q3_0|0','Q3_0|1','Q3_1|1',
                                                     'Q4_0|0','Q4_0|1','Q4_1|1',
                                                     'Q5_0|0','Q5_0|1','Q5_1|1',
                                                     'Q6_0|0','Q6_0|1','Q6_1|1'))
  
  g1 = ggplot(df_final, aes(fill=snp_ra, y=Gene, x=q_snp)) + 
    geom_violin(position = position_dodge(.5), trim=FALSE,linewidth = 0.25) + 
    geom_boxplot(position = position_dodge(.5), aes(middle = median(Gene)),width=0.1, linewidth = 0.25, outlier.size = .1) + 
    geom_smooth( mapping = aes(x = q_snp, y = Gene, group = q),formula = y~x, color = "gray", method='lm', size = .5, se =T,fill = alpha("gray", .5)) + 
    scale_x_discrete(labels = c('','Q1','','','Q2','','','Q3','','','Q4','','','Q5','','','Q6','' )) +
    labs(title="",x='', y = paste0("Scaled expression of ", goi))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(goi, ": ",soi, ": ", paste(chr_pos, "b38", sep = "_")))+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 7, color = 'black', face = 'plain'),
          axis.ticks.length.x = unit(rep(c(0,.1,0),6),"cm"),
          axis.text.x = element_text(color = "black", size = 6,face = "plain"),
          axis.text.y = element_text(color = "black", size = 6,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 6, face = "plain"),
          axis.title.y = element_text(color = "black", size = 6, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          
          legend.position="none")
  return(g1) 
  
}

sftpa2 = generate_graph_sf('SFTPA2','rs1617662')
ggsave(filename = file.path(plot_save_dir, 'sftpa2.png'), plot = sftpa2, width = 90, height = 60, units = "mm", dpi = 450, device = 'png')







