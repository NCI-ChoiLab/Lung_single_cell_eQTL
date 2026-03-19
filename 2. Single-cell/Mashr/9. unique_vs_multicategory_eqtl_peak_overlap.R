#Author: Thong Luong
#Date: May 29th 2025

library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Overlap_eQTL_cCRE/Shared_Unique/Mashr')

#load mashr defined unique eQTLs
unique_eqtl = readRDS('unique.rds')

gene = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
gene = gene[,c('phenotype_id','phenotype_name')]

snps = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/snps.rds')
snps$variant_id = paste(snps$V1, snps$V4, snps$V6, snps$V5, sep = '_')
snps = snps[,c('V1','V4','V2',"variant_id")]
colnames(snps) = c('variant_chr','POS','RSID','variant_id')

#Load multiome ATAC peaks
peak_info = read.csv('../peak_info_w_annot.tsv', sep = '\t')

#overlap variants to unique peaks
get_overlap_unique = function(c,cm){
  cell = subset(unique_eqtl, celltype == c)
  sub_snps = subset(snps, variant_id %in% cell$variant_id)
  sub_snps$colocalized_cCRE = NA
  
  #inclusive
  #peak_info$celltype = grepl(cm,peak_info$peak_called_in)
  #peak_info = subset(peak_info, celltype == 'TRUE')
  
  #exclusive
  peak_info_cs = subset(peak_info, peak_called_in == cm)
  
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_cs$peak[which(peak_info_cs$seqnames == sub_snps$variant_chr[i] & 
                                  as.numeric(peak_info_cs$start) <= as.numeric(sub_snps$POS[i]) &
                                  as.numeric(peak_info_cs$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  uni = sub_snps[complete.cases(sub_snps),]
  cell_uni = subset(cell, variant_id %in% uni$variant_id)
  cell_uni = merge(cell_uni, uni, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_uni$phenotype_id)
  cell_uni = merge(cell_uni, gene_sub, by = 'phenotype_id')
  cell_uni = cell_uni[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_uni, paste0('./top4_unique/',c,'_uni_specific.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
  
  #multi
  peak_info_multi = subset(peak_info, categories == 4)
  
  sub_snps$colocalized_cCRE = NA
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_multi$peak[which(peak_info_multi$seqnames == sub_snps$variant_chr[i] & 
                                     as.numeric(peak_info_multi$start) <= as.numeric(sub_snps$POS[i]) &
                                     as.numeric(peak_info_multi$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  multi = sub_snps[complete.cases(sub_snps),]
  cell_multi = subset(cell, variant_id %in% multi$variant_id)
  cell_multi = merge(cell_multi, multi, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_multi$phenotype_id)
  cell_multi = merge(cell_multi, gene_sub, by = 'phenotype_id')
  cell_multi = cell_multi[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_multi, paste0('./top4_unique/',c,'_multi.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
  
  #any
  peak_info_any = peak_info
  
  sub_snps$colocalized_cCRE = NA
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_any$peak[which(peak_info_any$seqnames == sub_snps$variant_chr[i] & 
                                        as.numeric(peak_info_any$start) <= as.numeric(sub_snps$POS[i]) &
                                        as.numeric(peak_info_any$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  any = sub_snps[complete.cases(sub_snps),]
  cell_any = subset(cell, variant_id %in% any$variant_id)
  cell_any = merge(cell_any, any, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_any$phenotype_id)
  cell_any = merge(cell_any, gene_sub, by = 'phenotype_id')
  cell_any = cell_any[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_any, paste0('./top4_unique/',c,'_any.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
}

#no ec_aero_cap or peri_fib or ec_ven_sus or sub_fib
get_overlap_unique('AT1','AT1')
get_overlap_unique('AT2','AT2')
get_overlap_unique('Club','Club')
get_overlap_unique('Goblet','Goblet')
get_overlap_unique('Multiciliated','Ciliated')
get_overlap_unique('Alv_mph','Macrophage')
get_overlap_unique('Mono_mph','Macrophage')
get_overlap_unique('Cla_mono','Monocyte')
get_overlap_unique('CD4','T')
get_overlap_unique('CD8','T')
get_overlap_unique('NK','NK')

get_overlap_unique('Noncla_mono','Monocyte')
get_overlap_unique('DC2','DC')
get_overlap_unique('DC1','DC')
get_overlap_unique('Lym_EC_mat','Lymphatic')
get_overlap_unique('Lym_EC_pro','Lymphatic')
get_overlap_unique('EC_art','Artery')
get_overlap_unique('EC_ven_pul','Vein')
get_overlap_unique('EC_gen_cap', 'Capillary')
get_overlap_unique('Adv_fib','Fibroblast')
get_overlap_unique('Alv_fib','Fibroblast')
get_overlap_unique('Alv_trans','AT1_AT2')

#load mashr defined multicategory eQTLs
multi_eqtl = readRDS('multi4.rds')

#overlap variants to multicategory peaks
get_overlap_multi = function(c,cm){
  cell = subset(multi_eqtl, celltype == c)
  sub_snps = subset(snps, variant_id %in% cell$variant_id)
  sub_snps$colocalized_cCRE = NA
  
  #inclusive
  #peak_info$celltype = grepl(cm,peak_info$peak_called_in)
  #peak_info = subset(peak_info, celltype == 'TRUE')
  
  #exclusive
  peak_info_cs = subset(peak_info, peak_called_in == cm)
  
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_cs$peak[which(peak_info_cs$seqnames == sub_snps$variant_chr[i] & 
                                     as.numeric(peak_info_cs$start) <= as.numeric(sub_snps$POS[i]) &
                                     as.numeric(peak_info_cs$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  uni = sub_snps[complete.cases(sub_snps),]
  cell_uni = subset(cell, variant_id %in% uni$variant_id)
  cell_uni = merge(cell_uni, uni, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_uni$phenotype_id)
  cell_uni = merge(cell_uni, gene_sub, by = 'phenotype_id')
  cell_uni = cell_uni[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_uni, paste0('./top4_multi/',c,'_uni_specific.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
  
  #multi
  peak_info_multi = subset(peak_info, categories == 4)
  
  sub_snps$colocalized_cCRE = NA
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_multi$peak[which(peak_info_multi$seqnames == sub_snps$variant_chr[i] & 
                                        as.numeric(peak_info_multi$start) <= as.numeric(sub_snps$POS[i]) &
                                        as.numeric(peak_info_multi$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  multi = sub_snps[complete.cases(sub_snps),]
  cell_multi = subset(cell, variant_id %in% multi$variant_id)
  cell_multi = merge(cell_multi, multi, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_multi$phenotype_id)
  cell_multi = merge(cell_multi, gene_sub, by = 'phenotype_id')
  cell_multi = cell_multi[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_multi, paste0('./top4_multi/',c,'_multi.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
  
  #any
  peak_info_any = peak_info
  
  sub_snps$colocalized_cCRE = NA
  for (i in 1:nrow(sub_snps)){
    cCRE = peak_info_any$peak[which(peak_info_any$seqnames == sub_snps$variant_chr[i] & 
                                      as.numeric(peak_info_any$start) <= as.numeric(sub_snps$POS[i]) &
                                      as.numeric(peak_info_any$end) >= as.numeric(sub_snps$POS[i])
    )]
    if (length(cCRE) == 1) {
      sub_snps$colocalized_cCRE[i] = cCRE
    }
  }
  
  any = sub_snps[complete.cases(sub_snps),]
  cell_any = subset(cell, variant_id %in% any$variant_id)
  cell_any = merge(cell_any, any, by ='variant_id')
  gene_sub = subset(gene, phenotype_id %in% cell_any$phenotype_id)
  cell_any = merge(cell_any, gene_sub, by = 'phenotype_id')
  cell_any = cell_any[,c("phenotype_id","phenotype_name","variant_id","RSID","start_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se","celltype","colocalized_cCRE")]
  write.table(cell_any, paste0('./top4_multi/',c,'_any.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
}

get_overlap_multi('AT1','AT1')
get_overlap_multi('AT2','AT2')
get_overlap_multi('Club','Club')
get_overlap_multi('Goblet','Goblet')
get_overlap_multi('Multiciliated','Ciliated')
get_overlap_multi('Alv_mph','Macrophage')
get_overlap_multi('Mono_mph','Macrophage')
get_overlap_multi('Cla_mono','Monocyte')
get_overlap_multi('CD4','T')
get_overlap_multi('CD8','T')
get_overlap_multi('NK','NK')

get_overlap_multi('Noncla_mono','Monocyte')
get_overlap_multi('DC2','DC')
get_overlap_multi('DC1','DC')
get_overlap_multi('Lym_EC_mat','Lymphatic')
get_overlap_multi('Lym_EC_pro','Lymphatic')
get_overlap_multi('EC_art','Artery')
get_overlap_multi('EC_ven_pul','Vein')
get_overlap_multi('EC_gen_cap', 'Capillary')
get_overlap_multi('Adv_fib','Fibroblast')
get_overlap_multi('Alv_fib','Fibroblast')
get_overlap_multi('Alv_trans','AT1_AT2')