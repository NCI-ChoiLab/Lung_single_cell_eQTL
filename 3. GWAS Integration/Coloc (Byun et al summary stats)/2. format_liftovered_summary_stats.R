#Author: Thong Luong
#Date: Jan 16th 2025

library(stringr)
library(dplyr)

setwd('/data/Choi_lung/TTL/Colocalization/Byun_gwas/')

total_lung = read.table('/data/Choi_lung/GWAS_Summary_Stats/Byun_et_al/GCST90134661_buildGRCh37.tsv', header = T, sep = '\t')
#luad = read.table('/data/Choi_lung/GWAS_Summary_Stats/Byun_et_al/LUAD/NG2022-ADE.meta.txt', header = T, sep = '\t')
#lusc = read.table('/data/Choi_lung/GWAS_Summary_Stats/Byun_et_al/LUSC/NG2022-SQC.meta.txt', header = T, sep = '\t')
#scc = read.table('/data/Choi_lung/GWAS_Summary_Stats/Byun_et_al/SCC/NG2022-SCC.meta.txt', header = T, sep = '\t')

liftover_function = function(m,l,f){
  liftover_snps = read.table(l)
  liftover_snps = liftover_snps[,c(1,3,4)]
  colnames(liftover_snps) = c('CHR','POS','SNP')
  
  meta = m
  meta$SNP = paste(meta$variant_id,meta$effect_allele, sep ='_')
  meta = subset(meta, SNP %in% liftover_snps$SNP)
  meta$chromosome = NULL
  meta$base_pair = NULL
  
  meta = merge(meta, liftover_snps, by = 'SNP')
  meta$SNP = NULL
  meta = meta[,c('variant_id','CHR','POS','effect_allele','other_allele','odds_ratio','standard_error','p_value')]
  write.table(meta, file = f, sep = '\t', col.names = T, row.names = F, quote = F)
}

liftover_function(total_lung,'hg38_total.bed','./hg38_total/hg38_total_meta.txt')

liftover_function_2 = function(m,l,f){
  liftover_snps = read.table(l)
  liftover_snps = liftover_snps[,c(1,3,4)]
  colnames(liftover_snps) = c('CHR','POS','SNP')
  
  meta = m
  meta$ref = str_split_fixed(meta$SNP,':',4)[,3]
  meta$alt = str_split_fixed(meta$SNP,':',4)[,4]
  meta$SNP = paste(meta$RS,meta$alt, sep ='_')
  meta = subset(meta, SNP %in% liftover_snps$SNP)
  meta$CHR = NULL
  meta$BP = NULL
  
  meta = merge(meta, liftover_snps, by = 'SNP')
  meta$SNP = NULL
  meta_rs =  meta[grepl('rs',meta$RS),]
  meta_not = meta[!grepl('rs',meta$RS),]
  meta_not$RS = paste(meta_not$CHR,meta_not$POS, sep = ':')
  meta = rbind(meta_not,meta_rs)
  meta = meta[,c('RS','CHR','POS','ref','alt','OR_FE','STD_FE','P_FE')]
  write.table(meta, file = f, sep = '\t', col.names = T, row.names = F, quote = F)
}

#liftover_function_2(luad,'hg38_luad.bed','./hg38_luad/hg38_luad_meta.txt')
#liftover_function_2(lusc,'hg38_lusc.bed','./hg38_lusc/hg38_lusc_meta.txt')
#liftover_function_2(scc,'hg38_scc.bed','./hg38_scc/hg38_scc_meta.txt')
