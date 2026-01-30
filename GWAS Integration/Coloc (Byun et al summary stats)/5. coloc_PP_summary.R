#Author: Thong Luong
#Date: Dec 15th 2024

library(coloc)
library(foreach)
library(doParallel)
library(parallel)
library(stringr)

setwd('/data/Choi_lung/TTL/Colocalization/Byun_gwas/hg38_total/coloc/')

ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

get_PP = function(d){
  setwd(d)
  files = list.files(pattern = '_SNP_coloc.rds')
  
  
  ldf = lapply(files, readRDS)
  genes = names(ldf[[1]])
  
  for (i in 1:length(genes)){
    gene_data_frame = data.frame()
    for (j in 1:length(ldf)){
      a = ldf[[j]][[genes[i]]]$summary
      gene_data_frame = rbind(gene_data_frame, a)
    }
    rownames(gene_data_frame) = ct
    gene_data_frame = gene_data_frame[,c(2:6)]
    colnames(gene_data_frame) = c('PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    write.table(gene_data_frame, paste0(genes[i],'.tsv'), sep = '\t', row.names = T, col.names = T, quote = F)
  }
  
}

dir = dir()[file.info(dir())$isdir]

for (d in 1:length(dir)){
  setwd('/data/Choi_lung/TTL/Colocalization/Byun_gwas/hg38_total/coloc/')
  get_PP(dir[d])
}


