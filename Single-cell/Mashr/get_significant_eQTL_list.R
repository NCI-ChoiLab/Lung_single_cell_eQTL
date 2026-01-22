#Author: Thong Luong
#Date: Jan 3rd 2024

setwd('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos')

#list of abbreviated cell type names
celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

myobject = list()

#extract significant eQTLs and make list
for (i in celltypes){
  myobject[[i]] = subset(read.table(paste0('./',i,'/',i,'.cis_qtl.txt'), sep = '\t', header = T), qval < 0.05)
}

saveRDS(myobject, file = 'full_sig_list.rds')

