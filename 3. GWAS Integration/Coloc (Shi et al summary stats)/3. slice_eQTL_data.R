#Author: Thong Luong
#Date: Jan 6th 2025
#Modified Jan 13th 2025 to include sig slice

setwd('/data/Choi_lung/TTL/files_for_JC_2501/')

library(dplyr)

gene_loc = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')
snps = readRDS('/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final_freq.rds')
colnames(snps)[1] = 'variant_id'
snps$af = NULL
eGenes = readRDS('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/full_sig_list.rds')
eGenes = sapply(eGenes, `[[`, 1)

locus = c('2p11.2','2p14','2p23.3','3q22.3','3q26.2','3q28',
          '4p13','4q32.1','4q32.2','5p15.33','6p21.33','6p21.32',
          '6p21.1','6p12.1','6q22.1','7q31.33','9p21.3','10q25.2',
          '10q26.13','11q12.2','11p13','11q23.3','12q13.13','14q13.2',
          '15q21.2','15q21.3','16q23.3','17q24.2','18q12.1','19p13.3')

chr = c('chr2','chr2','chr2','chr3','chr3','chr3',
        'chr4','chr4','chr4','chr5','chr6','chr6',
        'chr6','chr6','chr6','chr7','chr9','chr10',
        'chr10','chr11','chr11','chr11','chr12','chr14',
        'chr15','chr15','chr16','chr17','chr18','chr19')

pos = c(85666618,65268924,25534840,138851169,169764547,189636338,
        44172387,156973740,163148970,1286401,30801788,32606581,
        41515652,53525197,117464145,124733330,22160088,112749531,
        124635640,61814184,34507219,118237616,51954475,34823979,
        49465269,56162025,82119933,67964738,32342958,725066)

tensor_path = '/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/'

for (i in ct){
  
  output = read.table(paste0(tensor_path,'output/',i,'/qtl_results_all.txt'), sep = '\t', header = T)
  snps_sub = subset(snps, variant_id %in% output$variant_id)
  output = merge(output,gene_loc, by = 'phenotype_id')
  output = merge(output,snps_sub, by = 'variant_id')
  output = output[,c('phenotype_id','phenotype_name','phenotype_chr','phenotype_pos','rsid','variant_chr','variant_pos','variant_id','ref','alt','start_distance','af','ma_samples','ma_count',
                     'pval_nominal','slope','slope_se','celltype')]
  
  
  for (j in 1:length(locus)){
    dir.create(paste0('sliced_Shi/',locus[j]), showWarnings = F)
    output_save = subset(output, variant_chr == chr[j])
    output_save = subset(output_save, variant_pos >= (pos[j] - 100000) & variant_pos <= (pos[j] + 100000))
    
    write.table(output_save, file = paste0('sliced_Shi/',locus[j],'/',i,"_",locus[j],"_100kb.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    
    dir.create(paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Shi_sig/',locus[j]), showWarnings = F)
    output_save_sig = subset(output_save, phenotype_id %in% eGenes)
    write.table(output_save_sig, file = paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Shi_sig/',locus[j],'/',i,"_",locus[j],"_100kb.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
  }
  
}




