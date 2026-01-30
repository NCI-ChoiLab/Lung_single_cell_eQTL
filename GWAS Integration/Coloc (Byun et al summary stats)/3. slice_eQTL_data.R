#Author: Thong Luong
#Date: Jan 16th 2025


setwd('/data/Choi_lung/TTL/files_for_JC_2501/')

library(dplyr)

gene_loc = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')
snps = readRDS('/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final_freq.rds')
colnames(snps)[1] = 'variant_id'
snps$af = NULL
eGenes = readRDS('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/full_sig_list.rds')
eGenes = sapply(eGenes, `[[`, 1)
eGenes = unique(unlist(eGenes))


locus = c('1p31.1','2p14','2p16.1','2q34','3p22.1','3p25.3','3q28','5p15.33_1','5p15.33_2','6p21.32_1','6p21.32_2',
          '6p21.33','6p22.1_1','6p22.1_2','6p22.1_3','6p25.3_1','6p25.3_2','6p27','6q21','6q22.1','8p12','8p21.1','8p21.2',
          '9p21.3_1','9p21.3_2','10q24.3','10q24.33','10q25.2_1','10q25.2_2','11q22.3','11q23.3_1','11q23.3_2','12p13.33',
          '13q13.1_1','13q13.1_2','15q21.1_1','15q21.1_2','15q25.1_1','15q25.1_2','19q13.2','20q13.33','22q12.1')

chr = c('chr1','chr2','chr2','chr2','chr3','chr3','chr3','chr5','chr5','chr6','chr6',
        'chr6','chr6','chr6','chr6','chr6','chr6','chr6','chr6','chr6','chr8','chr8','chr8',
        'chr9','chr9','chr10','chr10','chr10','chr10','chr11','chr11','chr11','chr12','chr13',
        'chr13','chr15','chr15','chr15','chr15','chr19','chr20','chr22')

pos = c(77984833,65262608,58858891,213134686,42875555,9928389,189618055,1287079,1340986,32604329,
        32433302,31459618,28965208,29813272,29639324,396321,410848,166962978,109418898,117456641,
        32546055,27487202,27536740,21775493,21763749,103927874,103934543,112733139,112728166,108326169,
        118222832,118237616,889653,32394673,32394413,49091284,49084427,78565644,78590583,40847202,63682701,28725099)

tensor_path = '/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/'

for (i in ct){
  
  output = read.table(paste0(tensor_path,'output/',i,'/qtl_results_all.txt'), sep = '\t', header = T)
  snps_sub = subset(snps, variant_id %in% output$variant_id)
  output = merge(output,gene_loc, by = 'phenotype_id')
  output = merge(output,snps_sub, by = 'variant_id')
  output = output[,c('phenotype_id','phenotype_name','phenotype_chr','phenotype_pos','rsid','variant_chr','variant_pos','variant_id','ref','alt','start_distance','af','ma_samples','ma_count',
                     'pval_nominal','slope','slope_se','celltype')]
  
  for (j in 1:length(locus)){
    dir.create(paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Byun_sig/',locus[j]), showWarnings = F)
    output_save = subset(output, variant_chr == chr[j])
    output_save = subset(output_save, variant_pos >= (pos[j] - 100000) & variant_pos <= (pos[j] + 100000))
    output_save_sig = subset(output_save, phenotype_id %in% eGenes)
    write.table(output_save_sig, file = paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Byun_sig/',locus[j],'/',i,"_",locus[j],"_100kb.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
  }
  
}




