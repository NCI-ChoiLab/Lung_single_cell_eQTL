#Author: Thong Luong
#Date: December 3rd 2025

library(snpStats)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/')

#load_gene_information_from_eqtl_analysis
gene = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
gene = gene[,c("phenotype_id","phenotype_name")]
colnames(gene)[2] = 'symbol' 
top_eqtl = readRDS('../final_tensor.rds')
top_eqtl = top_eqtl[,c("phenotype_id","variant_id")]
top_eqtl$eqtl = paste(top_eqtl$phenotype_id,top_eqtl$variant_id)
top_eqtl = top_eqtl[!duplicated(top_eqtl$eqtl),]
top_eqtl = merge(top_eqtl,gene,by = 'phenotype_id')
top_eqtl = top_eqtl[,c("phenotype_id","variant_id","symbol")]
write.table(top_eqtl,'eqtl_list_top.tsv',sep = '\t', col.names = T, row.names = F, quote = F)

#load final genotype file
fam <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.fam"
bim <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bim"
bed <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bed"
sample <- read.plink(bed, bim, fam)
genotype_mtx <- sample$genotypes
genotypes<- genotype_mtx@.Data
rm(genotype_mtx)
rm(fam)
rm(bim)
rm(bed)

save_snp_matrix = function(l,f){
  geno_sub = genotypes[,l$variant_id]
  char_mat <- apply(geno_sub, MARGIN = 1:2, FUN = rawToChar)
  char_mat[char_mat == '\001'] <- 2
  char_mat[char_mat == '\002'] <- 1
  char_mat[char_mat == '\003'] <- 0
  
  num_matrix <- apply(char_mat, 2, as.numeric)
  rownames(num_matrix) = rownames(char_mat)
  #rownames(num_matrix) = paste0('0:',rownames(char_mat))
  saveRDS(num_matrix,f)
}

save_snp_matrix(top_eqtl, 'top_eqtl_snps.rds')


