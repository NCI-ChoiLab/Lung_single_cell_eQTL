#Author: Thong Luong
#Date: December 4th 2025

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/')

files = list.files('Q_pseudotime_output/',full.names = T)

full = data.frame()
for (i in 1:length(files)){
  t = read.csv(files[i])
  t = t[16,c("term","gene","snp","zvalue","pval")]
  full = rbind(full,t)
  rm(t)
}

full$fdr = p.adjust(full$pval, method = 'fdr')

t = read.csv('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/NBME_conQ_pseudotime_output_rmInficells.csv')

gene = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
gene = gene[,c("phenotype_id","phenotype_name")]
colnames(gene) = c('ID','gene')
full = merge(full, gene, by = 'gene')
full = full[,c("ID","gene","snp","zvalue","pval","fdr")]
colnames(full) = colnames(t)

write.csv(full, 'Q_pseudotime_output/NBME_top_snps.csv', row.names = F)
