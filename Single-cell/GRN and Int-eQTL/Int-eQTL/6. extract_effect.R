#Author: Thong Luong
#Date: December 4th 2025

files = list.files('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/EQTL_effect_byQuantile_Finite_time/', pattern = '_1.csv')

files = gsub('_1.csv','',files)
eqtl = gsub('result_epithelial_','',files)
data = data.frame(matrix(nrow = length(files), ncol = 6))
rownames(data) = eqtl

for (f in 1:length(files)){
  for (i in 1:6){
    t = read.csv(paste0('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/EQTL_effect_byQuantile_Finite_time/',files[f],'_',i,'.csv'))
    data[f,i] = t[2,'Estimate']
  }
}

for (f in 1:length(files)){
  t = read.csv(paste0('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/EQTL_effect_byQuantile_Finite_time/',files[f],'_1.csv'))
  rownames(data)[f] = paste(unique(t$gene),unique(t$snp), sep = '|')
}

colnames(data) = c('Q1','Q2','Q3','Q4','Q5','Q6')
saveRDS(data, '/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/int_eqtl_effect.rds')
