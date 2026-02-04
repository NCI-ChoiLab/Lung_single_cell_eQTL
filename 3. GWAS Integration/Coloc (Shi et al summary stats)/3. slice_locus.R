#Author: Thong Luong
#Date: Dec 16th 2024

setwd('/data/Choi_lung/TTL/Colocalization/jianxins/hg38_asn/')


gwas = read.table('asian.lung.meta.txt', sep = '\t', header = T)


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

for (i in 1:length(locus)){
  gwas_sub = subset(gwas, CHR == chr[i])
  gwas_sub = subset(gwas_sub, POS >= (pos[i] - 100000) & POS <= (pos[i] + 100000))
  write.table(gwas_sub, file = paste0(locus[i],'_100kb.tsv'), sep = "\t",
              quote = FALSE, row.names = FALSE)
}
