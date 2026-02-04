#Author: Thong Luong
#Date: Dec 16th 2024

setwd('/data/Choi_lung/TTL/Colocalization/jianxins/')

files = list.files()

asn = read.table(files[2], sep = '\t', header = T)

meta = asn

lead_snps = c('rs17038564','rs682888','rs137884934','rs2293607','rs55779747',
              'rs117715768','rs1373058','rs2736100','rs9380190','rs2760995',
              'rs9367106','rs531557','rs6937083','rs4268071','rs72658409',
              'rs11196089','rs10901793','rs174559','rs55768116','rs7962469',
              'rs1200399','rs71467682','rs764014','rs59956089','rs116863980')

meta_twas = c('rs1130866','rs2320614','rs766826',
              'rs34638657','rs638868')

final_snps = c(lead_snps,meta_twas)

liftover_function = function(m,l){
  liftover_snps = read.table(l)
  liftover_snps = liftover_snps[,c(1,3,4)]
  colnames(liftover_snps) = c('CHR','POS','SNP')
  
  meta = subset(m, SNP %in% liftover_snps$SNP)
  meta$CHR = NULL
  meta$BP = NULL
  
  meta = merge(meta, liftover_snps, by = 'SNP')
  meta = meta[,c('CHR','POS','SNP','A1','A2','N','P','P.R.','OR','OR.R.','Q','I')]
  return(meta)
}

asn = liftover_function(asn,'hg38_asn_snps.bed')

write.table(asn, './hg38_asn/asian.lung.meta.txt', sep = '\t', col.names = T, row.names = F, quote = F)
