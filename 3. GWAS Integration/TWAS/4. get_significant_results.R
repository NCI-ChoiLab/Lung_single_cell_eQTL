#Author: BL
#Modified: Thong Luong

setwd("/data/Choi_lung/TTL/FUSION/")

ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

rst.list <- list()
for (celltype in ct) {
  rst_path <- paste0("/data/Choi_lung/TTL/FUSION/", celltype,'/compute_weight_full/')
  dat.files <- list.files(rst_path, pattern = ".dat")
  if (length(dat.files) == 0) {
    rst.list[[celltype]] <- NULL
  }else{
    dat.list <- list()
    for (i in 1:length(dat.files)) {
      rst <- read.table(paste0(rst_path, "/",dat.files[i]), sep = "\t", header = TRUE)
      dat.list[[i]] <- rst
      
    }
    dat.rst <- Reduce(rbind, dat.list)
    dat.rst <- dat.rst[!is.na(dat.rst$TWAS.P),]
    dat.rst <- dat.rst[!is.na(dat.rst$BEST.GWAS.ID),]
    dat.rst$FDR <- p.adjust(dat.rst$TWAS.P, method = "BH" )
    dat.rst$Celltype <- celltype
    rst.list[[celltype]] <- dat.rst
  }
}
TWAS_sum <- sapply(rst.list, function(x) sum(x$FDR < 0.05))

full_cov = rst.list

gene_loc = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/gene_loc.rds')
gene_loc = gene_loc[,c("phenotype_id","phenotype_name")]
colnames(gene_loc) = c('ID','ID.NAME')

geno = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/snps.rds')
geno$full = paste(geno$V1,geno$V4,geno$V6,geno$V5,sep = '_')
geno = geno[,c("full",'V2')]
colnames(geno)[2] = 'RSID'

get_sig = function(l){
  df = data.frame()
  for (i in 1:length(l)){
    cell = as.data.frame(l[i]) 
    cell = cell[,c(3:22)]
    colnames(cell) = c('ID','CHR','P0','P1','HSQ','BEST.GWAS.ID','GWAS.Z','EQTL.ID','EQTL.R2','EQTL.Z',
                       'EQTL.GWAS.Z','NSNP','NWGT','MODEL','MODELCV.R2','MODELCV.PV','TWAS.Z','TWAS.P','FDR','CELLTYPE')
    gene = subset(gene_loc, ID %in% cell$ID)
    cell = merge(cell,gene, by = 'ID')
    gwas_geno = subset(geno, full %in% cell$BEST.GWAS.ID)
    colnames(gwas_geno) = c('BEST.GWAS.ID','BEST.GWAS.RSID')
    
    eqtl_geno = subset(geno, full %in% cell$EQTL.ID)
    colnames(eqtl_geno) = c('EQTL.ID','EQTL.RSID')
    
    cell = merge(cell,gwas_geno, by = 'BEST.GWAS.ID')
    cell = merge(cell,eqtl_geno, by = 'EQTL.ID')
    
    cell = cell[,c('ID','ID.NAME','CHR','P0','P1','HSQ','BEST.GWAS.ID','BEST.GWAS.RSID','GWAS.Z','EQTL.ID','EQTL.RSID','EQTL.R2','EQTL.Z',
                   'EQTL.GWAS.Z','NSNP','NWGT','MODEL','MODELCV.R2','MODELCV.PV','TWAS.Z','TWAS.P','FDR','CELLTYPE')]
    cell = subset(cell, FDR < 0.05)
    df = rbind(df,cell)
  }
  return(df)
}

one_kg = get_sig(rst.list)

full_cov_sig = get_sig(full_cov)

write.table(full_cov_sig,'./Covariates_comp/Our_TWAS.txt',quote = F, sep = '\t', col.names = T, row.names = F)



