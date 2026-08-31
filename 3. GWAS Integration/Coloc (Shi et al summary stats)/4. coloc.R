#Author: Thong Luong
#Date: Jan 15th 2025

library(coloc)
library(foreach)
library(doParallel)
library(parallel)
library(stringr)

setwd('/data/Choi_lung/TTL/Colocalization/jianxins/hg38_asn/sig')

gwas = read.table('/data/Choi_lung/TTL/Colocalization/jianxins/hg38_asn/asian.lung.meta.txt', header = T)


ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')
ss = c(96,92,127,128,124,126,128,63,127,126,127,128,57,123,58,128,82,
       120,117,126,111,128,82,101,48,125,128,128,123,82,127,114,89)
loci = list.files('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Shi_sig/')
not_sig = c('4p13','4q32.1','4q32.2','7q31.33','9p21.3')
loci = setdiff(loci,not_sig)
  
for (k in loci){
  dir.create(paste0(k),showWarnings = F)
  for (j in 1:length(ct)){
    eqtl = read.table(paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Shi_sig/',k,'/',ct[j],'_',k,'_100kb.tsv'), sep = '\t', header = T)
    eqtl$rsid = paste(eqtl$rsid,eqtl$alt,sep = '_')
    eqtl$varbeta_eqtl = (eqtl$slope_se)^2
    eqtl = eqtl[!(is.na(eqtl$varbeta_eqtl)),]
    
    gwas = read.table(paste0('/data/Choi_lung/TTL/Colocalization/jianxins/hg38_asn/',k,'_100kb.tsv'), sep = '\t', header = T)
    gwas$z = sign(log(gwas$OR.R.))*abs(qnorm((gwas$P)/2))
    gwas$se = (log(gwas$OR.R.))/(gwas$z)
    gwas$varbeta_gwas = (gwas$se)^2
    gwas = gwas[!(is.na(gwas$varbeta_gwas)),]
    gwas$rsid = paste(gwas$SNP,gwas$A1,sep = '_')
    gwas = subset(gwas, rsid %in% eqtl$rsid)
    
    #merging dataset
    merged_dataset <- merge(eqtl, gwas, by="rsid", all=FALSE)
    
    ##filtering snps to only include snps in the eQTL list and vice versa
    filt_gwas=gwas[gwas$rsid %in% eqtl$rsid,]
    filt_eqtl=merged_dataset[merged_dataset$rsid %in% filt_gwas$rsid,]
    
    gene_list=as.character(unique(filt_eqtl$phenotype_name))
    eqtl_genes=list()
    eqtl_genes=foreach(i=gene_list) %dopar%
      droplevels(filt_eqtl[filt_eqtl$phenotype_name==i,]) # Dropping unused levels makes a HUGE difference in time and final list size
    names(eqtl_genes)=gene_list
    
    ##note in the drop levels, it needs to be col 6 since that is where the rsIDs are for the eqtl data, otherwise you get an empty snp gene list.
    gwas_genes=list()
    gwas_genes=foreach(i=gene_list) %dopar%
      droplevels(filt_gwas[na.omit(match(as.character(eqtl_genes[[i]][,1]),as.character(filt_gwas$rsid))),]) 
    names(gwas_genes)=gene_list
    
    SNP_coloc=list()
    SNP_coloc=foreach(i=(1:length(eqtl_genes)),.packages='coloc') %dopar%
      coloc.abf(dataset1 = list(beta = eqtl_genes[[i]]$slope, varbeta = eqtl_genes[[i]]$varbeta_eqtl, snp=eqtl_genes[[i]]$rsid,pvalues=eqtl_genes[[i]]$pval_nominal,position=eqtl_genes[[i]]$variant_pos, N = ss[j],type="quant",MAF=eqtl_genes[[i]]$af),
                dataset2 = list(beta = log(gwas_genes[[i]]$OR.R.), varbeta = gwas_genes[[i]]$varbeta_gwas, snp=gwas_genes[[i]]$rsid,pvalues=gwas_genes[[i]]$P,position=gwas_genes[[i]]$POS,N=42315,type="cc",s=0.28))
    names(SNP_coloc)=gene_list
    
    output_genes=list()
    output_genes=foreach(i=gene_list) %dopar%
      capture.output(SNP_coloc[[i]][["summary"]])
    names(output_genes)=gene_list
    
    saveRDS(SNP_coloc, paste0(k,'/',ct[j],'_SNP_coloc.rds'))
    saveRDS(output_genes, paste0(k,'/',ct[j],'_output_genes.rds'))
  }
}






