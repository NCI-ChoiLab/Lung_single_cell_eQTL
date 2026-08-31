#Author: Thong Luong
#Date: Jan 17th 2025

library(coloc)
library(foreach)
library(doParallel)
library(parallel)
library(stringr)

setwd('/data/Choi_lung/TTL/Colocalization/Byun_gwas/hg38_total/coloc/')

ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')
ss = c(96,92,127,128,124,126,128,63,127,126,127,128,57,123,58,128,82,
       120,117,126,111,128,82,101,48,125,128,128,123,82,127,114,89)

loci = c('1p31.1','2p14','2p16.1','2q34','3p22.1','3p25.3','3q28','5p15.33_1','5p15.33_2','6p21.32_1','6p21.32_2',
          '6p21.33','6p22.1_1','6p22.1_2','6p22.1_3','6p25.3_1','6p25.3_2','6p27','6q21','6q22.1','8p12','8p21.1','8p21.2',
          '9p21.3_1','9p21.3_2','10q24.3','10q24.33','10q25.2_1','10q25.2_2','11q22.3','11q23.3_1','11q23.3_2','12p13.33',
          '13q13.1_1','13q13.1_2','15q21.1_1','15q21.1_2','15q25.1_1','15q25.1_2','19q13.2','20q13.33','22q12.1')

not_sig = c('6p25.3_1','6p25.3_2')
loci = setdiff(loci,not_sig)

for (k in loci){
  dir.create(paste0(k),showWarnings = F)
  for (j in 1:length(ct)){
    eqtl = read.table(paste0('/data/Choi_lung/TTL/files_for_JC_2501/sliced_Byun_sig/',k,'/',ct[j],'_',k,'_100kb.tsv'), sep = '\t', header = T)
    eqtl$rsid = paste(eqtl$rsid,eqtl$alt,sep = '_')
    eqtl$varbeta_eqtl = (eqtl$slope_se)^2
    eqtl = eqtl[!(is.na(eqtl$varbeta_eqtl)),]
    
    gwas = read.table(paste0('/data/Choi_lung/TTL/Colocalization/Byun_gwas/hg38_total/',k,'_100kb.tsv'), sep = '\t', header = T)
    gwas$z = sign(log(gwas$odds_ratio))*abs(qnorm((gwas$p_value)/2))
    gwas$se = (log(gwas$odds_ratio))/(gwas$z)
    gwas$varbeta_gwas = (gwas$se)^2
    gwas = gwas[!(is.na(gwas$varbeta_gwas)),]
    colnames(gwas)[1] = 'SNP'
    gwas$rsid = paste(gwas$SNP,gwas$effect_allele,sep = '_')
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
                dataset2 = list(beta = log(gwas_genes[[i]]$odds_ratio), varbeta = gwas_genes[[i]]$varbeta_gwas, snp=gwas_genes[[i]]$rsid,pvalues=gwas_genes[[i]]$p_value,position=gwas_genes[[i]]$POS,N=70156,type="cc",s=0.51))
    names(SNP_coloc)=gene_list
    
    output_genes=list()
    output_genes=foreach(i=gene_list) %dopar%
      capture.output(SNP_coloc[[i]][["summary"]])
    names(output_genes)=gene_list
    
    saveRDS(SNP_coloc, paste0(k,'/',ct[j],'_SNP_coloc.rds'))
    saveRDS(output_genes, paste0(k,'/',ct[j],'_output_genes.rds'))
  }
}






