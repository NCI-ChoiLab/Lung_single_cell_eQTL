#Author: BL
#Modified by TL

library(Seurat)
library(dplyr)
library(ggplot2)
library(future)
library(stringr)
library(stringi)
library(forcats)
library(Matrix)

setwd("/data/Choi_lung/TTL/FUSION/")
# prepare splicing QTL analysis
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/TTL/refdata-gex-GRCh38-2020-A/genes/genes.gtf'))
table(gff$type)
transcript_gff <- subset(gff, type == "gene")
rm(gff)
pos <- data.frame(X = transcript_gff$gene_id,
                  chr = transcript_gff$seqnames, 
                  start = transcript_gff$start, 
                  end = transcript_gff$end)
rownames(pos) <- pos$X

ct = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

for(i in 1:length(ct)){
  wgt_path <- paste0("/data/Choi_lung/TTL/FUSION/", ct[i], "/compute_weight_full/WEIGHTS/")
  wgt.files <- list.files(wgt_path)
  
  transcripts <- str_split_fixed(wgt.files, "\\.", 3)[,1]
  pos_sub <- pos[transcripts,]
  df <- data.frame(WGT = wgt.files,
                   ID = transcripts,
                   CHR = gsub("chr","",pos_sub$chr),
                   P0 = pos_sub$start,
                   P1 = pos_sub$end)
  write.table(df, file = paste0("/data/Choi_lung/TTL/FUSION/", ct[i], "/compute_weight_full/",  ct[i],".pos"), sep = " ",
              row.names = FALSE, quote = FALSE)
}


gwas = read.table('/data/Choi_lung/TTL/Colocalization/jianxins/hg38_asn/asian.lung.meta.txt', sep = '\t', header = T)
gwas$Z <- abs(qnorm(gwas$P/2))*sign(log(gwas$OR))
gwas$N <- NULL
gwas$CHR <- NULL
gwas$POS <- NULL
gwas$P.R. <- NULL
gwas$OR <- NULL
gwas$OR.R. <- NULL
gwas$Q <- NULL
gwas$I <- NULL
head(gwas)
write.table(gwas, file = "Shi_EAS_sumstats.txt", sep = " ",
            row.names = FALSE, quote = FALSE)
