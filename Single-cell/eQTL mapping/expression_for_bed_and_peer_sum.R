#Author: Thong Luong
#Date: Dec 6th 2024

setwd('/data/Choi_lung/TTL/tensor/Phenotype/Sum_Final')

library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(Matrix.utils)
library(data.table)
library(genio)
library(stringr)
library(edgeR)
library(scran)

#loading GTF file so we can match gene name to ensembl_id
gff = as.data.frame(rtracklayer::import('/data/Choi_lung/TTL/refdata-gex-GRCh38-2020-A/genes/genes.gtf'))

#load Seurat object
total_currently = readRDS('/data/Choi_lung/TTL/tensor/Cell_types/finalcombined.rds')
total_currently[['RNA']] = JoinLayers(total_currently[['RNA']])
CT = unique(total_currently$cell_types) 

#match gene symbol with ENSEMBLID, sce is your single-cell-experiment object that you just created, gff is your gtf
symbol2ensembl <- function(sce, gff){
  #subset gff so there's only genes
  gff = subset(gff, subset = type == 'gene')
  sceRowData <- as.data.frame(list(gene_name=rownames(sce))) %>% 
    left_join(gff, by="gene_name") %>%
    dplyr::select(gene_id=gene_id, symbol=gene_name, gene_chr=seqnames, 
                  gene_start=start, gene_end=end)
  
  dupEntries <- unique(sceRowData$symbol[which(duplicated(sceRowData$symbol))])
  message("# symbols mapping to multiple Ensembl IDs: ", length(dupEntries))
  sceRowData <- sceRowData[!sceRowData$symbol %in% dupEntries, ]
  sceRowData <- na.omit(sceRowData)
  
  sce <- sce[sceRowData$symbol, ]
  rowData(sce) <- sceRowData
  rownames(sce) <- rowData(sce)$gene_id
  return(sce)
}


#function to add tss to expression data
add_tss = function(m){
  #load tss
  pos = read.csv('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Annotation_final/tss.tsv', sep = '\t')
  pos2 = pos[,-1]
  rownames(pos2) = pos[,1]
  merged.df <- merge(pos2, m ,all=F,by='row.names')
  final = merged.df[,-1]
  colnames(final)[1] = '#chr'
  return(final)
}

aggr_sum = function(e,s,ag,r,b,p){
  #telling what to aggregate by
  groups = colData(e)[, c('Sample')]
  
  #samples are rows and genes are columns
  #sum aggregation
  aggr_sum = aggregate.Matrix(t(counts(e)), groupings = groups, fun = 'sum')
  
  #transposing so samples are columns and genes are row
  expression_sum = t(aggr_sum)
  
  #normalize
  expression_sum = NormalizeData(expression_sum)
  
  new_dataFrame <- apply(expression_sum, 1, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  new_dataFrame = t(new_dataFrame)
  
  final = add_tss(new_dataFrame)
  
  print(table(final$'#chr'))
  
  #finding location of weird chromosome, you can use table(AT2_mean$'#chr') to see which chr are present, get rid of any chromosome that isn't 1-22 or X. To do that first find location of those chr, then subset to not include. 
  chrM = which(final$`#chr` == 'chrM')
  GL = which(final$`#chr` == 'GL000219.1')
  
  expression_sum = expression_sum[-c(chrM,GL),]
  
  write.table(expression_sum, file = ag, sep='\t', quote=F, row.names=T, col.names = T, eol = '\n')
  
  #subset to remove unwanted chromosome. Notice the - sign, this is telling R to get rid of it. 
  final = final[-c(chrM,GL),]
  
  saveRDS(final, file = r)
  
  #writing the bed file
  write.table(final, file = b, sep='\t', quote=F, row.names=FALSE, col.names = T, eol = '\n')
  
  #making expression file for PEER
  make_expression_table = function(m){
    ta = m[,-c(1,2,3)]
    t2 = data.frame(t(ta[-1]))
    colnames(t2) = ta[,1]
    return(t2)
  }
  
  peer = make_expression_table(final)
  
  #saving expression file for PEER
  write.table(peer, file = p, sep = ',', quote = F, row.names = F, col.names = F)
  
}

load_and_make_expression = function(x){
  obj = x
  counts = obj@assays$RNA$counts
  meta = obj@meta.data
  
  #only keeping genes that is expressed in at least 10% of all cells
  counts = counts[rowSums(counts >0) >= (ncol(counts)/10),]
  #keeping genes with mean counts of greater than 0.1
  counts = counts[rowMeans(counts) > 0.1,]
  
  sce2 = SingleCellExperiment(assays = list(counts = counts),
                             colData = meta)
  
  sce2$Sample = factor(sce2$Sample)
  sce2$cell_types = factor(sce2$cell_types)

  #initial number of genes
  n_genes = nrow(sce2)
  #matching 
  sce_clone = symbol2ensembl(sce2, gff)
  message("Genes removed (not in GFF or multiple Ensembl IDs): ", 
          n_genes - nrow(sce_clone))
  
  for (i in CT){
    cell = sce_clone[,sce_clone$cell_types == i]
    sample = read.table(paste0('/data/Choi_lung/TTL/tensor/Genotype/',i,'_samples.txt'), sep = '\t')
    cell = cell[,cell$Sample %in% sample$V1]
    
    aggr_sum(cell,sample,paste0(i,'_agg.tsv'), paste0(i,'.rds'),paste0(i,".bed"),paste0(i,'_peer.csv'))
    
  }
}

load_and_make_expression(total_currently)






