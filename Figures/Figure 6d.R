#Author: Thong Luong
#Date: July 29th 2025

library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(Matrix)
library(Matrix.utils)
library(data.table)
library(textshape)
library(limma)
library(KEGGREST)
library(stringr)
library(org.Hs.eg.db)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/ROS1_DEG/')

sc = readRDS('/data/Choi_lung/TTL/Clinical_eQTL/total_with_preinfo.rds')
rel = subset(sc, cell_types %in% c('Alveolar Transitional Cells', 'AT2', 'Secretory Transitional Cells'))
alv_trans = subset(rel, cell_types == 'Alveolar Transitional Cells')
at2 = subset(rel, cell_types == 'AT2')
alv_cells = subset(rel, cell_types != 'Secretory Transitional Cells')
sec_trans = subset(rel, cell_types == 'Secretory Transitional Cells')

get_matrix = function(s){
  s[['RNA']] = JoinLayers(s[['RNA']])
  counts = s@assays$RNA$counts
  meta = s@meta.data
  
  #only keeping genes that is expressed in at least 10% of all cells
  counts = counts[rowSums(counts >0) >= (ncol(counts)/10),]
  #keeping genes with mean counts of greater than 0.1
  
  counts = counts[rowMeans(counts) > 0.1,]
  sce2 = SingleCellExperiment(assays = list(counts = counts),
                              colData = meta)
  
  sce2$orig.ident = factor(sce2$Sample)
  e = sce2
  #telling what to aggregate by
  groups = colData(e)[, c('Sample')]
  
  #samples are rows and genes are columns
  #sum aggregation
  aggr_sum = aggregate.Matrix(t(counts(e)), groupings = groups, fun = 'sum')
  
  
  #transposing so samples are columns and genes are row
  expression_sum = t(aggr_sum)
  expression_sum = as.data.frame(expression_sum)
  return(expression_sum)
}

rel_matrix = get_matrix(rel)
at2_matrix = get_matrix(at2)
alv_trans_matrix = get_matrix(alv_trans)
alv_cells_matrix = get_matrix(alv_cells)
sec_trans_matrix = get_matrix(sec_trans)

get_deg = function(mat){
  coldata = as.data.frame(NormalizeData(mat))
  coldata = as.data.frame(t(coldata['ROS1',]))
  coldata$q = coldata %>% mutate(q = ntile(ROS1,4)) %>% pull(q)
  coldata = subset(coldata , q %in% c('1','4'))
  coldata$q = factor(coldata$q, levels = c('4','1'))
  
  cts = mat[,rownames(coldata)]
  dds = DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                       colData = coldata,
                                       design = ~q)
  dds_finals = DESeq2::DESeq(dds)
  dds_results = DESeq2::results(dds_finals)
  gene_deg = as.data.frame(dds_results)
  gene_deg = gene_deg[order(gene_deg$log2FoldChange, decreasing = T),]
  return(gene_deg)
}

rel_deg = get_deg(rel_matrix)
alv_deg = get_deg(alv_trans_matrix)
at2_deg = get_deg(at2_matrix)
alv_cells_deg = get_deg(alv_cells_matrix)
sec_deg = get_deg(sec_trans_matrix)


get_deg2 = function(mat){
  coldata = as.data.frame(NormalizeData(mat))
  coldata = as.data.frame(t(coldata['ROS1',]))
  coldata$q = coldata %>% mutate(q = ntile(ROS1,4)) %>% pull(q)
  coldata = subset(coldata , q %in% c('1','4'))
  coldata$q = factor(coldata$q, levels = c('1','4'))
  
  cts = mat[,rownames(coldata)]
  dds = DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                       colData = coldata,
                                       design = ~q)
  dds_finals = DESeq2::DESeq(dds)
  dds_results = DESeq2::results(dds_finals)
  gene_deg = as.data.frame(dds_results)
  gene_deg = gene_deg[order(gene_deg$log2FoldChange, decreasing = T),]
  return(gene_deg)
}



#saveRDS(rel_deg, 'rel_deg.rds')
#saveRDS(alv_deg, 'alv_deg.rds')
#saveRDS(at2_deg,'at2_deg.rds')
#saveRDS(alv_cells_deg, 'alv_cells_deg.rds')
#saveRDS(sec_deg,'sec_deg.rds')


save_deg = function(lung.expr, f){
  up.genes <- lung.expr[lung.expr$log2FoldChange > .5 & lung.expr$padj < 0.05, ]
  dn.genes <- lung.expr[lung.expr$log2FoldChange < -.5 & lung.expr$padj < 0.05,]
  deg = rbind(up.genes,lung.expr['ROS1',])
  deg = rbind(deg,dn.genes)
  write.table(deg, f, sep = '\t', row.names = T, col.names = T, quote = F)
}

save_deg(rel_deg, 'pseudobulk_rel_deg.tsv')
save_deg(alv_deg, 'pseudobulk_alv_deg.tsv')
save_deg(at2_deg, 'pseudobulk_at2_deg.tsv')
save_deg(alv_cells_deg, 'pseudobulk_alv_cells_deg.tsv')
save_deg(sec_deg, 'pseudobulk_sec_deg.tsv')

rel_deg = readRDS('rel_deg.rds')
alv_deg = readRDS('alv_deg.rds')
at2_deg = readRDS('at2_deg.rds')
alv_cells_deg = readRDS('alv_cells_deg.rds')
sec_deg = readRDS('sec_deg.rds')

get_pathway = function(deg){
  lung.expr = deg
  bkgd.genes <- rownames(lung.expr)
  bkgd.genes.entrez <- clusterProfiler::bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  
  lung.expr$GeneID = rownames(lung.expr)
  
  lung.expr$fcsign <- sign(lung.expr$log2FoldChange)
  lung.expr$logfdr <- -log10(lung.expr$pvalue)
  lung.expr$sig <- lung.expr$logfdr/lung.expr$fcsign
  sig.lung.expr.entrez<-merge(lung.expr, bkgd.genes.entrez, by.x = "GeneID", by.y = "SYMBOL")
  gsea.sig.lung.expr <- sig.lung.expr.entrez[,10]
  names(gsea.sig.lung.expr) <- as.character(sig.lung.expr.entrez[,11])
  gsea.sig.lung.expr <- sort(gsea.sig.lung.expr,decreasing = TRUE)
  
  gwp.sig.lung.expr <- clusterProfiler::gseWP(
    gsea.sig.lung.expr,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05, #p.adjust cutoff
    organism = "Homo sapiens"
  )
  return(gwp.sig.lung.expr)
}

rel_path = get_pathway(rel_deg)
alv_path = get_pathway(alv_deg)
at2_path = get_pathway(at2_deg)
alv_cells_path = get_pathway(alv_cells_deg)
sec_path = get_pathway(sec_deg)

library(enrichplot)
library(DOSE)

edo <- pairwise_termsim(at2_path)
emapplot(edo)


edo2 <- gseDO(at2_path)

t = gseaplot2(edo, geneSetID = c('WP111','WP623','WP477','WP4324','WP4396'), pvalue_table = F, base_size = 5)
plot_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_6/'
ggsave(filename=file.path(plot_dir, 'ROS1_NES.pdf'), plot = t, width = 80, height = 50, units = "mm", dpi = 450, device = 'pdf')