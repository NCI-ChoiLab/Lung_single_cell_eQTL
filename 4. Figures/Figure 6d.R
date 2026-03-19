#Author: Thong Luong
#Date: July 29th 2025
#Modified: Sept 19th

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

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/ROS1_DEG/TCF7L2')

sc = readRDS('/data/Choi_lung/TTL/Clinical_eQTL/total_with_preinfo.rds')

at2 = subset(sc, cell_types == 'AT2')

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
  
  sce2$Sample = factor(sce2$Sample)
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


at2_matrix = get_matrix(at2)


get_deg = function(mat){
  coldata = as.data.frame(NormalizeData(mat))
  coldata = as.data.frame(t(coldata['TCF7L2',]))
  coldata$q = coldata %>% mutate(q = ntile(TCF7L2,4)) %>% pull(q)
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


at2_deg = get_deg(at2_matrix)

#saveRDS(at2_deg,'at2_deg.rds')

save_deg = function(lung.expr, f){
  up.genes <- lung.expr[lung.expr$log2FoldChange > .5 & lung.expr$padj < 0.05, ]
  dn.genes <- lung.expr[lung.expr$log2FoldChange < -.5 & lung.expr$padj < 0.05,]
  deg = rbind(up.genes,lung.expr['ROS1',])
  deg = rbind(deg,dn.genes)
  write.table(deg, f, sep = '\t', row.names = T, col.names = T, quote = F)
}


#save_deg(at2_deg, 'pseudobulk_at2_deg.tsv')



at2_deg = readRDS('at2_deg.rds')


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


at2_path = get_pathway(at2_deg)

library(enrichplot)
library(DOSE)

edo <- pairwise_termsim(at2_path)

#emapplot(edo)
#edo2 <- gseDO(at2_path)

t = gseaplot2(edo, geneSetID = c('WP2857','WP2853','WP428','WP3651','WP4336'), pvalue_table = F, base_size = 5)
plot_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_6/'
ggsave(filename=file.path(plot_dir, 'TCFL72_NES.pdf'), plot = t, width = 80, height = 50, units = "mm", dpi = 450, device = 'pdf')

data = as.data.frame(edo)
data = subset(data, ID %in% c('WP2857','WP2853','WP428','WP3651','WP4336'))

data$log_FDR = -log10(data$qvalue)
data = data[order(data$log_FDR),]
data$Description = factor(data$Description, levels = c('ncRNAs in Wnt signaling in hepatocellular carcinoma',
                                                       'Pathways affected in adenoid cystic carcinoma','Wnt signaling',
                                                       'Mesodermal commitment pathway','Endoderm differentiation'
                                                       ))
t2 = ggplot()+
  geom_point(data, mapping = aes(x= NES, y= Description, color = Description, size=log_FDR)) +
  scale_colour_manual(guide = 'none', values = c('Endoderm differentiation' = "indianred3", 'Mesodermal commitment pathway' = "darkolivegreen3", 
                                 'ncRNAs in Wnt signaling in hepatocellular carcinoma' = "seagreen3", 'Pathways affected in adenoid cystic carcinoma' = "turquoise4",
                                 'Wnt signaling' = 'mediumorchid3'))  +
  scale_size(range = c(.5,3),limits=c(3,8),breaks=c(3,4,7),labels=c(">=3",'>=4','>=7'),guide="legend", name = '-Log10(FDR)') + 
  labs(y = "Gene Set", x = "Normalized Enrichment Score", fill = "NES") +
  scale_x_continuous(limits = c(1.9, 2.2)) +
  scale_y_discrete(labels= c('ncRNAs in Wnt signaling \n in hepatocellular carcinoma',
                            'Pathways affected in \n adenoid cystic carcinoma','Wnt signaling',
                            'Mesodermal \n commitment pathway','Endoderm differentiation')) +
  ggtitle('') +
  theme(text = element_text(family = 'Arial'),
        axis.text=element_text(size=5.5),
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(3, "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        panel.background = element_rect(fill = 'white', colour = 'black', linetype = 'solid'),
        legend.key = element_rect(color = NA, fill = NA))


ggsave(filename=file.path(plot_dir, 'TCFL72_bubble.png'), plot = t2, width = 80, height = 60, units = "mm", dpi = 450, device = 'png')







