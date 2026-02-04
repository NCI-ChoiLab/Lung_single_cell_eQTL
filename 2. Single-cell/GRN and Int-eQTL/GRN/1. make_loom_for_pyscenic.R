#Author: Thong Luong
#Date: June 13th 2025

library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/GRN/')

#load Seurat object
sc = readRDS('/data/Choi_lung/TTL/Clinical_eQTL/total_with_preinfo.rds')
sc[['RNA']] = JoinLayers(sc[['RNA']])
Idents(sc) = 'cell_types'

epithelial_cells = c('Alveolar Transitional Cells', 'AT1', 'AT2', 'Basal','Club','Goblet',
                     'Multiciliated','Neuroendocrine','Secretory Transitional Cells')

sc = subset(sc, cell_types %in% epithelial_cells)

counts =  sc@assays$RNA$counts
cellInfo = sc@meta.data
cellInfo$CellID = rownames(cellInfo)
cellInfo = cellInfo[,c("nFeature_RNA","nCount_RNA",'cell_types')]
colnames(cellInfo)[1:2] = c('nGene','nUMI')

#only keeping genes that is expressed in at least 10% of all cells
counts = counts[rowSums(counts >0) >= (ncol(counts)/10),]
#keeping genes with mean counts of greater than 0.1
counts = counts[rowMeans(counts) > 0.1,]

add_cell_annotation <- function(loom, cellAnnotation){
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

loom <- build_loom("databases/Epi_filtered.loom", dgem=counts)
loom <- add_embedding(loom, embedding = sc@reductions$umap.harmony@cell.embeddings, name = 'umap')
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)


