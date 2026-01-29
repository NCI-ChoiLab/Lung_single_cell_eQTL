#======================
#Script used to overlay annotations from each subgroup onto the unfiltered Seurat object and remove unwanted cells before creating the final Seurat object with final annotations.
#======================

library(Seurat)
library(ggplot2)
library(harmony)

#Unfiltered Seurat Object that includes later removed cells
combined <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/final_no_spx.RDS')

combined$cell_types <- NA

# Annotating stromal cells - Reclustered with a resolution of 0.45
stromal <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/UMAPStromres0.45_w_removed')
stromal$cell_types <- NULL
stromal_cellTypeAnnotations <- c('Adventitial Fibroblasts', 'Alveolar Fibroblasts', 'Smooth Muscle', 'Smooth Muscle', 'Subpleural Fibroblasts', 'Peribronchial Fibroblasts', 'Mesothelium', 'Myofibroblasts')
stromal$cell_types <- stromal_cellTypeAnnotations[stromal$seurat_clusters]


# Annotating endothelial cells - Reclustered with a resolution of 0.25
endothelial <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/UMAPendo_w_removed')
endothelial$cell_types <- NULL
endothelial_cellTypeAnnotations <- c('Lymphatic EC Mature', 'EC Arterial', 'EC Venous Pulmonary', 'EC Venous Systemic', 'EC General Capillary', 'EC Aerocyte Capillary', 'Lymphatic EC Proliferating', 'Lymphatic EC Mature')
endothelial$cell_types <- endothelial_cellTypeAnnotations[endothelial$seurat_clusters]



# Annotating epithelial cells - Reclustered with a resolution of 0.60
epithelial <- readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/EpiImmResolutions/UMAPepiapp3.60_updated.rds')
epithelial$cell_types <- NULL
epithelial_cellTypeAnnotations <- c('AT1', 'Multiciliated', 'AT2', 'Club', 'Alveolar Transitional Cells', 'Secretory Transitional Cells', 'Basal', 'Goblet', 'Multiciliated', 'AT2',
                                    'Club', 'Secretory Transitional Cells', 'Goblet', 'Alveolar Transitional Cells', 'Alveolar Transitional Cells', 'Alveolar Transitional Cells', 'AT2', 'Multiciliated', 'Neuroendocrine')
epithelial$cell_types <- epithelial_cellTypeAnnotations[epithelial$seurat_clusters]

saveRDS(epithelial, '/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/epithelial_finalannotations0601.rds')


# Annotating immune cells - Reclustered with a resolution of 1.7  (CD3+ cells, which had been cluster 29 in the immune subgroup res 1.7 were separated out to be further reclustered hence the multiple groups here)

immune <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/filteredimmune1.7.rds')
immune$cell_types <- NULL
immune_cellTypeAnnotations <- c('Alveolar Mph', 'NK Cells', 'NK Cells', 'NK Cells', 'CD8+ T Cells', 'CD4+ T Cells', 'Non-classical Monocytes', 'NK Cells', 'Classical Monocytes', 'CD4+ T Cells', 
                                'NK Cells', 'CD4+ T Cells', 'CD8+ T Cells', 'Alveolar Mph', 'CD8+ T Cells', 'Classical Monocytes', 'Alveolar Mph', 'NK Cells', 'CD4+ T Cells', 'DC2', 
                                'Monocyte-derived Mph', 'Interstitial Mph', 'Mast Cells', 'Monocyte-derived Mph', 'Alveolar Mph', 'Alveolar Mph', 'NK Cells', 'Monocyte-derived Mph', 'Proliferating Macrophages', 'CD4+ T Cells', 
                                'B Cells', 'DC1', 'DC2', 'Migratory DCs', 'Plasma Cells', 'Plasmacytoid DCs', 'Alveolar Mph')
immune$cell_types <- immune_cellTypeAnnotations[immune$seurat_clusters] 


# Annotating immune CD3+ cells
immuneCD3 <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/UMAPCD3')
immuneCD3$cell_types <- NULL
immuneCD3_cellTypeAnnotations <- c('Proliferating T', 'Proliferating NK', 'Proliferating T')
immuneCD3$cell_types <- immuneCD3_cellTypeAnnotations[immuneCD3$seurat_clusters] 


# Annotating removed immune cells to later remove 
filter <- readRDS('/data/Choi_lung/Elelta/eQTLProjects/FinalizedCellAnnotation0412/removedimmunecells')
filter$cell_types <- rep('Remove', nrow(filter))

# Combining all annotations
combined$cell_types[Cells(stromal)] <- stromal$cell_types
combined$cell_types[Cells(endothelial)] <- endothelial$cell_types
combined$cell_types[Cells(epithelial)] <- epithelial$cell_types
combined$cell_types[Cells(immune)] <- immune$cell_types
combined$cell_types[Cells(immuneCD3)] <- immuneCD3$cell_types
combined$cell_types[Cells(filter)] <- filter$cell_types
combined <- subset(combined, subset = cell_types != 'Remove')

reorder <- c('Alveolar Transitional Cells','AT1', 'AT2', 'Club', 'Goblet', 'Basal', 'Multiciliated', 'Neuroendocrine', 'Secretory Transitional Cells',
             'Alveolar Mph', 'Proliferating Macrophages', 'NK Cells', 'Proliferating NK', 'CD8+ T Cells', 'CD4+ T Cells', 'Proliferating T', 'Interstitial Mph', 'Monocyte-derived Mph', 'Classical Monocytes', 'Non-classical Monocytes', 'DC1', 'DC2', 'Migratory DCs', 'Mast Cells', 'B Cells', 'Plasma Cells', 'Plasmacytoid DCs', 
             'Lymphatic EC Mature', 'Lymphatic EC Proliferating', 'EC Arterial', 'EC Venous Pulmonary', 'EC Venous Systemic', 'EC Aerocyte Capillary', 'EC General Capillary', 
             'Adventitial Fibroblasts', 'Smooth Muscle', 'Subpleural Fibroblasts', 'Peribronchial Fibroblasts', 'Alveolar Fibroblasts', 'Mesothelium', 'Myofibroblasts')
combined$cell_types <- factor(combined$cell_types, levels = reorder)


# Saving the annotated Seurat object
saveRDS(combined, file = "/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/finalcombined3.rds")



