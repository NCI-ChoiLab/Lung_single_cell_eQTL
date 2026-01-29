#Script used to generate a h5ad file to run with nsForest with R

library(Matrix)
library(future)
library(Seurat)
library(reticulate)

use_condaenv('leidenpackage')
library(zellkonverter)

setwd('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/')
plot_save_dir <- '/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/'


combined <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/finalcombined3.rds')

combined[['RNA']] <- JoinLayers(combined[['RNA']])
combined[['predicted.ann_level_1.score']] <- NULL
combined[['predicted.ann_level_2.score']] <- NULL
combined[['predicted.ann_level_3.score']] <- NULL
combined[['predicted.ann_level_4.score']] <- NULL
combined[['predicted.ann_level_5.score']] <- NULL
combined[['predicted.ann_finest_level.score']] <- NULL
sce_obj <- as.SingleCellExperiment(combined, assay = c("RNA"))

writeH5AD(sce_obj, file.path(plot_save_dir, "wholegroupjuly.h5ad"), X_name = 'counts')
