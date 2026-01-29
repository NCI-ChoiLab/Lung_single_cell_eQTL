#======================
# Script used to re-cluster cells within each cell category 
#======================

# set output directory
plot_save_dir <- '/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/'

library(Matrix)
library(future)
library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(harmony)
library(data.table)
library(ggplot2)
library(Azimuth)
library(reticulate)
use_condaenv('leidenpackage')

# size*1024^2, 3000*1024^2
options(future.globals.maxSize= 3145728000)

# loading the object
run = readRDS("/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/QC_and_clustering/res0.1_UMAP.RDS")

# getting rid of unnecessary assays
run@assays$prediction.score.ann_level_1 = NULL
run@assays$prediction.score.ann_level_2 = NULL
run@assays$prediction.score.ann_level_3 = NULL
run@assays$prediction.score.ann_level_4 = NULL
run@assays$prediction.score.ann_level_5 = NULL
run@assays$prediction.score.ann_finest_level = NULL

# copying over old cluster information 
run$old_harmony_clusters = run$harmony_clusters
run$harmony_clusters = NULL
run$old_seurat_clusters = run$seurat_clusters
run$seurat_clusters = NULL

# changing default assay to 'RNA' for SCTransform
DefaultAssay(run) = 'RNA'

# making a copy of old reduction information
run[['old_umap.harmony']] = run[['umap.harmony']]

# Subgroups (i.e., epithelial, immune, endothelial, and stromal cell categories) were originally grouped by their canonical surface markers. 

# subset based on main groups
imm = subset(run, idents = c(1,2,14))
epi = subset(run, idents = c(3,4,5,6,8,11))
endo = subset(run, idents = c(7,9,10))
stroma = subset(run, idents = c(12,13,15))


#=========================================================================================================================================================
# Iterative re-clustering within each cell category

# In some subgroups, a small fraction of cells or an entire cluster wrongly maps to a subgroup or cell type due to doublets or mixed subpopulations. 
# Azimuth predictions and visual inspection of UMAP and surface-marker expression can reveal such artifacts.
#=========================================================================================================================================================

#=========================================================================================================================================================
# Epithelial Cell Category
#=========================================================================================================================================================

obj = SCTransform(epi, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)

obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.60, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path <- file.path(plot_save_dir, 'UMAPepiapp3.60')
saveRDS(obj, file = file_path)

#=======
# Generate table and dim plot to view Azimuth fine-level annotations and UMAP layout
#=======
write.table(table(obj$seurat_clusters, obj$predicted.ann_finest_level), col.names = F, file.path(plot_save_dir, 'AzimuthPredictions_Epi_wo_col.csv'))

t2 = DimPlot(obj, group.by='predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(obj, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "Epithelial_AzimuthAnnotationsDimplot.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")

#======
# Generate feature plots to visualize separation of previously obscured cell types
#=======
feature_plot <- FeaturePlot(epiapp3.40, features = "MKI67", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MKI67res0.4.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(epiapp3.40, features = "MUC5B", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MUC5Bres0.4.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(immapp3res.60, features = "CCL3", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_CCL3res06.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(epiapp3.60, features = "MUC5AC", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MUC5ACres0.6.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(epiapp3.60, features = "MUC5B", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MUC5Bres0.6.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(epiapp3.60, features = "SCGB1A1", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_SCGB1A1res0.6.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")


#=========================================================================================================================================================
# Stromal Cell Category
#=========================================================================================================================================================

obj = SCTransform(stroma, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)

obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.45, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path <- file.path(plot_save_dir, 'UMAPstromapp3_updated')
saveRDS(obj, file = file_path)

#=======
# Re-clustering:
#=======
stromal= readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/UMAPstromapp3_updated.rds')
stromal$cell_types = NULL
stromal_cellTypeAnnotations = c('Adventitial Fibroblasts', 'Smooth Muscle', 'Alveolar Fibroblasts', 'Mesothelium', 'Myofibroblasts', 'Remove', 'Remove', 'Adventitial Fibroblasts')
stromal$cell_types = stromal_cellTypeAnnotations[stromal$seurat_clusters]
strom = subset(stromal, idents = c(1, 2, 3, 4, 5, 8))


obj = SCTransform(strom, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = 'harmony',
                      assay = "SCT", verbose = FALSE)
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.45, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path1 <- file.path(plot_save_dir, 'UMAPStromres0.45_w_removed')
saveRDS(obj, file = file_path1)

#=======
# Generate table and dim plot to view Azimuth fine-level cell type predictions and UMAP layout
#=======
write.table(table(obj$seurat_clusters, obj$predicted.ann_finest_level), col.names = F, file.path(plot_save_dir, 'AzimuthPredictions_Strom_wo_col.csv'))

t2 = DimPlot(obj, group.by='predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(obj, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "Stromal_AzimuthAnnotationsDimplot.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")



#=========================================================================================================================================================
# Endothelial Cell Category
#=========================================================================================================================================================

obj = SCTransform(endo, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)

obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.25, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path <- file.path(plot_save_dir, 'UMAPendoapp3_updated')
saveRDS(obj, file = file_path)

#=======
# Re-clustering:
#=======
endothelial= readRDS('/data/Choi_lung/Elelta/FinalAnnotations/Resolution/FinalResolutions/Approach2Final/UMAPendoapp3_updated.rds')
endothelial$cell_types = NULL
endothelial_cellTypeAnnotations = c('Lymphatic EC', 'EC Arterial', 'EC Venous Pulmonary', 'EC Venous Systemic', 'EC Arterial', 'EC Aerocyte Capillary', 'Lymphatic EC', 'Remove', 'Remove')
endothelial$cell_types = endothelial_cellTypeAnnotations[endothelial$seurat_clusters]
endo = subset(endothelial, idents = c(1, 2, 3, 4, 5, 6, 7))

recluster = function(endo, r, endoremovedclusters)
  obj = SCTransform(endo, vst.flavor = "v2")

obj = SCTransform(endo, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = 'harmony',
                      assay = "SCT", verbose = FALSE)
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.25, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path1 <- file.path(plot_save_dir, 'UMAPendo_w_removed')
saveRDS(obj, file = file_path1)

#=======
# Generate table and dim plot to view Azimuth fine-level cell type predictions and UMAP layout
#=======
write.table(table(obj$seurat_clusters, obj$predicted.ann_finest_level), col.names = F, file.path(plot_save_dir, 'AzimuthPredictions_Endo_wo_col.csv'))

t2 = DimPlot(obj, group.by='predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(obj, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "Endothelial_AzimuthAnnotationsDimplot.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")


#=========================================================================================================================================================
# Immune Cell Category
#=========================================================================================================================================================

obj = SCTransform(imm, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = 'harmony',
                       assay = "SCT", verbose = FALSE)

obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 1.7, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path <- file.path(plot_save_dir, 'UMAPimmune_reclustered1_7')
saveRDS(obj, file = file_path)

#=======
# Generate table and dim plot to view UMAP layout and Azimuth fine-level annotations
#=======
write.table(table(obj$seurat_clusters, obj$predicted.ann_finest_level), col.names = F, file.path(plot_save_dir, 'AzimuthPredictions_Immune_wo_col.csv'))

t2 = DimPlot(obj, group.by='predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(obj, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "Immune_AzimuthAnnotations.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")


#=======
# Visualizing changes to clustering with marker genes
#=======
immune1.7 = readRDS('./UMAPimmune_reclustered1_7')
t2 <- DimPlot(immune1.7, group.by = 'predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(immune1.7, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "immune1_7dimplot.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")

#=======
# Generate feature plots to look at distribution of expression levels of marker genes among cells of cluster 29 to see if multiple cell types naturally group in UMAP space
#=======
imm = subset(immune1.7, idents = c(29))
feature_plot <- FeaturePlot(imm, features = "TYMS", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_TYMS.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(imm, features = "CENPW", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_CENPW.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(imm, features = "MND1", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MND1.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(imm, features = "MKI67", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_MKI67.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

feature_plot <- FeaturePlot(imm, features = "CD3E", cols = c("lightgray", "red"), combine = TRUE)
feature_plot_save_path <- file.path(plot_save_dir, "featureplot_CD3E.png")
ggsave(feature_plot_save_path, feature_plot, width = 6, height = 4, units = "in")

#=======
# Separate out the proliferating NK and T cells based on NK and T- cell marker genes.

# Azimuth automatic annotations had annotated cluster 29 as multiple cell types, which could be due to a small fraction of cells wrongly map to a subgroup due to doublets or mixed subpopulations. Visual inspection of UMAP and surface-marker expression can reveal such cells. It was hard to separate out proliferating NK and T from the proliferating Mph even the cluster was clearly split on the UMAP layout and in marker gene expression. 
#=======
immune <- readRDS('/data/Choi_lung/Elelta/FinalizedCellAnnotation0412/UMAPimmune_reclustered1_7')

#===
# using marker gene expression to separate out proliferating CD3+ subtypes labelled by Azimuth Automatic Prediction and spotted to be grouped together on UMAP layout 
#===
# define the cluster and genes to filter
Immune_clusters_to_filter <- c(29)
Immune_genes_to_filter <- c('CD3E', 'CD3G', 'CD3D', 'NKG7', 'GNLY', 'SPON2')
                            
# selecting out cells of interest
immune_tofilter <- subset(immune, idents = Immune_clusters_to_filter)

# further subset to select cells based on gene expression thresholds
removeimmune <- subset(immune_tofilter, subset = CD3E > 0 | CD3G > 0 | CD3D > 0| NKG7 > 0 | GNLY > 0 | SPON2 > 0)


# saving the removed and remaining cells to separate RDS Files:
cells_to_remove <- Cells(removeimmune)
remaining_cells <- setdiff(Cells(immune), cells_to_remove)

immune_filtered <- subset(immune, cells = remaining_cells)
saveRDS(immune_filtered, file = file.path(plot_save_dir, "filteredimmune1.7.rds"))

immune_removed <- subset(immune, cells = cells_to_remove)
saveRDS(immune_removed, file = file.path(plot_save_dir, "removedimmune1.7.rds"))


#=======
# Re-clustering subset of proliferating T and NK cell types
#=======
immune_removed <- readRDS(file.path(plot_save_dir, "removedimmune1.7.rds"))
recluster = function(immune_removed, r)
obj = SCTransform(immune_removed, vst.flavor = "v2")
obj = RunPCA(obj, npcs = 50, verbose = F)
obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                      orig.reduction = "pca", new.reduction = 'harmony',
                      assay = "SCT", verbose = FALSE)
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, method = 'igraph', resolution = 0.1, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
file_path1 <- file.path(plot_save_dir, 'UMAPCD3') #This was then annotated for proliferating T and NK cell types
saveRDS(obj, file = file_path1)

#=======
# Generate table and dim plot to view UMAP layout and Azimuth fine-level annotations
#=======
write.table(table(obj$seurat_clusters, obj$predicted.ann_finest_level), col.names = F, file.path(plot_save_dir, 'AzimuthPredictions_ImmuneCD3_wo_col.csv'))

t2 = DimPlot(obj, group.by='predicted.ann_finest_level', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
t1 <- DimPlot(obj, group.by = 'seurat_clusters', raster=FALSE, raster.dpi=c(1024,1024), label = TRUE, repel = TRUE)
ggsave(filename = file.path(plot_save_dir, "ImmuneCD3_AzimuthAnnotations.tiff"), plot = t1 +t2, width = 20, height = 20, device = "tiff")



                    


