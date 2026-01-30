#Author: FQ

SCdata <- readRDS("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/Seurat_file/epithelial_cells.rds")
age_gpc_cov <- as.data.frame(t(read.table("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/Covariates/age_gpc_cov.tsv")))
ePC <- data.table::fread("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/Covariates/ePCs.tsv", sep=",")
ePC_trans <- t(ePC[, -1])
colnames(ePC_trans) <- ePC$V1

Countdata <- SCdata@assays$RNA@layers$counts
library(Matrix)
Expre_matrix <- as(Countdata, "dgCMatrix")
Meta <- SCdata@meta.data
gene_info <- as.data.frame(SCdata@assays$RNA@features@.Data)
gene_names <- row.names(as.data.frame(SCdata@assays$RNA@features@.Data))
gene_annotation <- data.frame(id=gene_names, gene_short_name=gene_names)
rownames(Expre_matrix) <- gene_names
colnames(Expre_matrix) <- rownames(Meta)
rownames(gene_annotation) <- gene_annotation$id

library(monocle3)
cds <- new_cell_data_set(Expre_matrix,
                         cell_metadata = Meta,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
#save(cds, file="/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/cds.RData")
#cds <- get(load("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/cds.RData"))

library(ggplot2)
pdf("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/Trajectory Lung epithelial.pdf")
plot_cells(cds,
           color_cells_by = "cell_types",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           trajectory_graph_color = "grey68",
           label_branch_points=F,
           group_label_size=4)
dev.off()

## We can sepecify the start point for pseudotime
cds <- order_cells(cds, root_pr_nodes="Y_1063")
save(cds, file="/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/epithelial_cells_pseudotime_root_Y_1063.RData")
#cds <- get(load("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/epithelial_cells_pseudotime_root_Y_1063.RData"))
pdf("/gpfs/gsfs12/users/qinf2/Jiyeon/Trajectory/Trajectory Lung epithelial pseudotime.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           trajectory_graph_color = "grey68",
           label_branch_points=F,
           group_label_size=4)
dev.off()
