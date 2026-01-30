#Author: Thong Luong
#Date: December 23rd 2025

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(circlize)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/')

load('epithelial_cells_pseudotime_root_Y_1063.RData')

#extracting pseudotime values from all epithelial cells
time = as.data.frame(cds@principal_graph_aux@listData$UMAP$pseudotime)
#extracing cell type
cell_type = as.data.frame(colData(cds)@listData$cell_types)
colnames(cell_type) = 'ct'
rownames(cell_type) = rownames(time)

#dataframe of cell type and pseudotime
epi_pseudotime = as.data.frame(cbind(time,cell_type))
colnames(epi_pseudotime)[1] = 'pseudotime'

#subsetting to just include the three relevant cell types
time_alv = subset(epi_pseudotime, ct %in% c('Alveolar Transitional Cells', 'AT2', 'AT1'))


#subsetting our initial "object" to only keep cells of interest
cds_sub = cds[,rownames(time_alv_finite)]
cds_sub_r <- reduce_dimension(cds_sub, umap.n_neighbors = 15)
umap_coords <- reducedDims(cds_sub_r)$UMAP
colnames(umap_coords) = c('UMAP1','UMAP2')
umap_coords = umap_coords[rownames(time_alv_finite),]

#removing cells with infinite pseudotime values
time_alv_finite = subset(time_alv, pseudotime != 'Inf')
#pseudotime values
time_alv_finite$pseudotime = as.numeric(time_alv_finite$pseudotime)

time_alv_finite = cbind(time_alv_finite,umap_coords)
time_alv_finite %>%
  mutate(Quantile = ntile(time_alv_finite$pseudotime, 6)) -> time_alv_finite
time_alv_finite$Quantile = paste0('Q',time_alv_finite$Quantile)

time_alv_finite$UMAP1 = -(time_alv_finite$UMAP1)

t = ggplot(time_alv_finite,aes(x=UMAP1,y=UMAP2)) + geom_point(aes(color = Quantile), size = .1) + 
  scale_color_manual(values = c("Q1" = "#008000",
                                "Q2" = "#2EB62C",
                                'Q3' = "#57C84D", 
                                'Q4' =  "#83D475", 
                                'Q5' = "#ABE098",
                                'Q6' = "#C5E8B7")) +
  xlab('UMAP1') +
  ylab('UMAP2') + 
  theme(text = element_text(family = 'Arial'),
        axis.text=element_text(size=5),
        axis.text.x = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.title = element_text(size = 6), # Change the title font size
        legend.text = element_text(size = 5),
        legend.key.size = unit(2, "mm")) + guides(colour = guide_legend(override.aes = list(size=1.5)))
plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_5'
ggsave(filename = file.path(plot_save_dir, 'pseudotime.png'), plot = t, width = 80, height = 55, units = "mm", dpi = 450, device = 'png')

