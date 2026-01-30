#Author: Thong Luong
#Date: June 4th 2025

library(stringr)
library(QTLExperiment)
library(Seurat)
library(patchwork)
library(Signac)
library(ggplot2)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Comparison/')

multiome = readRDS('/data/Choi_lung/ChiaHan/CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds')
Idents(multiome) = 'CellType'

#rs1176359
CoveragePlot(multiome, region = 'chr11-33934432-33934873', extend.upstream = 50000, extend.downstream = 50,
             region.highlight = StringToGRanges(c("chr11-33934485-33934685")),
             features = 'LMO2', assay = 'peaks', expression.assay = 'SCT', peaks = F, links = F, annotation = F)


lmo2 = CoveragePlot(multiome, region = "chr11-33934585-33934585", extend.upstream = 5000, extend.downstream = 5000,
                    region.highlight = StringToGRanges(c("chr11-33934535-33934635")),
                    features = 'LMO2', assay = 'peaks', expression.assay = 'SCT', peaks = F, links = F, annotation = F)

expr_plot <- ExpressionPlot(
  object = multiome,
  features = "LMO2",
  assay = "SCT"
) 


gene_plot <- AnnotationPlot(
  object = multiome,
  region = 'chr11-33934432-33934873', extend.upstream = 50000, extend.downstream = 50
)


lmo2_t = lmo2 + theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank())

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/'


ggsave(filename = file.path(plot_save_dir, 'Figure_3d1.1.pdf'), plot = lmo2_t, width = 90, height = 150, units = "mm", dpi = 450, device = 'pdf')
ggsave(filename = file.path(plot_save_dir, 'Figure_3d2.1.pdf'), plot = expr_plot, width = 60, height = 150, units = "mm", dpi = 450, device = 'pdf')
ggsave(filename = file.path(plot_save_dir, 'Figure_3d3.1.pdf'), plot = gene_plot, width = 60, height = 50, units = "mm", dpi = 450, device = 'pdf')

#rs3138486
CoveragePlot(multiome, region = 'chr9-89603007-89607153', extend.upstream = 1000, extend.downstream = 1000,
             region.highlight = StringToGRanges(c('chr9-89603409-89603429')),
             features = 'GADD45G',assay = 'peaks', expression.assay = 'SCT', peaks = F, links = TRUE)



gadd45g = CoveragePlot(multiome, region = "chr9-89603419-89603419", extend.upstream = 4500, extend.downstream = 5500,
                    region.highlight = StringToGRanges(c("chr9-89603369-89603469")),
                    features = 'GADD45G', assay = 'peaks', expression.assay = 'SCT', peaks = F, links = F, annotation = F)

expr_plot <- ExpressionPlot(
  object = multiome,
  features = "GADD45G",
  assay = "SCT"
) 

expr_plot + coord_cartesian(xlim = c(0,2))
gene_plot <- AnnotationPlot(
  object = multiome,
  region = 'chr9-89603419-89603419', extend.upstream = 200, extend.downstream = 3000,
)


gadd45g_t = gadd45g + theme(axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),axis.line.y = element_blank())

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/'


ggsave(filename = file.path(plot_save_dir, 'Figure_3d1.2.pdf'), plot = gadd45g_t, width = 90, height = 150, units = "mm", dpi = 450, device = 'pdf')
ggsave(filename = file.path(plot_save_dir, 'Figure_3d2.2.pdf'), plot = expr_plot, width = 60, height = 150, units = "mm", dpi = 450, device = 'pdf')
ggsave(filename = file.path(plot_save_dir, 'Figure_3d3.2.pdf'), plot = gene_plot, width = 60, height = 50, units = "mm", dpi = 450, device = 'pdf')


#rs6489721
CoveragePlot(multiome, region = 'chr12-6533220-6535753', extend.upstream = 1000, extend.downstream = 1000,
             region.highlight = StringToGRanges(c('chr12-6534140-6534160')),
             features = 'GAPDH',assay = 'peaks', expression.assay = 'SCT', peaks = TRUE, links = TRUE)

#rs3741918
CoveragePlot(multiome, region = 'chr12-6533220-6535753', extend.upstream = 1000, extend.downstream = 1000,
             region.highlight = StringToGRanges(c('chr12-6535080-6535100')),
             features = 'GAPDH',assay = 'peaks', expression.assay = 'SCT', peaks = TRUE, links = TRUE)

