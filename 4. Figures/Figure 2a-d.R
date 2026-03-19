#Author: Elelta Sisay
#Modified: Thong Luong

#===========================
# R script used to generate dim plots for Main Fig. 2a and Supplementary Fig. 2
# Dim plots were generated for the four cell categories and the final annotated seurat object
#===========================

# set output directory
plot_save_dir <- '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/'

library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(reshape2)
library(extrafont)

#===================
#Whole Group
#===================


#Figure 2a
# generate dim plot for Seurat object with initial clustering 
finalcombined = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/total_with_preinfo.rds')


# generate dim plot for final object with cell type annotations
 
# cell types listed just for reference when selecting colors and ordering 

listofcelltypes <- c('Alveolar Transitional Cells','AT1', 'AT2', 'Club', 'Goblet', 'Basal', 'Multiciliated', 'Neuroendocrine', 'Secretory Transitional Cells',
                     'Alveolar Mph', 'Proliferating Mph', 'NK Cells', 'Proliferating NK', 'CD8+ T Cells', 'CD4+ T Cells', 'Proliferating T', 'Interstitial Mph Perivascular', 'Monocyte-derived Mph', 'Classical Monocytes', 'Non-classical Monocytes', 'DC1', 'DC2', 'Migratory DCs', 'Mast Cells', 'B Cells', 'Plasma Cells', 'Plasmacytoid DCs', 
                    'Lymphatic EC Mature', 'Lymphatic EC Proliferating', 'EC Arterial', 'EC Venous Pulmonary', 'EC Venous Systemic', 'EC Aerocyte Capillary', 'EC General Capillary', 
                     'Adventitial Fibroblasts', 'Smooth Muscle', 'Subpleural Fibroblasts', 'Peribronchial Fibroblasts', 'Alveolar Fibroblasts', 'Mesothelium', 'Myofibroblasts')

cols <- c('darkseagreen1' ,'lawngreen', 'green2','seagreen4', 'olivedrab', 'green1', 'palegreen2', 'limegreen', 'darkgreen', 
'lightblue', 'royalblue1', 'slategray3', 'cyan1','steelblue', 'aquamarine4', 'cyan3', 'lightblue', 'darkturquoise', 'cornflowerblue', 'dodgerblue4', 'aquamarine', 'cadetblue', 'slateblue1','cornflowerblue','cadetblue3','lightsteelblue3','cyan2',
'purple1','orchid1','mediumpurple3','mediumorchid2','darkviolet','hotpink','maroon',
'orange', 'goldenrod2', 'gold', 'tan', 'darkorange3', 'sandybrown', 'lightsalmon') 

colors = data.frame(cbind(listofcelltypes,cols))
colnames(colors)[1] = 'Cell_Type'

DimPlot(object = finalcombined, cols = cols , raster = FALSE, label = TRUE, repel = FALSE) + NoLegend()
  ggsave(filename = file.path(plot_save_dir, "Dimplot_Wholegroup_labelled.tiff"), width = 20, height = 20, device = "tiff") 

Idents(finalcombined) = finalcombined$cell_types
Idents(finalcombined) = factor(Idents(finalcombined), levels = listofcelltypes)
  
t = DimPlot(object = finalcombined, cols = cols, raster = T, label = F) + NoLegend() + NoAxes() + theme(panel.background = element_rect(fill = "transparent",colour = NA_character_))
ggsave(filename = file.path(plot_save_dir, "Figure2a_raster_t.png"), plot = t, units = 'mm', width = 100, height = 100, dpi = 300, device = "png", bg = 'transparent') 

#Figure 2b

group_data = data.frame(table(finalcombined$cell_types))
colnames(group_data) = c('Cell_Type','cells')
rownames(group_data) = group_data$Cell_Type
group_data = group_data[listofcelltypes,]
group_data$Subgroup = c(rep('Epithelial',9),rep('Immune',18),rep('Endothelial',7), rep('Stromal',7))
epi = subset(group_data,Subgroup == 'Epithelial')
epi = epi[order(epi$cells),]
imm = subset(group_data,Subgroup == 'Immune')
imm = imm[order(imm$cells),]
endo = subset(group_data,Subgroup == 'Endothelial')
endo = endo[order(endo$cells),]
strom = subset(group_data,Subgroup == 'Stromal')
strom = strom[order(strom$cells),]
group_data = rbind(strom,endo)
group_data = rbind(group_data,imm)
group_data = rbind(group_data,epi)
order = as.character(group_data[,'Cell_Type'])

group_data$Cell_Type = factor(group_data$Cell_Type, levels = order)

ymax = 35000

t = ggplot(group_data, aes(x = cells, y = Cell_Type,fill = Cell_Type)) + 
  geom_col(position = 'dodge') +
  geom_text(aes(label = cells), vjust = .5, hjust = -.05, size = 5, size.unit = 'pt') + 
  coord_cartesian(xlim = c(0,43000)) + 
  xlab('Cell Type') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x = element_text(color = "black", size = 5, angle = 0, hjust = .5, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = .5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        panel.grid.major = element_blank(), # Removes major grid lines
        panel.grid.minor = element_blank(), # Removes minor grid lines
        legend.position = 'none') +
  scale_x_continuous(expand = c(0,0), labels = scales::label_comma()) +
  scale_fill_manual(values = colors$cols)  


ggsave(filename = file.path(plot_save_dir, 'Figure_2b.tiff'), plot = t, width = 75, height = 90, units = "mm", dpi = 300, device = 'tiff')


#Figure 2c

data = data.frame(cell_types = c('Multiciliated','AT2','AT1','Club','Alv_trans','Sec_trans','Goblet','Basal',
                                 'Alv_mph','NK','Cla_mono','Noncla_mono','CD4','Mono_mph','CD8','DC2','Int_mph_peri','DC1','Mig_DC','Bcells','Mast',
                                 'Lym_EC_mat','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_pro','EC_gen_cap','EC_aero_cap',
                                 'SM','Adv_fib','Alv_fib','Peri_fib','Sub_fib'),
                  Subgroup = c(rep("Epithelial", 8), rep("Immune", 13), rep("Endothelial", 7), rep("Stromal", 5))
)

#generate dataframe from tensor results
data$tensorqtl = NA
for (c in data$cell_types){
  t = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  data[which(data$cell_types == c),'tensorqtl'] = length(unique(t$phenotype_id))
  rm(t)
}

cols <- c('darkseagreen1' , 'green4', 'green1', 'seagreen4','olivedrab', 'palegreen2', 'limegreen', 'darkgreen',
          'royalblue1', 'dodgerblue4', 'steelblue4', 'blue', 'aquamarine4', 'lightblue', 'darkturquoise', 
          'cornflowerblue', 'dodgerblue', 'skyblue', 'cadetblue', 'midnightblue', 'steelblue',
          'purple1',
          'orchid1','mediumpurple3','mediumorchid2','darkviolet','magenta','orchid4', 
          'tan', 'lightsalmon', 'sienna', 'gold', 'goldenrod3')

melt_and_graph = function(d,l){
  # Melt the data frame and include 'Subgroup' in id.vars
  melt_d = melt(d, id.vars = c('cell_types', 'Subgroup'), variable.name = 'tensorqtl', value.name = 'eGenes')
  
  # Reorder cell types in a specific order
  melt_d$cell_types = factor(melt_d$cell_types, 
                             levels = c('Multiciliated','AT2','AT1','Club','Alv_trans','Sec_trans','Goblet','Basal',
                                        'Alv_mph','NK','Cla_mono','Noncla_mono','CD4','Mono_mph','CD8','DC2','Int_mph_peri','DC1','Mig_DC','Bcells','Mast',
                                        'Lym_EC_mat','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_pro','EC_gen_cap','EC_aero_cap',
                                        'SM','Adv_fib','Alv_fib','Peri_fib','Sub_fib'))
  
  # Define color assignment for each cell type based on its Subgroup
  subgroup_colors = list(
    Epithelial = cols[1:8],    # Colors for Epithelial subgroup
    Immune = cols[9:21],       # Colors for Immune subgroup
    Endothelial = cols[22:28], # Colors for Endothelial subgroup
    Stromal = cols[29:33]      # Colors for Stromal subgroup
  )
  
  # Create a named color vector directly
  color_map <- c(
    setNames(subgroup_colors$Epithelial, d$cell_types[d$Subgroup == "Epithelial"]),
    setNames(subgroup_colors$Immune, d$cell_types[d$Subgroup == "Immune"]),
    setNames(subgroup_colors$Endothelial, d$cell_types[d$Subgroup == "Endothelial"]),
    setNames(subgroup_colors$Stromal, d$cell_types[d$Subgroup == "Stromal"])
  )
  
  # Create the plot
  my_plot <- melt_d %>%
    ggplot(aes(x = cell_types, y = eGenes, fill = cell_types)) + 
    geom_col(position = 'dodge') +
    geom_text(aes(label = eGenes), vjust = -0.25, size = 5, size.unit = 'pt') + 
    coord_cartesian(ylim = c(0, l * 1.1)) + 
    ylab('# Significant eGenes (q < 0.05)') + 
    theme(text = element_text(family = 'Arial'),
          axis.text.x = element_text(color = "black", size = 5, angle = 45, hjust = 1, vjust = 1, face = "plain"),
          axis.text.y = element_text(color = "black", size = 5, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          panel.grid.major = element_blank(), # Removes major grid lines
          panel.grid.minor = element_blank(), # Removes minor grid lines
          legend.position = 'none') +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = color_map)   # Map the correct fill color
  return(my_plot) 
}

# Call the function with your data
t = melt_and_graph(data, 800)
ggsave(filename = file.path(plot_save_dir, 'Figure_2d.png'), plot = t, width = 110, height = 50, units = "mm", dpi = 300, device = 'png')




