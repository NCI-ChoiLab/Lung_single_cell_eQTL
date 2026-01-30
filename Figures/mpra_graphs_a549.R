#Author: Thong Luong

library(ggplot2)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/MPRA')

load('/data/Choi_lung/MPRA_Lung/regression_A549/Lung_MPRA_data.RData')


snps_of_interest = c('6:117780158','6:117782634','6:117784658')
MPRA_data = subset(MPRA_data, Coord_hg19 %in% snps_of_interest)
MPRA_data = as.data.frame(MPRA_data)

MPRA_data$Strand = factor(MPRA_data$Strand, levels = c('fwd','rev'))
MPRA_data_cleaned = na.omit(MPRA_data)
MPRA_data_cleaned = subset(MPRA_data_cleaned, group == 'DMSO')
MPRA_data_cleaned$class = factor(MPRA_data_cleaned$class, levels = c('fwd_alt','fwd_ref','rev_alt','rev_ref'))


#6:117780158	rs2180811	2.01E-02
#6:117782634	rs9374663	7.53E-04
#6:117784658	rs5879422	1.49E-03

mpra_graph = function(coord, snp, fdr){
  g1 = ggplot(MPRA_data_cleaned[which(MPRA_data_cleaned$Coord_hg19 == coord),], aes(fill=Type, y = Ratio, x=class)) + 
    geom_violin(position = position_dodge(), trim=FALSE,linewidth = 0.25) + 
    geom_boxplot(position = position_dodge(), aes(middle = mean(Ratio)),width=0.1, linewidth = 0.25) + 
    scale_x_discrete(labels = c('alt','ref','alt','ref')) +
    labs(x = '', y = 'RNA TPM/DNA TPM') +
    ggtitle(snp) +
    geom_text(
      x = Inf, y = Inf, # Position at the top right corner of the plot area
      label = paste0('FDR = ',fdr), size = 6, size.unit = 'pt',
      hjust = 1.05,      # Right-align the text
      vjust = 1.2       # Top-align the text
    ) +
    theme(axis.text.x = element_text(color = "black", size = 5,face = "plain"),
          axis.text.y = element_text(color = "black", size = 5,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 6, face = "plain"),
          axis.title.y = element_text(color = "black", size = 6, face = "plain"),
          plot.subtitle = element_text(size = 5),
          plot.title = element_text(hjust = 0.5, size = 7),
          panel.background = element_rect(fill='transparent'),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
          plot.background = element_rect(fill='transparent'),
          plot.caption = element_text(size = 5),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="None")
  return(g1)
}
g1 = mpra_graph('6:117784658','rs5879422','1.49E-03')
plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'rs5879422_mpra.png'), plot = g1, width = 50, height = 50, units = "mm", dpi = 300, device = 'png')

g2 = mpra_graph('6:117782634','rs9374663','7.53E-04')
plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'rs9374663_mpra.png'), plot = g2, width = 50, height = 50, units = "mm", dpi = 300, device = 'png')

