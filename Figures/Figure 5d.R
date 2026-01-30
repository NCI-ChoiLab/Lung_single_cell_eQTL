#Author: Thong Luong
#Date: August 19th 2025
#Modified: December 19th 2025

library(ggplot2)


pseudotime = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/pseudotime_finite.rds')
expmat = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/exprMat.rds')
expmat = expmat[,rownames(pseudotime)]
expmat_norm = NormalizeData(expmat)
#expmat_norm2 =  apply(expmat_norm, 1, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
expmat_norm2 = as.data.frame(t(expmat_norm))
expmat_norm2 = expmat_norm2[rownames(pseudotime),]
expmat_norm2$Quantile = pseudotime$q

regulon = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/regulonAUC.rds')
regulon = as.data.frame(regulon@assays@data@listData$AUC)
regulon_scaled = t(scale(t(regulon), center = T, scale=T))
regulon_scaled = as.data.frame(t(regulon_scaled))
regulon_scaled = regulon_scaled[rownames(pseudotime),]
regulon_scaled$Quantile = pseudotime$q


library(wesanderson)
library(magrittr)

cebpd = ggplot(regulon_scaled, aes(x=Quantile, y=regulon_scaled[,'CEBPD(+)'], color=Quantile)) +
  geom_violin(trim=FALSE,linewidth = 0.5)+
  #geom_point(position=position_jitterdodge(),aes(group=snp_ra)) +
  geom_boxplot(aes(middle=median(regulon_scaled[,'CEBPD(+)'])),
               width=0.1, linewidth = 0.5, outlier.size = .1)+
  #geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
  #method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
  labs(title="",x='', y = "Scaled Regulon Activity")+
  scale_color_manual(values=wes_palette(n=6, name="AsteroidCity2")) +
  ggtitle('CEBPD(+) Regulon Activity')+
  theme(text = element_text(family = 'Arial'),
        plot.title = element_text(size = 7, color = 'black', face = 'plain'),
        axis.text.x = element_text(color = "black", size = 6,face = "plain"),
        axis.text.y = element_text(color = "black", size = 6,face = "plain"),  
        axis.title.x = element_text(color = "black", size = 6, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, face = "plain"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent'),
        axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
        legend.position="none")

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_5/'
ggsave(filename = file.path(plot_save_dir, 'cebpd.png'), plot = cebpd, width = 70, height = 45, units = "mm", dpi = 450, device = 'png')

graph_expression = function(goi){
  df = expmat_norm2[,c('Quantile',goi)]
  df[,goi] = scale(df[,goi], center = T, scale = T)
  ggplot(df, aes(x=Quantile, y=df[,goi], color=Quantile)) +
    geom_violin(trim=FALSE,linewidth = 0.25)+
    #geom_point(position=position_jitterdodge(),aes(group=snp_ra)) +
    geom_boxplot(aes(middle=median(df[,goi])),
                 width=0.1, linewidth = 0.25, outlier.size = .1)+
    #geom_smooth(mapping = aes(x = snp_ra, y = isoform, group = 1),formula = y~x, color = "gray",
    #method='lm', size = 1, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="FOXP2",x='', y = paste0("Scaled Expression of ",goi))+
    scale_color_manual(values=wes_palette(n=6, name="IsleofDogs1")) +
    ggtitle(goi)+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 7, color = 'black', face = 'plain'),
          axis.text.x = element_text(color = "black", size = 6,face = "plain"),
          axis.text.y = element_text(color = "black", size = 6,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 6, face = "plain"),
          axis.title.y = element_text(color = "black", size = 6, face = "plain"),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent'),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")
  
}

foxp2 = graph_expression('FOXP2')
ggsave(filename = file.path(plot_save_dir, 'foxp2.png'), plot = foxp2, width = 70, height = 45, units = "mm", dpi = 450, device = 'png')
