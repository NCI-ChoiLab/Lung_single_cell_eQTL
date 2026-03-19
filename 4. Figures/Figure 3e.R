#Author: Thong Luong
#Date: May 29th 2025

library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Overlap_eQTL_cCRE/Shared_Unique/Mashr/')

unique_top = readRDS('unique.rds')
unique_top$chr = str_split_fixed(unique_top$variant_id, '_',4)[,1]
unique_top$pos = as.numeric(str_split_fixed(unique_top$variant_id, '_',4)[,2])
unique_top = unique_top[!(duplicated(unique_top$qtl)),]


multi_top = readRDS('multi4.rds')
multi_top$chr = str_split_fixed(multi_top$variant_id, '_',4)[,1]
multi_top$pos = as.numeric(str_split_fixed(multi_top$variant_id, '_',4)[,2])
multi_top = multi_top[!(duplicated(multi_top$qtl)),]

excl_df = data.frame(matrix(vector(), nrow(unique_top),2, 
                            dimnames = list(c(1:nrow(unique_top)), c('Exclusive', 'Distance'))),
                     stringsAsFactors = F)

excl_df$Exclusive = 'Unique'
excl_df$Distance = unique_top$start_distance

non_df = data.frame(matrix(vector(), nrow(multi_top),2, 
                           dimnames = list(c(1:nrow(multi_top)), c('Exclusive', 'Distance'))),
                    stringsAsFactors = F)

non_df$Exclusive = 'Multi-categories'
non_df$Distance = multi_top$start_distance

dist_df = rbind(excl_df,non_df)
dist_df$Exclusive = factor(dist_df$Exclusive, levels = c('Unique','Multi-categories'))

t = ggplot(data=dist_df, aes(x=Distance, group=Exclusive, fill=Exclusive)) +
  geom_density(adjust=.5, alpha=.2, linewidth = .1) + ylab('Density') +
  theme(text = element_text(family = 'Arial'),
        axis.text.x =element_text(size=5),
        axis.text.y =element_text(size=5),
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(2, "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.background = element_rect(fill = 'white', colour = 'black', linetype = 1),
        panel.grid.major = element_line(color = "grey", size = .1, linetype = 1),
        legend.key = element_rect(color = NA, fill = NA))

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/'

ggsave(filename = file.path(plot_save_dir, 'Figure_3e.png'), plot = t, width = 65, height = 45, units = "mm", dpi = 450, device = 'png')

print_dist = function(d){
  dist = subset(d, Distance < 100 & Distance > -100)
  print(nrow(dist)/nrow(d))
  dist = subset(d, Distance < 500 & Distance > -500)
  print(nrow(dist)/nrow(d))
  dist = subset(d, Distance < 1000 & Distance > -1000)
  print(nrow(dist)/nrow(d))
  dist = subset(d, Distance < 5000 & Distance > -5000)
  print(nrow(dist)/nrow(d))
  dist = subset(d, Distance < 10000 & Distance > -10000)
  print(nrow(dist)/nrow(d))
  dist = subset(d, Distance < 50000 & Distance > -50000)
  print(nrow(dist)/nrow(d))
}

print_dist(excl_df)
print_dist(non_df)
