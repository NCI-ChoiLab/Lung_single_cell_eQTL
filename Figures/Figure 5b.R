#Author: Thong Luong
#Date: December 23rd 2025

library(ComplexHeatmap)

int_effect = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/int_eqtl_effect.rds')
int_effect = t(scale(t(int_effect)))

get_sig_time = function(d){
  tested_df = d
  tested_df = as.data.frame(t(tested_df))
  tested_df$Pseudotime = c(1,2,3,4,5,6)
  df = as.data.frame(matrix(nrow = ncol(tested_df)-1, ncol = 2))
  rownames(df) = colnames(tested_df)[1:(ncol(tested_df)-1)]
  colnames(df) = c('r.squared','pval')
  df_quad = df
  
  for (i in 1:(ncol(tested_df)-1)){
    linear = lm(formula = tested_df[,i] ~ Pseudotime, data = tested_df)
    df[colnames(tested_df)[i],'r.squared'] = summary(linear)$r.squared
    df[colnames(tested_df)[i],'pval'] = summary(linear)$coefficients[,4][['Pseudotime']]
  }
  
  sig_ln_df = subset(df, pval < 0.05)
  print(nrow(sig_ln_df))
  
  for (i in 1:(ncol(tested_df)-1)){
    quad = lm(formula = tested_df[,i] ~ Pseudotime + I(Pseudotime^2), data = tested_df)
    df_quad[colnames(tested_df)[i],'r.squared'] = summary(quad)$r.squared
    df_quad[colnames(tested_df)[i],'pval'] = summary(quad)$coefficients[,4][['Pseudotime']]
  }
  
  sig_quad_df = subset(df_quad, pval < 0.05)
  print(nrow(sig_quad_df))
  
  both = intersect(rownames(sig_quad_df),rownames(sig_ln_df))
  linear_name = setdiff(rownames(sig_ln_df),both)
  quad_name = setdiff(rownames(sig_quad_df),both)
  signif_name = unique(c(rownames(sig_quad_df),rownames(sig_ln_df)))
  
  sig_ln = as.data.frame(subset(d, rownames(d) %in% linear_name))
  sig_quad = as.data.frame(subset(d, rownames(d) %in% quad_name))
  both_df = as.data.frame(subset(d, rownames(d) %in% both))
  non_sig = as.data.frame(subset(d, !(rownames(d) %in% signif_name)))
  
  sig_ln$Model = 'Linear'
  sig_quad$Model = 'Quadratic'
  both_df$Model = 'Both'
  non_sig$Model = 'Neither'
  print(nrow(both_df))
  print(nrow(non_sig))
  
  tested_df = rbind(sig_ln,both_df)
  tested_df = rbind(tested_df,sig_quad)
  tested_df = rbind(tested_df,non_sig)
  return(tested_df)
}


int_effect = get_sig_time(int_effect)

sig_effect = subset(int_effect, Model != 'Neither')


sig_effect$Model = factor(sig_effect$Model, levels = c('Linear','Both','Quadratic'))

ha2 = rowAnnotation(bar = sig_effect$Model, 
                    col = list(bar = c("Linear" = "blue", "Both" = "green", "Quadratic" = "yellow")),
                    annotation_legend_param = list(title = "Model"),
                    show_annotation_name= F)

final = Heatmap(as.matrix(sig_effect[,c(1:6)]), name = 'Relative Effect', cluster_columns = F, show_row_dend = F,
                width = unit(60,'mm'), height = unit(80,'mm'),
                show_column_dend = F, show_row_names = F,cluster_row_slices = F, row_split = sig_effect$Model, left_annotation = ha2, row_title = c('','',''),
                column_names_gp = gpar(fontsize = 5))

c = draw(final)
w = ComplexHeatmap:::width(c)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(c)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)
pdf('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_5/heatmap.pdf', w = w, h = h)
draw(final)
dev.off()

