
library(e1071)
library(cluster)
library(factoextra)

regulon = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/regulons_qt.rds')
regulon = t(scale(t(regulon)))
int_effect = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/NBME_trajectory/Q_pseudotime_output/int_eqtl_effect.rds')
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

regulon = get_sig_time(regulon)
regulon = subset(regulon, Model == 'Linear')

int_effect = get_sig_time(int_effect)
int_effect = subset(int_effect, Model == 'Linear')
int_effect$Gene = str_split_fixed(rownames(int_effect),'\\|',2)[,1]

incidmat = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/Trajectory/incidMat.rds')
incidmat = incidmat[rownames(regulon),intersect(int_effect$Gene,colnames(incidmat))]
incidmat = incidmat[rowSums(incidmat) > 0, colSums(incidmat) > 0] 

regulon_sub = regulon[rownames(incidmat),]
regulon_sub$Model = NULL
fuzzy_cluster = function(d){
  data_for_clustering <- d
  
  set.seed(123)
  n_cluster <- 5
  m <- 2
  result <- cmeans(data_for_clustering, centers = n_cluster, m = m)
  
  fuzzy_membership_matrix <- result$membership
  
  initial_centers <- result$centers
  final_centers <- t(result$centers)
  
  cluster_membership <- as.data.frame(result$membership)
  max_col_names <- as.data.frame(apply(cluster_membership, 1, function(x) names(x)[which.max(x)]))
  colnames(max_col_names) = 'fuzzy'
  
  data = cbind(d,max_col_names)
  data = data[order(data$fuzzy),]
  data$fuzzy = as.numeric(data$fuzzy)
  
  mat = as.matrix(data)
  return(mat)
}

regulon_sub = fuzzy_cluster(regulon_sub)
mat_qt = as.matrix(regulon_sub)

Activity = colorRamp2(c(-2,0,2), c('darkblue','white','red'))

h1 = Heatmap(mat_qt[,c(1:6)], name = 'Relative Activity', col = Activity, show_row_dend = F, cluster_rows = F, cluster_columns = F, column_title = NULL, row_title = NULL,
             show_column_names = T, row_split = mat_qt[,7], width = unit(60,'mm'), height = unit(100, 'mm'), show_row_names = T, 
             row_names_side = 'left',
             column_names_gp = gpar(fontsize = 5),
             row_names_gp = gpar(fontsize=5),
             show_heatmap_legend = F)

library(circlize)


lgd = Legend(
  col_fun = Activity, 
  title = "Relative Activity", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)


b = draw(h1)
w = ComplexHeatmap:::width(b)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(b)
h = convertY(h, "inch", valueOnly = TRUE)
h2 = ComplexHeatmap:::height(lgd)
h2 = convertY(h2, "inch", valueOnly = TRUE)
h = h+h2
c(w, h)
png('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/Fig_9.png', w = w, h = (h+.2), units = 'in', res = 450)
draw(h1)
draw(
  lgd, 
  x = unit(0.2, "npc"),  # Adjust x coordinate (e.g., 5% from left)
  y = unit(1, "npc"),  # Adjust y coordinate (e.g., 95% from bottom/5% from top)
  just = c("left", "top") # Align the left-top corner of the legend with the x,y coords
)
dev.off()







