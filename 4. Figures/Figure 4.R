#Author: Thong Luong
#Date: April 14th 2025

library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(circlize)
library(Seurat)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/TWAS_heatmap/')

our_coloc = read.csv('our_coloc.txt', sep = '\t')
our_coloc = our_coloc[c(1:44),]
our_coloc_genes = unique(our_coloc$Gene)

byun_coloc_genes = c('CYP2A6','HLA-J','ZNRD1','CYP2A6','CYP2A6','SECISBP2L','FUBP1','HLA-DRB5','HLA-DRB6','SECISBP2L','HLA-DQB1','HLA-DQB1-AS1','CCHCR1','SKIV2L','FLOT1','TP63','TP63',
                     'IRF4','MPZL3','HLA-DQA1','HLA-DQB1','HLA-DQB2','HLA-DQB1-AS1','C4A','MPZL3','TNXA','CYP21A1P','SFTA2','DCBLD1','NOTCH4','SMPD2','XBP1','LY6G5B','XXbac-BPG154L12.4',
                     'CHRNA3','CTA-292E10.8','CTC-490E21.11','VARS2','STN1','NRG1','LINC00243')
shi_coloc_genes = c('FADS1','TP63','ACVR1B')
zhu_coloc_genes = c('DCAF16','DCBLD1','CBL','MPZL3','NPHP4','ATR','GYPE','PARD3')
fu_coloc_genes = c('SPRED2','CASP8', 'ALS2CR12','RP11-362K14.6','TP63','RP3-512B11.3', 'DSP','WDR46','OARD1', 'RP11-328M4.2','ROS1', 'DCBLD1','AQP3','PRPF19', 'CD6', 'DAGLA', 
                   'FADS1', 'NXF1','JAML', 'MPZL2', 'SLC37A4','ACVR1B','SECISBP2L','BPTF', 'C17orf58','CSRNP2','RP11-73E17.2')

published_coloc_genes = unique(c(byun_coloc_genes,shi_coloc_genes,zhu_coloc_genes))
published_un_coloc_genes = unique(c(byun_coloc_genes,shi_coloc_genes,zhu_coloc_genes,fu_coloc_genes))

our_twas = read.csv('Our_TWAS.txt', sep = '\t')

shi_twas = c('ELF5','FADS1')
zhu_twas = c('DCAF16','CBL','DCBLD1','MPZL3','ATR','GYPE','PARD3','NPHP4','MYNN','TIGD2','AHRR','CEP72','LPCAT1','KIF24','TLDC2','RBL1','CLPTM1L','BPTF')
bosse1_twas = c('IREB2', 'CHRNA5', 'HYKK', 'PSMA4','SECISBP2L','APOM', 'ZNRD1', 'PRRCA', 'ZFP57', 'FLOT1', 'NOTCH4', 'TUBB', 'HLA-A', 'HLA-F', 'HLA-G', 'HCG8', 'HLA-J', 'TRIM38', 'HCP5', 'SFTA2', 'HLA-L', 'HLA-DQB1', 'ZKSCAN3', 
                'PSORS1C3', 'CCHCR1','FGFR10P','JAML','IREB2', 'CHRNA5', 'HYKK', 'PSMA4', 'SECISBP2L', 'GALK2', 'FAM227B','TP63','JAML','NRG1','AQP3','IREB2','APOM', 'NOTCH4', 'PRRC2A, HCP5', 'HLA-DQB1', 'FLOT1', 'ZNRD1', 'MICB', 
                'TRIM38', 'TUBB', 'HLA-A', 'ZKSCAN3', 'HCG8', 'ZNF192P1', 'ZSCAN26', 'HLA-L', 'ZFP67', 'HLA-F', 'HCG14','BLOC1S2','IREB2','HIST1H2BD','TMA16','IREB2', 'CHRNA5', 'HYKK', 'PSMA4','APOM', 'ZNRD1', 'PRR2A', 
                'ZFP57', 'NOTCH4', 'HLA-G', 'TUBB', 'HLA-F', 'HCP5', 'HLA-J', 'HLA-A', 'FLOT1', 'HLA-DQB1', 'CCHCR1','SECISBP2L','NEXN')
bosse2_twas = c('NPHP4','PIGK','DCST2','NUP35','ORMDL1','SUMF1','TP63','PACRGL','SLC6A3','XRCC4','SLC22A5','SPINK1','TRIM38','APOM','DCBLD1','FGFR1OP','CLU','NRG1','MTAP','PPP1R3C','TMEM180',
                'VTI1A','JAML','RAD52','SLC11A2','GOLGA2P5','ATXN2','MIPEP','N4BP2L2','MIR17HG','SECISBP2L','IREB2','BPTF','APCDD1','GAREM','C19orf54','PIGU','CBLN4','XBP1','MTMR3')
bosse_twas = unique(c(bosse1_twas,bosse2_twas))

fu_twas = c('CHIA','ANGEL2','RP11-488L18.10','AC016747.3','IFT122','A4GNT','FAM13A','RP5-874C20.3','ELF5','ZG16B','FAM101B','VTN','GSDMA','NAPSB','ADA','DDT','GST22','ST7L','CLDN18',
            'CYP39A1','PPP6C','OBFC1','CAT','RBM7','PSMD3','DLX3','AC079325.6','JAKMIP1','HMGN3','FLJ00104')


published_twas_genes = unique(c(shi_twas,zhu_twas,bosse_twas))
published_un_twas_genes = unique(c(shi_twas,zhu_twas,bosse_twas,fu_twas))

published_genes = unique(c(published_coloc_genes,published_twas_genes))

our_twas = read.csv('test_twas.txt', sep = '\t')
our_twas = our_twas[c(1:56),]
#our_twas[53,'Novel'] = 'Yes'

celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')
sc = readRDS('/data/Choi_lung/TTL/Figs_for_manuscript/total_with_preinfo.rds')
numb = as.data.frame(table(sc$cell_types))
rownames(numb) = numb$Var1
numb$Var1 = NULL
rownames(numb) = c('Adv_fib','Alv_fib','Alv_mph','Alv_trans','AT1','AT2','Bcells','Basal',
                   'CD4','CD8','Cla_mono','Club','DC1','DC2','EC_aero_cap','EC_art',
                   'EC_gen_cap','EC_ven_pul','EC_ven_sys','Goblet','Int_mph_peri','Lym_EC_mat','Lym_EC_pro','Mast',
                   "Mesothelium",'Mig_DC','Mono_mph','Multiciliated','Myo_fib','Neuroendocrine','NK','Noncla_mono',
                   'Peri_fib','Plasma','Pla_DC','Pro_mph','Pro_NK','Pro_T','Sec_trans','SM','Sub_fib')


epithelial = c("Alv_trans",'AT1','AT2','Basal','Club','Goblet','Multiciliated','Sec_trans')
immune = c("Alv_mph","Bcells","CD4","CD8","Cla_mono",'DC1',"DC2","Int_mph_peri","Mast","Mig_DC",
           "Mono_mph","NK","Noncla_mono")
endothelial = c("EC_aero_cap","EC_art","EC_gen_cap","EC_ven_pul","EC_ven_sys","Lym_EC_mat","Lym_EC_pro")
stromal = c("Adv_fib","Alv_fib","Peri_fib","SM","Sub_fib")


coloc_df = reshape2::dcast(our_coloc, Cell.Type + Study ~ Gene , value.var = "PP.H4")
coloc_df[is.na(coloc_df)] = 0
for (s in unique(coloc_df$Study)){
  study = subset(coloc_df, Study == s)
  rownames(study) = study$Cell.Type
  study = study[,c(3:14)]
  study = study[,c('ACTR2','HLA-DQA1','ROS1','NRG1','TCF7L2','CAT',
                   'JAML','MPZL3','ACVR1B','FAM227B','SECISBP2L','XBP1')]
  assign(s, as.matrix(study))
}
#Figure 4a
PP = colorRamp2(c(0,.8,1), c('white','deeppink','red'))

Shi_EAS = Shi_EAS[c('Alv_trans','AT1','AT2','Club','Multiciliated','Sec_trans', 'Alv_mph','Mono_mph'),]
coloc1 = Heatmap(Shi_EAS, name = 'PP.H4', col = PP, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(Shi_EAS)*3.3),'mm'), height = unit((nrow(Shi_EAS)*3.3),'mm'),
                 column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                 row_names_gp = gpar(fontsize = 5,col = c(rep('green',6),rep('blue',2))),
                 column_names_gp = gpar(fontsize = 5),
                 show_heatmap_legend = F,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(Shi_EAS[i, j] > 0) {
                     grid::grid.text(sprintf("%.2f", Shi_EAS[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                   }
                 })

Byun_Total = Byun_Total[c('Alv_trans','AT2','Club','Multiciliated','Sec_trans', 'Alv_mph','CD4','CD8','Mono_mph','NK'),]
coloc2 = Heatmap(Byun_Total, name = 'PP.H4', col = PP, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(Byun_Total)*3.3),'mm'), height = unit((nrow(Byun_Total)*3.3),'mm'),
                 column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                 row_names_gp = gpar(fontsize = 5, col = c(rep('green',5),rep('blue',5))),
                 column_names_gp = gpar(fontsize = 5), 
                 show_heatmap_legend = F,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(Byun_Total[i, j] > 0) {
                     grid::grid.text(sprintf("%.2f", Byun_Total[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                   }
                 })

Byun_LUAD = Byun_LUAD[c('Alv_trans','AT1','AT2','Multiciliated','Sec_trans', 'Alv_mph','Mono_mph'),]
coloc3 = Heatmap(Byun_LUAD, name = 'PP.H4', col = PP, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(Byun_LUAD)*3.3),'mm'), height = unit((nrow(Byun_LUAD)*3.3),'mm'),
                 column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                 row_names_gp = gpar(fontsize = 5, col = c(rep('green',5),rep('blue',2))),
                 column_names_gp = gpar(fontsize = 5), 
                 show_heatmap_legend = F,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(Byun_LUAD[i, j] > 0) {
                     grid::grid.text(sprintf("%.2f", Byun_LUAD[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                   }
                 })

coloc4 = Heatmap(Byun_LUSC, name = 'PP.H4', col = PP, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(Byun_LUSC)*3.3),'mm'), height = unit((nrow(Byun_LUSC)*3.3),'mm'),
                 column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                 row_names_gp = gpar(fontsize = 5, col = 'blue'),
                 column_names_gp = gpar(fontsize = 5), 
                 show_heatmap_legend = FALSE,
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if(Byun_LUSC[i, j] > 0) {
                     grid::grid.text(sprintf("%.2f", Byun_LUSC[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                   }
                 })

lgd = Legend(
  col_fun = PP, 
  title = "PP.H4", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)

draw(coloc1 %v% coloc2 %v% coloc3 %v% coloc4, ht_gap = unit(.5, "mm"))
draw(
  lgd, 
  x = unit(0.44, "npc"),  # Adjust x coordinate (e.g., 5% from left)
  y = unit(0.81, "npc"),  # Adjust y coordinate (e.g., 95% from bottom/5% from top)
  just = c("left", "top") # Align the left-top corner of the legend with the x,y coords
)

a = draw(coloc1 %v% coloc2 %v% coloc3 %v% coloc4, ht_gap = unit(.5, "mm"))
w = ComplexHeatmap:::width(a)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(a)
h = convertY(h, "inch", valueOnly = TRUE)
h2 = ComplexHeatmap:::height(lgd)
h2 = convertY(h2, "inch", valueOnly = TRUE)
h = h+h2
c(w, h)
png('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_4/Figure_4a.png', w = w, h = (h+.2), units = 'in', res = 450)
draw(coloc1 %v% coloc2 %v% coloc3 %v% coloc4, ht_gap = unit(.5, "mm"))
draw(
  lgd, 
  x = unit(0.23, "npc"),  # Adjust x coordinate (e.g., 5% from left)
  y = unit(1, "npc"),  # Adjust y coordinate (e.g., 95% from bottom/5% from top)
  just = c("left", "top") # Align the left-top corner of the legend with the x,y coords
)
dev.off()

#Figure 4b
our_twas = our_twas[complete.cases(our_twas),]
our_twas_sorted = our_twas[with(our_twas, order(CHR,P0,ID.NAME)),]
twas_df = reshape2::dcast(our_twas_sorted, ID.NAME + Novel ~ CELLTYPE , value.var = "TWAS.Z")
sorted_genes = unique(our_twas_sorted$ID.NAME)
rownames(twas_df) = twas_df$ID.NAME
twas_df = twas_df[sorted_genes,]

twas_df[is.na(twas_df)] = 0
for (n in unique(twas_df$Novel)){
  study = subset(twas_df, Novel == n)
  study = study[,c(3:23)]
  assign(n, as.matrix(study))
}

No1 = subset(No, rownames(No) %in% published_genes)
No2 = subset(No, !(rownames(No) %in% published_genes))
No1 = t(No1)
No2 = t(No2)
Yes = t(Yes)
No1 = as.data.frame(No1)
No2 = as.data.frame(No2)
Yes = as.data.frame(Yes)
No1$ACVR1B = Yes$ACVR1B
No2$CAT = Yes$CAT
Yes$ACVR1B = NULL
Yes$CAT = NULL
No1 = No1[,c("ZNRD1","HLA-DQB1",'HLA-DQA1',"NRG1","JAML","MPZL3","ACVR1B","SECISBP2L","FAM227B")]
No2 = No2[,c('DTNB','PPP1R18','HLA-C','HLA-B','CLIC1','HLA-DRB1','HLA-DQA2','HLA-DPA1', 'CUTA','ROS1', 'TCF7L2', 'FAM53B','CAT','TMEM258')]
No1 = as.matrix(No1)
No2 = as.matrix(No2)
Yes = as.matrix(Yes)

Z.color = colorRamp2(c(min(twas_df[,c(3:23)]),0,max(twas_df[,c(3:23)])),c('blue','white','red'))

No1 = No1[c('EC_aero_cap','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_mat',
           'Alv_trans','AT1','AT2','Basal','Club','Multiciliated','Sec_trans',
           'Alv_mph','CD8','Cla_mono','DC1','Mono_mph','NK','Noncla_mono',
           'Alv_fib','SM'),]
no_hm1 = Heatmap(No1, name = 'TWAS.Z', col = Z.color, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(No1)*3.3),'mm'), height = unit((nrow(No1)*3.3),'mm'),
        column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
        row_names_gp = gpar(fontsize = 5, col = c(rep('purple',5),rep('green',7),rep('blue',7),rep('orange',2))),
        column_names_gp = gpar(fontsize = 5),
        show_heatmap_legend = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(No1[i, j] != 0) {
            grid::grid.text(sprintf("%.1f", No1[i, j]), x, y, gp = grid::gpar(fontsize = 5))
          }
        })

No2 = No2[c('EC_aero_cap','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_mat',
          'Alv_trans','AT1','AT2','Basal','Club','Multiciliated','Sec_trans',
          'Alv_mph','CD8','Cla_mono','DC1','Mono_mph','NK','Noncla_mono',
          'Alv_fib','SM'),]
no_hm2 = Heatmap(No2, name = 'TWAS.Z', col = Z.color, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(No2)*3.3),'mm'), height = unit((nrow(No2)*3.3),'mm'),
                column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                row_names_gp = gpar(fontsize = 5,col = c(rep('purple',5),rep('green',7),rep('blue',7),rep('orange',2))),
                column_names_gp = gpar(fontsize = 5),
                show_heatmap_legend = F,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(No2[i, j] != 0) {
                    grid::grid.text(sprintf("%.1f", No2[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                  }
                })

Yes = Yes[c('EC_aero_cap','EC_art','EC_ven_pul','EC_ven_sys','Lym_EC_mat',
             'Alv_trans','AT1','AT2','Basal','Club','Multiciliated','Sec_trans',
             'Alv_mph','CD8','Cla_mono','DC1','Mono_mph','NK','Noncla_mono',
             'Alv_fib','SM'),]
yes_hm = Heatmap(Yes, name = 'TWAS.Z', col = Z.color, show_row_dend = F, cluster_rows = F,  cluster_columns = F, width = unit((ncol(Yes)*3.3),'mm'), height = unit((nrow(Yes)*3.3),'mm'),
                column_title = NULL, show_column_names = T, show_row_names = T, border = T, row_names_side = 'left',
                row_names_gp = gpar(fontsize = 5, col = c(rep('purple',5),rep('green',7),rep('blue',7),rep('orange',2))),
                column_names_gp = gpar(fontsize=5),
                show_heatmap_legend = F,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(Yes[i, j] != 0) {
                    grid::grid.text(sprintf("%.1f", Yes[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                  }
                })

lgd2 = Legend(
  col_fun = Z.color, 
  title = "TWAS.Z", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)

b = draw(no_hm1 + no_hm2 + yes_hm, ht_gap = unit(.5, "mm"))
w = ComplexHeatmap:::width(b)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(b)
h = convertY(h, "inch", valueOnly = TRUE)
h2 = ComplexHeatmap:::height(lgd2)
h2 = convertY(h2, "inch", valueOnly = TRUE)
h = h+h2
c(w, h)
png('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_4/Figure_4b.png', w = w, h = (h+.2), units = 'in', res = 450)
draw(no_hm1 + no_hm2 + yes_hm, ht_gap = unit(.5, "mm"))
draw(
  lgd2, 
  x = unit(0.12, "npc"),  # Adjust x coordinate (e.g., 5% from left)
  y = unit(1, "npc"),  # Adjust y coordinate (e.g., 95% from bottom/5% from top)
  just = c("left", "top") # Align the left-top corner of the legend with the x,y coords
)
dev.off()

#Figure 4c
cellnumber = data.frame(matrix(vector(), length(celltypes), 1, 
                         dimnames = list(celltypes, 'Cell_number')),
                  stringsAsFactors = F)

for (c in celltypes){
  cellnumber[c,'Cell_number'] = numb[c,'Freq']
}

epi = as.data.frame(subset(cellnumber, rownames(cellnumber) %in% epithelial))
epi$Lineage = 'Epithelial'
epi = epi[order(epi$Cell_number, decreasing = T),]

imm = as.data.frame(subset(cellnumber, rownames(cellnumber) %in% immune))
imm$Lineage = 'Immune'
imm = imm[order(imm$Cell_number, decreasing = T),]

endo = as.data.frame(subset(cellnumber, rownames(cellnumber) %in% endothelial))
endo$Lineage = 'Endothelial'
endo = endo[order(endo$Cell_number, decreasing = T),]

strom = as.data.frame(subset(cellnumber, rownames(cellnumber) %in% stromal))
strom$Lineage = 'Stromal'
strom = strom[order(strom$Cell_number, decreasing = T),]


cellnumber = rbind(endo,epi,imm,strom)


data = data.frame(matrix(vector(), 7, length(celltypes), 
                         dimnames = list(c('Sig_loci','eGenes', 'Coloc_rep','Coloc_novel','TWAS_rep','TWAS_novel_gene','TWAS_novel_loci'), celltypes)),
                  stringsAsFactors = F)

for (c in celltypes){
  egenes = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  egenes = length(unique(egenes$phenotype_id))
  data['eGenes',c] = egenes
  cell_coloc_genes = subset(our_coloc, Cell.Type == c)
  cell_twas_genes = subset(our_twas, CELLTYPE == c)
  uni_twas_genes = subset(cell_twas_genes, !(ID.NAME %in% cell_coloc_genes$Gene))
  data["Sig_loci",c] = length(unique(cell_coloc_genes$Locus)) + nrow(uni_twas_genes)
  data["Coloc_rep",c] = length(intersect(unique(cell_coloc_genes$Gene),published_coloc_genes))
  data["Coloc_novel",c] = length(setdiff(unique(cell_coloc_genes$Gene), published_coloc_genes))
  data["TWAS_rep",c] = length(intersect(unique(cell_twas_genes$ID.NAME),published_twas_genes))
  data["TWAS_novel_gene",c] = length(setdiff(unique(cell_twas_genes$ID.NAME),published_twas_genes))
  uni_twas_genes = subset(uni_twas_genes, Novel == 'Yes')
  data['TWAS_novel_loci',c] = nrow(uni_twas_genes)
}

data = data[,rownames(cellnumber)]
data = as.matrix(data)


col_sig_loci = colorRamp2(c(0,7,15), c('white','grey80','grey40'))

t = as.matrix(t(data['Sig_loci',]))
data_sig_loci = Heatmap(t, name = '# of Significant Loci (Coloc & TWAS)', col = col_sig_loci, show_row_dend = F, cluster_rows = F, cluster_columns = F, column_title = NULL,
                        show_column_names = T, width = unit(ncol(data)*5, 'mm'), height = unit(5,'mm'), border = T, column_split = cellnumber$Lineage, show_row_names = T, row_labels = 'Signficant Loci',
                        column_names_gp = gpar(fontsize = 6, col = c(rep('purple',7),rep('green',8),rep('blue',13),rep('orange',5))),
                        row_names_gp = gpar(fontsize=6),
                        show_heatmap_legend = F,
                        cell_fun = function(j, i, x, y, w, h, fill) {
                          if(t[i,j] > 0) {
                            grid::grid.text(sprintf("%.0f", t[i,j]), x, y, gp = grid::gpar(fontsize = 5))
                          }
                        })

lgd_sig_loci = Legend(
  col_fun = col_sig_loci, 
  title = "# of Signficant Loci (Coloc & TWAS)", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)


coloc_genes = colorRamp2(c(0,2,4), c('white','aquamarine','darkgreen'))

coloc_matrix = as.matrix(data[c(4,5),])
data_coloc = Heatmap(coloc_matrix, name = 'Colocalized Genes', col = coloc_genes, show_row_dend = F, cluster_rows = F, cluster_columns = F, column_title = NULL,
                     show_column_names = T, column_split = cellnumber$Lineage, width = unit(ncol(coloc_matrix)*5, 'mm'), height = unit(nrow(coloc_matrix)*5, 'mm'), show_row_names = T, border = T, row_labels = c('Novel','Replicated'),
                     bottom_annotation = columnAnnotation(eGenes = anno_barplot(data['eGenes',], axis_param= list(gp = gpar(fontsize = 6))) ,annotation_name_gp = gpar(fontsize = 6)),
                     column_names_gp = gpar(fontsize = 6, col = c(rep('purple',7),rep('green',8),rep('blue',13),rep('orange',5))),
                     row_names_gp = gpar(fontsize=6),
                     show_heatmap_legend = F,
                     cell_fun = function(j, i, x, y, w, h, fill) {
                       if(coloc_matrix[i,j] > 0) {
                         grid::grid.text(sprintf("%.0f", coloc_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                       }
                     })

lgd_coloc = Legend(
  col_fun = coloc_genes, 
  title = "Colocalized Genes", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)

twas_genes = colorRamp2(c(0,4,8), c('white','darkslategray1','deepskyblue4'))

twas_matrix = as.matrix(data[c(6,5),])
data_twas = Heatmap(twas_matrix, name = 'TWAS Genes', col = twas_genes, show_row_dend = F, cluster_rows = F, cluster_columns = F, column_title = NULL,
                    show_column_names = T, column_split = cellnumber$Lineage, width = unit(ncol(twas_matrix)*5, 'mm'), height = unit(nrow(twas_matrix)*5, 'mm'), show_row_names = T, border = T, 
                    row_labels = c('Novel','Replicated'),
                    column_names_gp = gpar(fontsize = 6, col = c(rep('purple',7),rep('green',8),rep('blue',13),rep('orange',5))),
                    row_names_gp = gpar(fontsize=6),
                    show_heatmap_legend = F, 
                    cell_fun = function(j, i, x, y, w, h, fill) {
                      if(twas_matrix[i,j] > 0) {
                        grid::grid.text(sprintf("%.0f", twas_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 5))
                      }
                    })

lgd_twas = Legend(
  col_fun = twas_genes, 
  title = "TWAS Genes", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)

twas_loci = colorRamp2(c(0,3,6), c('white','darkgoldenrod1','darkorange4'))
twas_novel = as.matrix(t(data['TWAS_novel_loci',]))

data_twas_loci = Heatmap(twas_novel, name = 'Novel Loci (TWAS)', col = twas_loci, show_row_dend = F, cluster_rows = F, cluster_columns = F, column_title = NULL,
                         show_column_names = T, column_split = cellnumber$Lineage, width = unit(ncol(twas_novel)*5, 'mm'), height = unit(5, 'mm'), show_row_names = T, border = T, row_labels = 'Novel Loci',
                         column_names_gp = gpar(fontsize = 6, col = c(rep('purple',7),rep('green',8),rep('blue',13),rep('orange',5))),
                         row_names_gp = gpar(fontsize=6),
                         show_heatmap_legend = F,
                         cell_fun = function(j, i, x, y, w, h, fill) {
                           if(twas_novel[i,j] > 0) {
                             grid::grid.text(sprintf("%.0f", twas_novel[i,j]), x, y, gp = grid::gpar(fontsize = 5))
                           }
                         })

lgd_twas_loci = Legend(
  col_fun = twas_loci, 
  title = "Novel Loci (TWAS)", 
  direction = "horizontal",
  title_position = "topleft",
  labels_gp = gpar(fontsize = 5),
  title_gp = gpar(fontsize=5),
  grid_width = unit(3, "mm"),
  grid_height = unit(2, "mm")# You can also set title position within the legend box itself
)

pd =packLegend(lgd_sig_loci, lgd_twas_loci, lgd_twas, lgd_coloc, direction = 'horizontal')




c = draw(data_sig_loci %v% data_twas_loci %v% data_twas %v% data_coloc)
w = ComplexHeatmap:::width(c)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(c)
h = convertY(h, "inch", valueOnly = TRUE)
h2 = ComplexHeatmap:::height(pd)
h2 = convertY(h2, "inch", valueOnly = TRUE)
h = h+h2
c(w, h)
png('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_4/Figure_4c.png', w = w, h = (h+.2), units = 'in', res = 450)
draw(data_sig_loci %v% data_twas_loci %v% data_twas %v% data_coloc, ht_gap = unit(.5,'mm'))
draw(
  pd, 
  x = unit(0.05, "npc"),  # Adjust x coordinate (e.g., 5% from left)
  y = unit(.98, "npc"),  # Adjust y coordinate (e.g., 95% from bottom/5% from top)
  just = c("left", "top") # Align the left-top corner of the legend with the x,y coords
)
dev.off()






