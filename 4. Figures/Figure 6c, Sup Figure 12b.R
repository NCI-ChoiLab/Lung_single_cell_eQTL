library(dplyr)
library(ggplot2)
library(stringr)
library(stringi)
library(reshape2)
library(Matrix)
library(Hmisc)
library(wesanderson)
library(magrittr)
library(rtracklayer)
library(geomtextpath)
library(snpStats)

fam <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.fam"
bim <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bim"
bed <- "/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/final.bed"
sample <- read.plink(bed, bim, fam)
genotype_mtx <- sample$genotypes
genotypes<- genotype_mtx@.Data
snp_info <- sample$map

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/')

celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

snp_info = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/snps.rds')
rownames(snp_info) = paste(snp_info$V1,snp_info$V4,snp_info$V6,snp_info$V5,sep='_')

expr_list = list()
for (c in celltypes){
  exp = read.csv(paste0('/data/Choi_lung/TTL/tensor/Phenotype/Sum_Final/',c,'.bed'), sep = '\t')
  expr_list[[c]] = exp
}

eQTL_plot_pub <- function(celltype,rs,transcript,gene_name,p,b){
  exprs <- expr_list[[celltype]][,-c(1:4)]
  rownames(exprs) <- expr_list[[celltype]]$gene_id
  samples <- colnames(exprs)
  
  datainput2 <- data.frame(snp = genotypes[samples,rs],
                           expression = as.numeric(exprs[transcript,]))
  
  datainput2$snp_ra <- "0|0"
  datainput2$snp_ra[which(datainput2$snp == "01")] <- "1|1"
  datainput2$snp_ra[which(datainput2$snp == "02")] <- "0|1"
  df <- as.vector(table(datainput2$snp_ra))
  datainput2$snp_ra <- paste0("0|0\n", "n=", df[1])
  datainput2$snp_ra[which(datainput2$snp == "01")] <- paste0("1|1\n", "n=", df[3])
  datainput2$snp_ra[which(datainput2$snp == "02")] <- paste0("0|1\n", "n=", df[2])
  p2 <- ggplot(datainput2, aes(x=snp_ra, y=expression, color=snp_ra)) +
    geom_violin(trim=FALSE,linewidth = 0.25)+
    #geom_point(position=position_jitterdodge(),aes(group=snp_ra)) +
    geom_boxplot(aes(middle=median(expression)),
                 width=0.1, linewidth = 0.25)+
    geom_smooth(mapping = aes(x = snp_ra, y = expression, group = 1),formula = y~x, color = "gray",
                method='lm', size = .5, se =TRUE,fill = alpha("gray", .5) ) +
    labs(title="",x=snp_info[rs,'V2'], y = paste0("Normalized expression of ", gene_name))+
    scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ggtitle(paste0(celltype, "\n",gene_name , ": ",snp_info[rs,'V2'], ": ", paste(rs, "b38", sep = "_"),'\n',
                   'pval = ',p))+
    theme(text = element_text(family = 'Arial'),
          plot.title = element_text(size = 6.5, color = 'black', face = 'plain'),
          axis.text.x = element_text(color = "black", size = 6,face = "plain"),
          axis.text.y = element_text(color = "black", size = 6,face = "plain"),  
          axis.title.x = element_text(color = "black", size = 6, face = "plain"),
          axis.title.y = element_text(color = "black", size = 6, face = "plain"),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(linewidth = .5, colour = "black", linetype=1),
          legend.position="none")

  return(p2)
}

p3 = eQTL_plot_pub('AT2','chr6_117464145_A_T','ENSG00000047936','ROS1','6.23e-08')
plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Sup_Figs/'
ggsave(filename = file.path(plot_save_dir, 'ROS1_box.png'), plot = p3, width = 60, height = 60, units = "mm", dpi = 300, device = 'png')

p4 = eQTL_plot_pub('AT2','chr10_112749531_T_C','ENSG00000148737','TCF7L2','1.43e-05')
plot_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_6/'
ggsave(filename = file.path(plot_save_dir, 'tcf7l2_box.png'), plot = p4, width = 60, height = 60, units = "mm", dpi = 300, device = 'png')
