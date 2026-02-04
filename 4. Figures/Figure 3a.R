#Author: Thong Luong
#Date: April 28th 2025
#Modfiied: December 31 2025

library(AnnotationHub)
library(rtracklayer)
library(annotatr)
library(ggplot2)
library(qvalue)
library(rstatix)

setwd('/data/Choi_lung/TTL/Figs_for_manuscript/Enrichment/LD_prune_strict')

#get LD prune genotype
ld_prune = read.csv('/data/Choi_lung/TTL/tensor/Genotype/PLINK/chr_pos/LD_prune_strict/prune_geno.bim', header = F, sep = '\t')
eqtl_snps = read.table('/data/Choi_lung/TTL/SNP_cCRE/eqtl_snps.bed', sep = '\t')
eqtl_snps = subset(eqtl_snps, V4 %in% ld_prune$V2)

celltypes = readRDS('/data/Choi_lung/TTL/files_for_JC_2410/ct.rds')

egenes = data.frame()

for (c in celltypes){
  cell = read.csv(paste0('/data/Choi_lung/TTL/tensor/Output_Sum_Final_chr_pos/output/',c,'/top_qtl_results_all_FDR0.05.txt'), sep = '\t')
  egenes = rbind(egenes,cell)
}

snps = unique(egenes$variant_id)
esnp = subset(eqtl_snps, V4 %in% snps)
not_esnp = subset(eqtl_snps, !(V4 %in% snps))
write.table(esnp, './Bed_files/full_esnp.bed',sep = '\t', col.names = F, row.names = F, quote= F)
write.table(not_esnp,'./Bed_files/full_not_esnp.bed',sep = '\t', col.names = F, row.names = F, quote= F)

ah = AnnotationHub()

E096 <- query(ah , c("EpigenomeRoadMap", "E096", "segmentations"))
chainfiles = query(ah, c('hg38','hg19','chainfile'))
chain = chainfiles[['AH14150']]
gtf = E096[["AH46949"]]
grch38 = liftOver(gtf,chain)
genome(grch38) = 'hg38'
grch38 = unlist(grch38)
type = unique(grch38$name)

get_overlap = function(c){
  erange = read_regions(con = paste0('./Bed_files/',c,'_esnp.bed') , genome = 'hg38', format = 'bed')
  eintersected = as.data.frame(annotate_regions(erange, grch38, ignore.strand = T, quiet = F))
  eintersected = eintersected[!duplicated(eintersected$name),]
  
  not_erange = read_regions(con = paste0('./Bed_files/',c,'_not_esnp.bed') , genome = 'hg38', format = 'bed')
  not_eintersected = as.data.frame(annotate_regions(not_erange, grch38, ignore.strand = T, quiet = F))
  not_eintersected = not_eintersected[!duplicated(not_eintersected$name),]
  
  df = data.frame(matrix(vector(), 15, 2,
                         dimnames=list(c(type), c('esnp','not_esnp'))),
                  stringsAsFactors=F)
  
  
  df[unique(eintersected$annot.name),'esnp'] = table(eintersected$annot.name)
  df[unique(not_eintersected$annot.name),'not_esnp'] = table(not_eintersected$annot.name)
  df$sig_background = nrow(as.data.frame(erange))
  df$background = nrow(as.data.frame(not_erange))
  write.table(df, paste0(c,'_enrichment.tsv'), sep = '\t', row.names = T, col.names = T, quote =F)
}

get_overlap('full')

test = read.csv('full_enrichment.tsv', sep = '\t')
test[is.na(test)] <- 0
states = rownames(test)

Fish_enrich = data.frame(matrix(vector(), 15, 2,
                                dimnames=list(c(rownames(test)), c('Fisher_pval','Odds_ratio'))),
                         stringsAsFactors=F)
Fish_enrich$State = rownames(Fish_enrich)

data = data.frame()
df = read.csv(paste0('full_enrichment.tsv'),sep= '\t')
df[is.na(df)] <- 0
for (j in states){
    dat <- data.frame(
      "Sig" = c(df[j,'esnp'], df[j,'sig_background']),
      "Not_sig" = c(df[j,'not_esnp'], df[j,'background']),
      row.names = c("Overlap", "Background"),
      stringsAsFactors = FALSE
    )
    Fish_enrich[j,'Fisher_pval'] = fisher.test(dat)$p.value
    Fish_enrich[j,'Odds_ratio'] = fisher.test(dat)$estimate[['odds ratio']]
}
rownames(Fish_enrich) = NULL
data = Fish_enrich
data$FDR = p.adjust(data$Fisher_pval, method = 'BH')

data$log_odds = log(data$Odds_ratio)
data[data==0] = min(data[data$FDR>0,"FDR"])
data$log_FDR = -log10(data$FDR)
data[12,'log_odds'] = 0
data[13,'log_odds'] = 0
data$Significance = ifelse(data$log_FDR >= 1.3, 'Yes', 'No')

order = c('Active TSS','Flanking Active TSS',"Transcr. at gene 5' and 3'",'Strong transcription','Weak transcription',
          'Genic enhancers','Enhancers','ZNF genes & repeats','Heterochromatin','Bivalent/Poised TSS','Flanking Bivalent TSS/Enh',
          'Bivalent Enhancer', 'Repressed PolyComb', 'Weak Repressed PolyComb','Quiescent/Low')
order = rev(order)
data$State = factor(data$State, levels = order)


t = ggplot()+
  geom_point(data, mapping = aes(x= log_odds, y=State, color = log_odds, size=log_FDR, shape = Significance)) +
  scale_colour_gradient(low="blue", high="red", name = "Log2(Odds Ratio)") +
  scale_shape_manual(values=c(1,16), name = 'Significant') +
  scale_size(range=c(.5,2.5),breaks=c(0,0.5,1.3,5,10),labels=c(">=0",'>=0.5',">=1.3",'>=5','>=10'),guide="legend", name = '-Log10(FDR)') + 
  labs(y = "Chromatin State", x = "Log Odds", fill = "Odds Ratio") +
  coord_cartesian(xlim = c(-8, 8)) +
  ggtitle('') +
  theme(text = element_text(family = 'Arial'),
        axis.text=element_text(size=5),
        axis.title.x = element_text(color = "black", size = 6, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 6, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        legend.key.size = unit(2, 'mm'),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(3, "pt"), 
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        panel.background = element_rect(fill = 'white', colour = 'black', linetype = 'solid'),
        legend.key = element_rect(color = NA, fill = NA))

plot_save_dir = '/data/Choi_lung/TTL/Figs_for_manuscript/Final_Figures/Fig_3/'

ggsave(filename = file.path(plot_save_dir, 'Figure_3a.png'), plot = t, width = 80, height = 65, units = "mm", dpi = 450, device = 'png')


  
