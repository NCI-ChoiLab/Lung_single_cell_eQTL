setwd('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/QC_and_clustering/')

library(Matrix)
library(future)
library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(SeuratData)
library(harmony)
library(data.table)
library(ggplot2)
library(Azimuth)
library(reticulate)
use_condaenv('ld_clustering')

#size*1024^2, 3000*1024^2
options(future.globals.maxSize= 3145728000)

nci1257=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci1257.RDS')
nci1257_II=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci1257_II.RDS')
nci9_16=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci9_16.RDS')
nci17_22=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci17_22.RDS')
nci23_29=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci23_29.RDS')
nci30_35=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci30_35.RDS')
nci36_41=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci36_41.RDS')
nci42_47=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci42_47.RDS')
nci42_47 = subset(nci42_47, subset = Sample != 'SC845484_PC62791_G05')
nci48_54=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci48_54.RDS')
nci56_61=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci56_61.RDS')
nci62_68=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci62_68.RDS')
nci69_74=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci69_74.RDS')
nci75_80=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci75_80.RDS')
nci81_86=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci81_86.RDS')
nci87_92=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci87_92.RDS')
nci93_98=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci93_98.RDS')
nci99_104=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci99_104.RDS')
nci105_110=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci105_110.RDS')
nci111_116=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci111_116.RDS')
nci117_123=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci117_123.RDS')
nci124_129=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci124_129.RDS')
nci130_135=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci130_135.RDS')
nci136_140=readRDS('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref/nci136_140.RDS')


qc_subset = function(r){
  
  #scrublet and demuxlet = vireo
  r3 = subset(r, subset = cell_status == 'cell')
  r3 = subset(r3, subset = Is_doublet == 'True', invert = TRUE)
  r3 = subset(r3, subset = demux.doublet.call == 'DBL' & donor_id == 'doublet', invert = TRUE)
  r3 = subset(r3, subset = demux.doublet.call == 'DBL' & donor_id == 'unassigned', invert = TRUE)
  
  lf_f3 = log10(r3@meta.data$nFeature_RNA + 1)
  lmed_f3 = median(lf_f3)
  lmsd_f3 =mad(lf_f3, center=lmed_f3)
  nFeat_u3 = 10^(lmed_f3 + 2*lmsd_f3)
  nFeat_l3 = 10^(lmed_f3 - 2*lmsd_f3)
  
  r3 = subset(r3, subset = nFeature_RNA > nFeat_l3 & nFeature_RNA < nFeat_u3 & percent.mt < 20)
  
  return(r3)
}

nci1257_s = qc_subset(nci1257)
nci1257_II_s = qc_subset(nci1257_II)

nci9_16_s = qc_subset(nci9_16)
nci17_22_s = qc_subset(nci17_22)
nci23_29_s = qc_subset(nci23_29)
nci30_35_s = qc_subset(nci30_35)
nci36_41_s = qc_subset(nci36_41)
nci42_47_s = qc_subset(nci42_47)
nci48_54_s = qc_subset(nci48_54)
nci56_61_s = qc_subset(nci56_61)
nci62_68_s = qc_subset(nci62_68)
nci69_74_s = qc_subset(nci69_74)
nci75_80_s = qc_subset(nci75_80)

nci81_86_s = qc_subset(nci81_86)
nci87_92_s = qc_subset(nci87_92)
nci93_98_s = qc_subset(nci93_98)
nci99_104_s = qc_subset(nci99_104)
nci105_110_s = qc_subset(nci105_110)
nci111_116_s = qc_subset(nci111_116)
nci117_123_s = qc_subset(nci117_123)
nci124_129_s = qc_subset(nci124_129)
nci130_135_s = qc_subset(nci130_135)
nci136_140_s = qc_subset(nci136_140)

after =  merge(nci1257_s,y=c(nci1257_II_s,nci9_16_s,nci17_22_s, nci23_29_s,nci30_35_s,nci36_41_s,nci42_47_s,nci48_54_s,nci56_61_s,nci62_68_s,nci69_74_s,nci75_80_s,
                            nci81_86_s,nci87_92_s,nci93_98_s,nci99_104_s,nci105_110_s,nci111_116_s,nci117_123_s,nci124_129_s,nci130_135_s,nci136_140_s),project='nci1_140_after')

int = function(o, s){
  obj = JoinLayers(o)
  obj = RunAzimuth(obj, reference ='lungref')
  
  
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
  
  # run sctransform
  obj = SCTransform(obj, vst.flavor = "v2")
  obj = RunPCA(obj, npcs = 50, verbose = FALSE)
  
  # one-liner to run Integration
  obj = IntegrateLayers(object = obj, method = HarmonyIntegration,
                        orig.reduction = "pca", new.reduction = 'harmony',
                        assay = "SCT", verbose = FALSE)
  
  obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
  obj2 = FindClusters(obj, method = 'igraph', resolution = 0.1, algorithm = 4, verbose = F, cluster.name = 'harmony_clusters')
  obj2 = RunUMAP(obj2, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  saveRDS(obj2, file = s)
  return(obj2)
}


a = int(after, 'res0.1_UMAP.RDS')
rm(after)
rm(a)





