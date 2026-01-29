setwd('/data/Choi_lung/TTL/Seurat/Final_Batch_Lift/Data_vireo/same_ref')

library(Seurat)
library(Rsamtools)
library(GenomicRanges)
library(DropletQC)
library(ggplot2)
library(patchwork)
library(dplyr)

#loading and integrating data from scrublet, demuxlet, vireo, dropletqc

load_d = function(d,s,n,nf,dm,vr,dm2) {
  run.data = Read10X(data.dir=d)
  run = CreateSeuratObject(counts= run.data, project = n, min.cells=3, min.features =200)
  
  run[["percent.mt"]] = PercentageFeatureSet(run, pattern = "^MT-")
  run[['percent.rb']] = PercentageFeatureSet(run, pattern = '^RP[SL]')
  
  doublets = read.table(s, header= F, row.names=1)
  colnames(doublets) = c('Doublet_score', 'Is_doublet')
  run = AddMetaData(run, doublets)
  
  nf2_umi_file = read.csv(nf, header = TRUE, row.names=1)
  colnames(nf2_umi_file) = c('nuclear_fraction')
  run = AddMetaData(run, nf2_umi_file)
  data = data.frame(nf=run$nuclear_fraction, umi=run$nCount_RNA)
  droplet_test = identify_empty_drops(nf_umi = data, include_plot = FALSE)
 
  run[['cell_status']] =  droplet_test$cell_status
  
  demux = read.table(dm,row.names = 1,header=TRUE)
  
  run = AddMetaData(run, demux)
  
  vireo = read.table(vr,row.names = 1,header=TRUE)
  
  run = AddMetaData(run, vireo)
  
  columns.to.remove <- c('RD.TOTL',	'RD.PASS',	'RD.UNIQ',	'N.SNP',	'SNG.LLK1',	'SNG.2ND',	'SNG.LLK2',	'SNG.LLK0',	'DBL.1ST',	'DBL.2ND',	'ALPHA',	'LLK12',	'LLK1',	'LLK2',	'LLK10',	'LLK20',	'LLK00', 
                         'n_vars',  'best_doublet')
  for(i in columns.to.remove) {
    run[[i]] <- NULL
  }
  
  run$BEST = sub("^DBL.*", "DBL", run$BEST)
  run$BEST = sub('^SNG.*', 'SNG', run$BEST)
  run$BEST = sub('^AMB.*', 'AMB', run$BEST)
  run$demux.doublet.call = run$BEST
  run$BEST = NULL
  run$Sample = run$SNG.1ST
  run$SNG.1ST = NULL
  
  demux2 = read.table(dm2,row.names = 1,header=TRUE)
  
  run = AddMetaData(run, demux2)
  
  run$BEST = sub("^DBL.*", "DBL", run$BEST)
  run$BEST = sub('^SNG.*', 'SNG', run$BEST)
  run$BEST = sub('^AMB.*', 'AMB', run$BEST)
  run$demux.doublet.call_2 = run$BEST
  run$BEST = NULL
  run$Sample_2 = run$SNG.1ST
  run$SNG.1ST = NULL
  
  return(run)
} 

nci1257=load_d('/data/Choi_lung/scRNA_eQTL/nci1257/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci1257_scrublet.tsv','nci1257', '/data/Choi_lung/TTL/Seurat/DropletQC/nci1257_nf.csv',
               '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci1257.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci1257_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci1257.tsv')

nci1257_II=load_d('/data/Choi_lung/scRNA_eQTL/nci1257_II/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci1257_II_scrublet.tsv','nci1257_II','/data/Choi_lung/TTL/Seurat/DropletQC/nci1257_II_nf.csv',
                 '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci1257_II.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci1257_II_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci1257_II.tsv')

nci9_16=load_d('/data/Choi_lung/scRNA_eQTL/nci9_16/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci9_16_scrublet.tsv','nci9_16','/data/Choi_lung/TTL/Seurat/DropletQC/nci9_16_nf.csv',
               '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci9_16.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci9_16_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci9_16.tsv')

nci17_22=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci17_22_scrublet.tsv','nci17_22','/data/Choi_lung/TTL/Seurat/DropletQC/nci17_22_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci17_22.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci17_22_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci17_22.tsv')

nci23_29=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci23_29_scrublet.tsv','nci23_29','/data/Choi_lung/TTL/Seurat/DropletQC/nci23_29_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci23_29.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci23_29_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci23_29.tsv')

nci30_35=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci30_35_scrublet.tsv','nci30_35','/data/Choi_lung/TTL/Seurat/DropletQC/nci30_35_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci30_35.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci30_35_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci30_35.tsv')

nci36_41=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci36_41_scrublet.tsv','nci36_41','/data/Choi_lung/TTL/Seurat/DropletQC/nci36_41_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci36_41.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci36_41_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci36_41.tsv')

nci42_47=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci42_47_scrublet.tsv','nci42_47','/data/Choi_lung/TTL/Seurat/DropletQC/nci42_47_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci42_47.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci42_47_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci42_47.tsv')

nci48_54=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci48_54_scrublet.tsv','nci48_54','/data/Choi_lung/TTL/Seurat/DropletQC/nci48_54_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci48_54.tsv','/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci48_54_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci48_54.tsv')

nci56_61=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci56_61_scrublet.tsv','nci56_61','/data/Choi_lung/TTL/Seurat/DropletQC/nci56_61_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci56_61.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci56_61_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci56_61.tsv')

nci62_68=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci62_68_scrublet.tsv','nci62_68','/data/Choi_lung/TTL/Seurat/DropletQC/nci62_68_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci62_68.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci62_68_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci62_68.tsv')

nci69_74=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci69_74_scrublet.tsv','nci69_74','/data/Choi_lung/TTL/Seurat/DropletQC/nci69_74_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci69_74.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci69_74_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci69_74.tsv')

nci75_80=load_d('/data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci75_80_scrublet.tsv','nci75_80','/data/Choi_lung/TTL/Seurat/DropletQC/nci75_80_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci75_80.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci75_80_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci75_80.tsv')

nci81_86=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci81_86_scrublet.tsv','nci81_86','/data/Choi_lung/TTL/Seurat/DropletQC/nci81_86_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci81_86.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci81_86_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci81_86.tsv')

nci87_92=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci87_92_scrublet.tsv','nci87_92','/data/Choi_lung/TTL/Seurat/DropletQC/nci87_92_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci87_92.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci87_92_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci87_92.tsv')

nci93_98=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci93_98_scrublet.tsv','nci93_98','/data/Choi_lung/TTL/Seurat/DropletQC/nci93_98_nf.csv',
                '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci93_98.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci93_98_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci93_98.tsv')

nci99_104=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci99_104_scrublet.tsv','nci99_104','/data/Choi_lung/TTL/Seurat/DropletQC/nci99_104_nf.csv',
                 '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci99_104.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci99_104_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci99_104.tsv')

nci105_110=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci105_110_scrublet.tsv','nci105_110','/data/Choi_lung/TTL/Seurat/DropletQC/nci105_110_nf.csv',
                  '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci105_110.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci105_110_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci105_110.tsv')

nci111_116=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci111_116_scrublet.tsv','nci111_116','/data/Choi_lung/TTL/Seurat/DropletQC/nci111_116_nf.csv',
                  '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci111_116.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci111_116_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci111_116.tsv')

nci117_123 =load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci117_123_II_scrublet.tsv','nci117_123','/data/Choi_lung/TTL/Seurat/DropletQC/nci117_123_II_nf.csv',
                   '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci117_123_II.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci117_123_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci117_123.tsv')

nci124_129=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci124_129_scrublet.tsv','nci124_129','/data/Choi_lung/TTL/Seurat/DropletQC/nci124_129_nf.csv',
                  '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci124_129.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci124_129_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci124_129.tsv')

nci130_135=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci130_135_scrublet.tsv','nci130_135','/data/Choi_lung/TTL/Seurat/DropletQC/nci130_135_nf.csv',
                  '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci130_135.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci130_135_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci130_135.tsv')

nci136_140=load_d('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/filtered_feature_bc_matrix/','/data/Choi_lung/TTL/Seurat/Scrublet/nci136_140_scrublet.tsv','nci136_140','/data/Choi_lung/TTL/Seurat/DropletQC/nci136_140_nf.csv',
                  '/data/Choi_lung/TTL/Demuxlet/Lift_work/nci136_140.tsv', '/data/Choi_lung/TTL/vireoSNP/Final/vireo_nci136_140_wSNP/donor_ids.tsv', '/data/Choi_lung/TTL/Genotype/liftover_chr_files/demuxlet__test/Output/nci136_140.tsv')

