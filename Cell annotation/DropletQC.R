setwd('/data/Choi_lung/TTL/Seurat/DropletQC')

library(Rsamtools)
library(GenomicRanges)
library(DropletQC)
library(ggplot2)
library(patchwork)
library(dplyr)

#Create a reference to the BAM file

nf=function(x,y,z){
  
  nf2 = nuclear_fraction_tags(bam=x,
                              barcodes= y,
                              cores = 1,
                              verbose= FALSE)
  
  write.csv(nf2, z)
}

nf('/data/Choi_lung/scRNA_eQTL/nci1257/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/nci1257/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci1257_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/nci1257_II/outs/possorted_genome_bam.bam',
  '/data/Choi_lung/scRNA_eQTL/nci1257_II/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci1257_II_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/nci9_16/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/nci9_16/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', './nci9_16_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci17_22_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci23_29_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci30_35_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci36_41_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci42_47_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci48_54_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci56_61_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci62_68_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci69_74_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci75_80_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci81_86_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci87_92_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci93_98_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci99_104_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci105_110_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci111_116_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci117_123_II_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci124_129_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci130_135_nf.csv')

nf('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/possorted_genome_bam.bam',
   '/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/filtered_feature_bc_matrix/barcodes.tsv.gz','./nci136_140_nf.csv')


