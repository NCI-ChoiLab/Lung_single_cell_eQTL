#!/bin/sh

source /data/Choi_lung/TTL/mambaforge/etc/profile.d/conda.sh

conda activate vireo_demulti

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/4FACS/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/4FACS/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci1257_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/4FACS_2/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/4FACS_2/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci1257_II_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/8FACS/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/8FACS/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci9_16_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci17_22_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci23_29_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci30_35_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci36_41_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci42_47_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci48_54_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci56_61_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci62_68_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci69_74_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci75_80_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci81_86_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci87_92_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci93_98_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci99_104_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci105_110_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci111_116_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci117_123_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci124_129_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci130_135_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

cellsnp-lite -s /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/possorted_genome_bam.bam -b /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/filtered_feature_bc_matrix/barcodes.tsv -O /data/Choi_lung/TTL/vireoSNP/Final/nci136_140_wSNP -R /data/Choi_lung/TTL/vireoSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz  -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

