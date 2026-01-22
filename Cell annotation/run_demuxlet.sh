#!/bin/sh

cd /data/Choi_lung/TTL/Demuxlet/Lift_work/

/data/lib14/Software/bin/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci1257_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci1257_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/nci1257/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci1257

/data/lib14/Software/bin/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci1257_II_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci1257_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/nci1257_II/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci1257_II


/data/lib14/Software/bin/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci9_16_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci9_16_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/nci9_16/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci9_16

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci17_22_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci17_22_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci17_22

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci23_29_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci23_29_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci23_29

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci30_35_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci30_35_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci30_35

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci36_41_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci36_41_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci36_41

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci42_47_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci42_47_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci42_47

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci48_54_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci48_54_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci48_54

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci56_61_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci56_61_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci56_61

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci62_68_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci62_68_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci62_68

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci69_74_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci69_74_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci69_74

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci75_80_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci75_80_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci75_80

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci81_86_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci81_86_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci81_86


/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci87_92_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci87_92_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci87_92

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci93_98_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci93_98_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci93_98

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci99_104_chr.bam\
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci99_104_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci99_104

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci105_110_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci105_110_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci105_110

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci111_116_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci111_116_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci111_116

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci117_123_I_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci117_123_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-I/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci117_123_I

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci117_123_II_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci117_123_sorted.vcf  \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci117_123_II

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci124_129_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci124_129_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci124_129

/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci130_135_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci130_135_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci130_135
                                                         
/data/Choi_lung/scRNA_eQTL/demuxlet/bin/centos/demuxlet  --sam /data/Choi_lung/TTL/Demuxlet/Reheader_Bam/nci136_140_chr.bam \
                                                         --vcf /data/Choi_lung/TTL/Demuxlet/0.9/nci136_140_sorted.vcf \
                                                         --group-list /data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/filtered_feature_bc_matrix/barcodes.tsv \
                                                         --alpha 0 --alpha 0.5 \
                                                         --field GT \
                                                         --out nci136_140




