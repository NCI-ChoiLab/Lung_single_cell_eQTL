import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt


#d is file directory, e is expected_doublet_rate, h is histogram name, u is UMAP cluster, z is scrublet tsv
def scb(d,e,h,u,z):
    data = sc.read_10x_mtx(d, cache=False)
    data.var_names_make_unique()
    scrub = scr.Scrublet(data.X, expected_doublet_rate=e)
    data.obs['doublet_scores'], data.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2,
                                                                                      min_cells=3,
                                                                                      min_gene_variability_pctl=85,
                                                                                      n_prin_comps=30)
    scrub.plot_histogram()
    plt.savefig(h)
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)
    plt.savefig(u)

    data.obs['predicted_doublets'].value_counts()
    pd.DataFrame(data.obs).to_csv(z, sep='\t', header=False)

#adjust for expected doublet rates based on estimated number of cells
scb('/data/Choi_lung/scRNA_eQTL/nci1257/outs/filtered_feature_bc_matrix/',0.09,'nci1257.tiff','nci1257_umap.tiff','nci1257_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/nci1257_II/outs/filtered_feature_bc_matrix/',0.08,'nci1257_2.tiff','nci1257_2_umap.tiff','nci1257_II_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/nci9_16/outs/filtered_feature_bc_matrix/',0.14,'nci9_16_his.tiff','nci9_16_umap.tiff','nci9_16_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI17_22/outs/filtered_feature_bc_matrix/',0.25,'nci17_22_his.tiff','nci17_22_umap.tiff','nci17_22_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI23_29/outs/filtered_feature_bc_matrix/',0.41,'nci23_29_his.tiff','nci23_29_umap.tiff','nci23_29_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI30_35/outs/filtered_feature_bc_matrix/',0.18,'nci30_35_his.tiff','nci30_35_umap.tiff','nci30_35_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI36_41/outs/filtered_feature_bc_matrix/',0.22,'nci36_41_his.tiff','nci36_41_umap.tiff','nci36_41_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI42_47/outs/filtered_feature_bc_matrix/',0.26,'nci42_47_his.tiff','nci42_47_umap.tiff','nci42_47_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI48_54/outs/filtered_feature_bc_matrix/',0.23,'nci48_54_his.tiff','nci48_54_umap.tiff','nci48_54_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI56_61/outs/filtered_feature_bc_matrix/',0.21,'nci56_61_his.tiff','nci56_61_umap.tiff','nci56_61_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI62_68/outs/filtered_feature_bc_matrix/',0.25,'nci62_68_his.tiff','nci62_68_umap.tiff','nci162_68_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI69_74/outs/filtered_feature_bc_matrix/',0.23,'nci69_74_his.tiff','nci69_74_umap.tiff','nci69_74_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/CS033335/NCI75_80/outs/filtered_feature_bc_matrix/',0.27,'nci75_80_his.tiff','nci75_80_umap.tiff','nci75_80_scrublet.tsv')

scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI81-86/outs/filtered_feature_bc_matrix/',0.26,'nci81_86_his.tiff','nci81_86_umap.tiff','nci81_86_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI87-92/outs/filtered_feature_bc_matrix/',0.22,'nci87_92_his.tiff','nci87_92_umap.tiff','nci87_92_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI93-98/outs/filtered_feature_bc_matrix/',0.24,'nci93_98_his.tiff','nci93_98_umap.tiff','nci93_98_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI99-104/outs/filtered_feature_bc_matrix/',0.25,'nci99_104_his.tiff','nci99_104_umap.tiff','nci99_104_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI105-110/outs/filtered_feature_bc_matrix/',0.20,'nci105_110_his.tiff','nci105_110_umap.tiff','nci105_110_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI111-116/outs/filtered_feature_bc_matrix/',0.26,'nci111_116_his.tiff','nci111_116_umap.tiff','nci111_116_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI117-123-II/outs/filtered_feature_bc_matrix/',0.17,'nci117_123_II_his.tiff','nci117_123_II_umap.tiff','nci117_123_II_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI124-129/outs/filtered_feature_bc_matrix/',0.15,'nci124_129_his.tiff','nci124_129_umap.tiff','nci124_129_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI130-135/outs/filtered_feature_bc_matrix/',0.23,'nci130_135_his.tiff','nci130_135_umap.tiff','nci130_135_scrublet.tsv')
scb('/data/Choi_lung/scRNA_eQTL/230131_A00423_0226_AH2JGKDSX5_ttl/run_count_NCI136-140/outs/filtered_feature_bc_matrix/',0.18,'nci136_140_his.tiff','nci136_140_umap.tiff','nci136_140_scrublet.tsv')







