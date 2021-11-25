import hail as hl
import argparse
import pandas as pd

hl.init(default_reference="GRCh38")


def run_pca(mt: hl.MatrixTable, out_prefix: str, overwrite: bool = True):
    """
    Run PCA on a dataset
    :param mt: dataset to run PCA on
    :param out_prefix: directory and filename prefix for where to put PCA output
    :return:
    """
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    pca_mt = mt.annotate_rows(pca_af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)

    pca_scores.write(out_prefix + 'scores.ht', overwrite)
    pca_scores = hl.read_table(out_prefix + 'scores.ht')
    pca_scores = pca_scores.transmute(**{f'PC{i}': pca_scores.scores[i - 1] for i in range(1, 21)})
    pca_scores.export(out_prefix + 'scores.txt')  # individual-level PCs
    pca_loadings.write(out_prefix + 'loadings.ht', overwrite)  # PCA loadings


pheno = hl.import_table('gs://dsge-covid19-data/phenotypes_main_2021_04_11_wo_mismatch.tsv', impute=True)
pheno = pheno.key_by('ID')

mt = hl.read_matrix_table('gs://dsge-covid19-data/15042021/qc_step_5/to_imputation/cov_15042021.qc.merged_controls.mt')
mt = mt.annotate_cols(pheno=pheno[mt.s])

print(mt.aggregate_cols(hl.agg.counter(mt.pheno.Cohort)))

out_prefix = 'gs://dsge-covid19-data/15042021/qc_step_2-pca_cas_con/'

cohortList = ['belgium', 'egypt', 'germany', 'iran', 'italy', 'sweden']

print("\n all: " + str(mt.count()))
run_pca(mt, out_prefix + 'all_')

for i in cohortList:
    mtco = mt.filter_cols(mt.pheno.Cohort == i, keep=True)

    print("\n {}: ".format(i) + str(mtco.count()))

    run_pca(mtco, out_prefix + '{}_'.format(i))
