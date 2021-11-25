import hail as hl
import argparse
import sys
import subprocess
import pkg_resources

required = {'matplotlib', 'pandas'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing])

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def pca_plot(pcs_data: str = None, ref_data: str = None, ref_info_data: str = None, include_ref: bool = None,
             cohort: str = None, out_dir: str = None, prob: float = None, phenotypes_file: str = None, pheno: str = None):

    ref = pd.read_table(ref_data, header=0, sep='\t', compression='gzip')
    ref_info = pd.read_table(ref_info_data, header=0, sep='\t')
    ref_info.rename(columns={'Sample': 's'}, inplace=True)
    ref_update = pd.merge(ref, ref_info, how='left', on=['s'])

    pcs = pd.read_table(pcs_data, header=0, sep='\t')

    pheno_file = pd.read_table(phenotypes_file, header=0, sep='\t')
    pheno_file.rename(columns={'ID': 's'}, inplace=True)

    if cohort:
        pcs_df = pd.merge(pcs, pheno_file, how='left', on=['s'])
        pcs_df = pcs_df[pcs_df.Cohort == cohort]
        # Only retain samples with phenotype
        pcs_df = pcs_df[pcs_df[pheno].notna()]
    else:
        pcs_df = pcs
        pcs_df = pd.merge(pcs_df, pheno_file, how='left', on=['s'])
        # Only retain samples with phenotype
        pcs_df = pcs_df[pcs_df[pheno].notna()]

    cbPalette = {'AFR': "#999999", 'AMR': "#E69F00", 'EAS': "#56B4E9", 'EUR': "#009E73",
             'oth': "#F0E442", 'SAS': "#0072B2"}

    if include_ref is True:
        cbPalette = cbPalette
        if cohort:
            plt_title = f'{cohort} projected against 1KG (p > {prob})'
            outfile = f'{cohort}.1KG.PC1-PC10.prob_{prob}.png'
        else:
            plt_title = f'All cohorts projected against 1KG (p > {prob})'
            outfile = f'all.cohorts.1KG.PC1-PC10.prob_{prob}.png'
    else:
        plt_title = f'{cohort} (p > {prob})'
        outfile = f'{cohort}.PC1-PC10.prob_{prob}.png'
        cbPalette = {k: cbPalette[k] for k in list(pcs_df['pop'].value_counts().index)}

    fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(25, 25))

    for i in range(1,10,2):
        if i == 1 or i == 3:
            axisx = 0
        elif i == 5 or i == 7:
            axisx = 1
        else:
            axisx = 2

        xpc = f'PC{i}'
        ypc = f'PC{i+1}'

        if i == 1 or i == 5 or i == 9:
            axisy = 0
        else:
            axisy = 1

        xpc = f'PC{i}'
        ypc = f'PC{i+1}'

        if include_ref:
            axs[axisx, axisy].scatter(ref_update[xpc], ref_update[ypc], c=ref_update['SuperPop'].map(cbPalette),
                                      s=5, alpha=0.1)
        axs[axisx, axisy].scatter(pcs_df[xpc], pcs_df[ypc], c=pcs_df['pop'].map(cbPalette), s=5, alpha=1)
        axs[axisx, axisy].set_xlabel(xlabel=xpc, fontsize=15)
        axs[axisx, axisy].set_ylabel(ylabel=ypc, fontsize=15)
        axs[axisx, axisy].tick_params(axis='both', which='major', labelsize=12)
    handles = []
    for key in cbPalette:
        # manually define a new patch
        data_key = Line2D([0], [0], marker='o', color='w', label=key,
                        markerfacecolor=cbPalette[key], markersize=10)
        handles.append(data_key)

    fig.legend(handles=handles, title='Populations', bbox_to_anchor=(0.7, 0.2), loc='lower left', frameon=False)
    fig.delaxes(axs[2][1])
    fig.suptitle(plt_title, fontsize=25)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.close()
    fig.savefig(f'/tmp/{outfile}', dpi=300, facecolor='w')
    hl.hadoop_copy(f'file:///tmp/{outfile}', f'{out_dir}plots/{outfile}')


def main(args):
    rf_thresh = [0.5, 0.8]

    dir = args.dirname
    ref_scores = args.dirname + 'qc_step_2/1000G_scores.txt.bgz'
    ref_info = args.dirname + 'qc_step_3/G1000_pops.txt'
    outdirectory = args.dirname + 'qc_step_3/'
    for threshold in rf_thresh:
        cohorts = args.cohort

        # PCA plot for all cohorts projected against 1KG
        pca_plot(pcs_data = '{}qc_step_3/pca_sup_pops_{}_probs.txt'.format(dir, threshold), ref_data=ref_scores,
        ref_info_data=ref_info, include_ref=True, out_dir=outdirectory, prob=threshold,
        phenotypes_file=args.phenotype_file, pheno=args.phenotype)

        for i in cohorts:
            # PCA for each cohort separately projected against 1KG
            pca_plot(pcs_data = '{}qc_step_3/pca_sup_pops_{}_probs.txt'.format(dir, threshold), ref_data=ref_scores,
            ref_info_data=ref_info, include_ref=True, cohort = i, out_dir=outdirectory, prob=threshold,
            phenotypes_file=args.phenotype_file, pheno=args.phenotype)

            # # PCA for each cohort separately
            pca_plot(pcs_data = '{}qc_step_3/pca_sup_pops_{}_probs.txt'.format(dir, threshold), ref_data=ref_scores,
            ref_info_data=ref_info, include_ref=False, cohort = i, out_dir=outdirectory, prob=threshold,
            phenotypes_file=args.phenotype_file, pheno=args.phenotype)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--ref_scores', type=str, required=True)
    # parser.add_argument('--ref_info', type=str, required=True)
    # parser.add_argument('--pca_dir', type=str, required=True)
    parser.add_argument('--phenotype_file', type=str, required=True)
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--phenotype', default='Analysis_C2')
    parser.add_argument('--cohort', action='append', help="Name of the cohort you want to generate PCA plots for")
    # parser.add_argument('--out_dir', type=str, required=True)

    args = parser.parse_args()
    main(args)

# TODO: Plots are upside down
# hailctl dataproc submit qc pca.py
# --dirname gs://lindo/qc_tests/
# --phenotype_file gs://dsge-covid19-data/phenotypes_main_2021_04_11_wo_mismatch.tsv
# --cohort italy --cohort belgium --cohort egypt --cohort iran --cohort sweden --cohort germany
