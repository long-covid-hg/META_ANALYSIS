import hail as hl
from hail import hadoop_open
import argparse
import sys
import subprocess
import pkg_resources

required = {'matplotlib'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing])

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

hl.init(default_reference='GRCh38')

def mahalanobis(x=None, data=None, cov=None):
    """Compute the Mahalanobis Distance between each row of x and the data
    x       : vector or matrix of data with, say, p columns.
    data    : ndarray of the distribution from which Mahalanobis distance of each observation of x is to be computed.
    cov     : covariance matrix (p x p) of the distribution. If None, will be computed from data.
    return  : Mahalanobis distance in a numpy array
    """
    x_minus_mu = x - np.mean(data)
    if not cov:
        cov = np.cov(data.values.T)
    inv_covmat = np.linalg.inv(cov)
    x_minus_mu_arr = x_minus_mu.to_numpy() # convert the x_minus_mu df to a numpy array
    x_minus_mu_bm = hl.linalg.BlockMatrix.from_numpy(x_minus_mu_arr)
    inv_covmat_bm = hl.linalg.BlockMatrix.from_numpy(inv_covmat)
    left_term = x_minus_mu_bm @ inv_covmat
    mahal = left_term @ x_minus_mu_bm.T
    return mahal.diagonal().to_numpy().ravel() # to_numpy returns ndarray and .ravel() converts the ndarray to array

def mahalanobis_distance(mt, dirname, cohortname):
    """Compute the Mahalanobis Distance between each row of x and the data
    mt          : Hail MatrixTable
    dirname     : path to where the outputs will be saved
    cohortname  : the name of a cohort e.g. Sweden, ita etc.
    return      : filtered Hail MatrixTable
    """

    # PART 1
    gnomad_ht = hl.read_table('gs://dsge-covid19-data/gnomad/gnomad_nfe_AF_100p.ht')
    #print("The are {} SNPs in the gnomAD".format(gnomad_ht.count()))

    mt_input = hl.variant_qc(mt)
    ht = mt_input.rows()
    #print("The are {} SNPs in the data".format(ht.count()))

    # this will join the the two tables by locus, so there'll be SNPs with the same locus but different alleles (we remove them later on)
    table_joined = ht.key_by('locus').join(gnomad_ht.key_by('locus'))
    #print("The are {} SNPs common between the data and gnomAD".format(table_joined.count()))

    # remove the key so we can be able to select only a few columns
    # because both the gnomad and our data have the 'locus' and 'allele' cols, gnomad will be recoded as locus_1 and alleles_1
    table_result = table_joined.key_by()
    table_result = table_result.select(table_result.locus, data_allele=table_result.alleles, gnomad_alleles=table_result.alleles_1,
                                       data_rsid=table_result.rsid, data_AF=table_result.variant_qc.AF[1], gnomad_AF=table_result.nfe_AF)

    #flip = hl.case().when(table_result.data_allele[1] == table_result.gnomad_alleles[0], True).when(mt.allele == mt.alleles[1], False).or_missing()
    table_result = table_result.annotate(to_swap=hl.cond((table_result.gnomad_alleles[0]==table_result.data_allele[1]) & (table_result.gnomad_alleles[1]==table_result.data_allele[0]), True, False))
    #table_result = table_result.annotate(to_swap=hl.cond((str(table_result.gnomad_alleles[0])==str(table_result.data_allele[1]) and str(table_result.gnomad_alleles[1])==str(table_result.data_allele[0])), True, False))
    raw_ht = table_result.annotate(final_AF=hl.cond(table_result.to_swap == True, 1-table_result.gnomad_AF, table_result.gnomad_AF))

    # match keep only snps with either: (1) matching allele (our data and gnomad); or (2) 'swapped' alleles (to_swap is True)
    filtered_ht = raw_ht.annotate(keep=hl.cond((raw_ht.data_allele==raw_ht.gnomad_alleles) | (raw_ht.to_swap == True), True, False))
    #filtered_ht.show()
    filtered = filtered_ht.filter(filtered_ht.keep == True)
    filtered.export(dirname + '{}/'.format(cohortname) + cohortname + '.AlleleFreqCheck.tsv')

    # PART 2 (Mahalanobis Distance)
    file = pd.read_csv(dirname + '{}/'.format(cohortname) + cohortname + '.AlleleFreqCheck.tsv', sep='\t')
    file = file.dropna() # remove any NAs in the df

    df_x = file[['data_AF', 'final_AF']] # select only the allele frequency columns
    df_x['mahala'] = mahalanobis(x=df_x, data=file[['data_AF', 'final_AF']]) # comput MD and create new col
    md30 = df_x[df_x['mahala'] > 30] # select variants with MD > 30

    # plot
    plt.scatter(file['data_AF'], file['final_AF'], color = 'black', s = 0.3)
    plt.scatter(md30['data_AF'], md30['final_AF'], color = 'red', s = 0.3)
    plt.xlabel('data_AF')
    plt.ylabel('gnomad_AF')
    plt.savefig('/tmp/{}.alleleFreqCheck.png'.format(cohortname), dpi=300)
    hl.hadoop_copy('file:///tmp/{}.alleleFreqCheck.png'.format(cohortname), dirname + '{}/'.format(cohortname) + cohortname + '.alleleFreqCheck.png')
    #plt.savefig(dirname + 'alleleFreqCheck.png', dpi=300)

    rsids = pd.merge(md30, file, how='inner', on=['data_AF','final_AF']) # merge the outliers with the original dataset wo we can get rsids for outliers (MD > 30)
    rsids_list = list(set(rsids['data_rsid'])) # convert to list and remove any duplicates using set()

    mt = mt_input.filter_rows(hl.literal(rsids_list).contains(mt_input['rsid']), keep=False) # filter out the outliers
    print("The are {} SNPs in the final data after filtering allele checks".format(mt.count_rows()))

    return mt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str)

    parser.add_argument('--mt', type=str, help="Path to the matrix table you want to use for HWE")
    parser.add_argument('--plink', type=str, help="Path to plink files (prefix only)")
    parser.add_argument('--pcasamples', type=str,
                        help="File from PCA containing list of samples to be used in HWE filtering")
    parser.add_argument('--phenoType', type=str)
    parser.add_argument('--cohortsFile', type=str, help="Path to file containing list of cohorts you want to run HWE for")
    args = parser.parse_args()

    mt = hl.import_plink(bim=args.plink + ".bim",
                         fam=args.plink + ".fam",
                         bed=args.plink + ".bed",
                         min_partitions=100)

    mt = mahalanobis_distance(mt, args.dirname, "")


if __name__ == '__main__':
    main()


# EXAMPLE
# Running Allele Frequency and HWE checks
# hailctl dataproc submit hail 2allele_hwe_checks.py
# --dirname gs://lindo/test/
# --basename plus_cov_illumina_14082020.chr0.pos0.removed
# --mt gs://lindo/test/plus_cov_illumina_14082020.chr0.pos0.removed.qc.mt
# --pcasamples gs://lindo/test/eur.samples.to.keep
# --phenoType Control
# --cohortsFile gs://lindo/test/cohorts.txt