import hail as hl
import argparse


def recode_GP(mt):
    """
    Recode GP entry field for males to be diploid
    :param mt: MatrixTable to recode
    :return: MatrixTable with GP field recoded as [GP(0|0), GP(0|1), 0] for Males samples (where GP has 2 elements)
    """
    mt = mt.annotate_entries(GP_new=hl.cond(hl.len(mt.GP) == 2,
                                            [mt.GP[0], mt.GP[1], 0],
                                            mt.GP))
    mt = mt.annotate_entries(GP=mt.GP_new)
    mt = mt.drop(mt.GP_new)
    return mt


def recode_GP2(mt):
    """
    Recode GP entry field for males
    :param mt: MatrixTable to recode
    :return: MatrixTable with GP field recoded as [GP(0|0), 0, GP(0|1)] for Males samples (where GP has 2 elements)
    """
    mt = mt.annotate_entries(GP_new=hl.cond(hl.len(mt.GP) == 2,
                                            [mt.GP[0], 0, mt.GP[1]],
                                            mt.GP))
    mt = mt.annotate_entries(GP=mt.GP_new)
    mt = mt.drop(mt.GP_new)
    return mt


def main(args):
    chrX = hl.import_vcf(args.dirname + args.vcfname,
                         force_bgz=True,
                         reference_genome='GRCh38',
                         entry_float_type=hl.dtype('float32'),
                         min_partitions=100)

    chrX = recode_GP(chrX).drop('HDS')

    chrX.describe()

    basename = args.vcfname.replace(".vcf.gz", "")
    hl.export_vcf(chrX,
                  args.dirname + basename + '.rec.vcf.bgz',
                  metadata=hl.get_vcf_metadata(args.dirname + args.vcfname))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type = str, required=True)
    parser.add_argument('--vcfname', type = str, required=True)

    args = parser.parse_args()

    print(args)
    main(args)


# Run from command line
# gcloud dataproc jobs submit pyspark recode_chrX.py \
# --cluster hail \
# --region europe-west1 \
# -- --dirname gs://dsge-covid19-data/ita_14082020/data_imputed/ --vcfname chrX.dose.vcf.gz