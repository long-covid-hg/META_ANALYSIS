import hail as hl
import argparse

hl.init(default_reference='GRCh38')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--outname', type=str)

    parser.add_argument('--mt', type=str, help="Path to the matrix table you want to use")
    parser.add_argument('--pcasamples', type=str,
                        help="File from PCA containing list of samples to be filtered")

    args = parser.parse_args()

    mt = hl.read_matrix_table(args.mt)
    keep = hl.import_table(args.pcasamples, no_header=True)
    # keep = hl.import_table(samplestokeep, key='01947TG2')
    keep = keep.key_by('f0')
    # mt = mt.filter_cols(hl.is_defined(keep.index(mt['s'])))
    mt = mt.semi_join_cols(keep)
    print("Number of samples in file: {}".format(mt.count_cols()))

    hl.export_plink(mt,
                    output=args.dirname + args.outname,
                    is_female=mt.is_female)


if __name__ == '__main__':
    main()
