import hail as hl

hl.init(default_reference="GRCh38")

mt_19032021 = hl.import_plink(bed='gs://dsge-covid19-data/15042021/1903/PLUS_COV_ILLUMINA_19032021.chr0.pos0.samples.removed.bed',
                              bim='gs://dsge-covid19-data/15042021/1903/PLUS_COV_ILLUMINA_19032021.chr0.pos0.samples.removed.bim',
                              fam='gs://dsge-covid19-data/15042021/1903/PLUS_COV_ILLUMINA_19032021.chr0.pos0.samples.removed.fam')
print("mt_19032021")
print(mt_19032021.count_rows(), mt_19032021.count_cols())

mt_15122020 = hl.import_plink(bed='gs://dsge-covid19-data/cov_illumina_15122020/PLUS_COV_ILLUMINA_15122020.chr0.pos0.samples.removed.bed',
                              bim='gs://dsge-covid19-data/cov_illumina_15122020/PLUS_COV_ILLUMINA_15122020.chr0.pos0.samples.removed.bim',
                              fam='gs://dsge-covid19-data/cov_illumina_15122020/PLUS_COV_ILLUMINA_15122020.chr0.pos0.samples.removed.fam')
print("mt_15122020")
print(mt_15122020.count_rows(), mt_15122020.count_cols())


mt_sweden = hl.import_plink(bed='gs://dsge-covid19-data/cov_uppsala_15042021/TC-2500_210205.GRCh38.chr0.pos0.removed.ID.updated.bed',
                            bim='gs://dsge-covid19-data/cov_uppsala_15042021/TC-2500_210205.GRCh38.chr0.pos0.removed.ID.updated.bim',
                            fam='gs://dsge-covid19-data/cov_uppsala_15042021/TC-2500_210205.GRCh38.chr0.pos0.removed.ID.updated.fam')
print("Sweden")
print(mt_sweden.count_rows(), mt_sweden.count_cols())

merged = mt_15122020.union_cols(mt_19032021).union_cols(mt_sweden)
print("Merged")
print(merged.count_rows(), merged.count_cols())

merged = merged.repartition(n_partitions=100, shuffle=True)
merged.write("gs://dsge-covid19-data/15042021/cov_15042021.merged.mt", overwrite=True)

hl.export_plink(merged, "gs://dsge-covid19-data/15042021/cov_15042021.merged")
