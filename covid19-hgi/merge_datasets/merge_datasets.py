import hail as hl
hl.init(default_reference="GRCh38")

# # # # # # # # # LOAD CONTROLS # # # # # # # # #
print("*******GERMANY CONTROLS*********")
ctl_de = hl.read_matrix_table("gs://dsge-covid19-data/controls_DE/PsyCourse_Controls/germany.hwe_filtered_qced.chrN-updated.mt")
ctl_de = ctl_de.select_cols(ctl_de.is_female)
print(ctl_de.count())
print("")


print("*******SWEDEN CONTROLS*********")
ctl_swe = hl.read_matrix_table("gs://dsge-covid19-data/swe_controls_ANGISE/ANGISE_controls_EUR_ensembl-updated.Sweden.hwe_filtered_qced.mt")
ctl_swe = ctl_swe.select_cols(ctl_swe.is_female)
print(ctl_swe.count())
print("")

print("*******BELGIUM CONTROLS*********")
ctl_broad_be = hl.read_matrix_table("gs://dsge-covid19-data/controls_broad/controls_merged_EUR_ensembl-updated.belgium.hwe_filtered_qced.belgium.hwe_filtered_removed_palindromic_qced.mt")
ctl_broad_be = ctl_broad_be.select_cols(ctl_broad_be.is_female)
print(ctl_broad_be.count())
print("")

print("*******ITALY CONTROLS*********")
ctl_broad_ita = hl.read_matrix_table("gs://dsge-covid19-data/controls_broad/controls_merged_EUR_ensembl-updated.ita.hwe_filtered_qced.mt")
ctl_broad_ita = ctl_broad_ita.select_cols(ctl_broad_ita.is_female)
print(ctl_broad_ita.count())
print("")

loci_to_remove = [hl.parse_locus('chr19:42779471', reference_genome='GRCh38')]
ctl_broad_ita = ctl_broad_ita.filter_rows(hl.literal(loci_to_remove).contains(ctl_broad_ita.locus), keep=False)

# # # # # # # # # LOAD CASES # # # # # # # # #
basedir = "gs://dsge-covid19-data/15042021/2/qc_step_5/"

print("*******BELGIUM CASES*********")
cases_be = hl.import_plink(bim=basedir + "belgium.hwe_filtered_qced_N-updated.bim",
                           bed=basedir + "belgium.hwe_filtered_qced_N-updated.bed",
                           fam=basedir + "belgium.hwe_filtered_qced_N-updated.fam")
cases_be = cases_be.select_cols(cases_be.is_female)
cases_be = cases_be.repartition(n_partitions=100, shuffle=True)
print(cases_be.count())
print("")

print("*******ITA CASES*********")
cases_ita = hl.import_plink(bim=basedir + "italy.hwe_filtered_qced_N-updated.bim",
                            bed=basedir + "italy.hwe_filtered_qced_N-updated.bed",
                            fam=basedir + "italy.hwe_filtered_qced_N-updated.fam")
cases_ita = cases_ita.select_cols(cases_ita.is_female)
cases_ita = cases_ita.repartition(n_partitions=100, shuffle=True)
print(cases_ita.count())
print("")

print("*******DE CASES*********")
cases_de = hl.import_plink(bim=basedir + "germany.hwe_filtered_qced_N-updated.bim",
                           bed=basedir + "germany.hwe_filtered_qced_N-updated.bed",
                           fam=basedir + "germany.hwe_filtered_qced_N-updated.fam")
cases_de = cases_de.select_cols(cases_de.is_female)
cases_de = cases_de.repartition(n_partitions=100, shuffle=True)
print(cases_de.count())
print("")

print("*******IRAN CASES*********")
cases_iran = hl.import_plink(bim=basedir + "iran.hwe_filtered_qced_N-updated.bim",
                             bed=basedir + "iran.hwe_filtered_qced_N-updated.bed",
                             fam=basedir + "iran.hwe_filtered_qced_N-updated.fam")
cases_iran = cases_iran.select_cols(cases_iran.is_female)
cases_iran = cases_iran.repartition(n_partitions=100, shuffle=True)
print(cases_iran.count())
print("")

print("*******EGYPT CASES*********")
cases_egypt = hl.import_plink(bim=basedir + "egypt.hwe_filtered_qced_N-updated.bim",
                              bed=basedir + "egypt.hwe_filtered_qced_N-updated.bed",
                              fam=basedir + "egypt.hwe_filtered_qced_N-updated.fam")
cases_egypt = cases_egypt.select_cols(cases_egypt.is_female)
cases_egypt = cases_egypt.repartition(n_partitions=100, shuffle=True)
print(cases_egypt.count())
print("")

print("*******SWE CASES*********")
cases_swe = hl.import_plink(bim=basedir + "sweden.hwe_filtered_qced_N-updated.bim",
                            bed=basedir + "sweden.hwe_filtered_qced_N-updated.bed",
                            fam=basedir + "sweden.hwe_filtered_qced_N-updated.fam")
cases_swe = cases_swe.select_cols(cases_swe.is_female)
cases_swe = cases_swe.repartition(n_partitions=100, shuffle=True)
print(cases_swe.count())
print("")

print("*******SWE CASES UPPSALA*********")
basedir = "gs://dsge-covid19-data/cov_uppsala_15042021/qc_step_5/"
cases_swe2 = hl.import_plink(bim=basedir + "sweden.hwe_filtered_qced_N-updated.bim",
                             bed=basedir + "sweden.hwe_filtered_qced_N-updated.bed",
                             fam=basedir + "sweden.hwe_filtered_qced_N-updated.fam")
cases_swe2 = cases_swe2.select_cols(cases_swe2.is_female)
cases_swe2 = cases_swe2.repartition(n_partitions=100, shuffle=True)
print(cases_swe2.count())
print("")

# # # # # # # # # MERGE DATASETS # # # # # # # # #
print("*******MERGING CASES FIMM+CASES UPPSALA*********")
cases = cases_be.union_cols(cases_ita).union_cols(cases_de).union_cols(cases_iran).union_cols(cases_egypt).union_cols(cases_swe).union_cols(cases_swe2)
print(cases.count())
print("")

cases.write("gs://dsge-covid19-data/15042021/data_chr3/cases.main.ancestry.mt", overwrite=True)

cases = hl.read_matrix_table("gs://dsge-covid19-data/15042021/data_chr3/cases.main.ancestry.mt")
hl.export_plink(cases, "gs://dsge-covid19-data/15042021/data_chr3/cases.main.ancestry")

print("*******MERGING CASES FIMM+CASES UPPSALA+CONTROLS WO/ GERMANY*********")
mt = ctl_broad_ita.union_cols(ctl_broad_be).union_cols(ctl_swe).union_cols(cases)
print(mt.count())
print("")

mt.write("gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.wo.germany.mt", overwrite=True)

mt = hl.read_matrix_table("gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.wo.germany.mt")

hl.export_plink(mt, "gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.wo.germany")
print(mt.count())
print("")

print("*******MERGING CASES FIMM+CASES UPPSALA+CONTROLS+CONTROLS GERMANY*********")
mt2 = mt.union_cols(ctl_de)
print(mt2.count())
print("")

mt2.write("gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.w.germany.mt", overwrite=True)

mt2 = hl.read_matrix_table("gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.w.germany.mt")
hl.export_plink(mt2, "gs://dsge-covid19-data/15042021/qc_step_5/cov_15042021.controls.w.germany")
