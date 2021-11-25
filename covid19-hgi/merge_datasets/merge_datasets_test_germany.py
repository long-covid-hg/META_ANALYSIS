import hail as hl
hl.init(default_reference="GRCh38")

# # # # # # # # # LOAD CONTROLS # # # # # # # # #
print("*******GERMANY CON*********")
ctl = hl.read_matrix_table("gs://dsge-covid19-data/controls_DE/df7/qc_step_4/germany/germany.hwe_filtered_qced.mt")
ctl = ctl.select_cols(ctl.is_female)
print(ctl.count())
print("")

print("*******GERMANY CAS*********")
cas = hl.read_matrix_table("gs://dsge-covid19-data/20211030/qc_step_4/germany/germany.hwe_filtered_qced.mt")
cas = cas.select_cols(cas.is_female)
print(cas.count())
print("")


print("*******MERGE GERMANY CAS + GERMANY CON*********")
mt = cas.union_cols(ctl)
print(mt.count())
print("")

mt.write("gs://dsge-covid19-data/controls_DE/20211030/qc_step_4/germany/germany_cases_controls.mt", overwrite=True)
hl.export_plink(mt, "gs://dsge-covid19-data/20211030/qc_step_4/germany/germany_cases_controls")
