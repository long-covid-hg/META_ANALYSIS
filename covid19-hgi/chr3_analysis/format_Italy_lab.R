setwd("/home/tomoko/01.data.QC")
library(xlsx)
library(dplyr)
library(tidyr)
a <- read.xlsx("../covid19-hgi-clinical-values/hgi_italy/Clinical_characteristics_1001_COVID19_patients.xlsx", sheetName = "Foglio1")#1002
a <- a[!is.na(a$Sample.ID),]#1001
b <- read.xlsx("../covid19-hgi-clinical-values/hgi_italy/Clinical_characteristics_37_COVID19_patients_to substitute.xlsx", sheetName = "Foglio1")#1002

for(i in b$Sample.ID){
  a[a$Sample.ID == i,] <- b[b$Sample.ID == i,]
}
a <- a %>% rename(anonymized_patient_id = Sample.ID, age_at_diagnosis = Age, sex = Gender, ancestry = Ethnicity)
a <- a %>% mutate(lab_result_date = 0, lab_fibrinogen = as.numeric(Fibrinogen..mg.dl..n.v..Klauss.200.400..min.value.)/100,
                  lab_d_dimer = as.numeric(D.dimer..ng.ml..n.v...500..max.value.)/1000) 
#fibrinogen value something wrong!!!
a <- a %>% rename(lab_cd4 = CD4.lymphocytes..cell.ul..mm.3...n.v..400.1500..worst.value., lab_nk = NK.cells..cell.ul..mm..3...n.v..90.590..worst.value.,
                  lab_n_l_ratio = Neutrophils..min.value..Lymphocytes..min.value.., lab_il_6 = IL6..pg.ml....n.v...basal.value.,
                  lab_crp = CRP..mg.dl..n.v...0.5..max.value., lab_ldh = LDH..UI.l..n.v..M.135.225..F.135.214..max.value.,
                  )

blood_values <- c("lab_result_date",
                  "lab_wbc",
                  "lab_lymphocytes",
                  "lab_cd4",
                  "lab_cd8",
                  "lab_neutrophils",
                  "lab_monocytes",
                  "lab_platelets",
                  "lab_eosinophils",
                  "lab_basophils",
                  "lab_crp",
                  "lab_trop_t",
                  "lab_ast",
                  "lab_alt",
                  "lab_bilirubin",
                  "lab_ldh",
                  "lab_ggt",
                  "lab_alp",
                  "lab_d_dimer",
                  "lab_il_6",
                  "lab_serum_ferritin",
                  "lab_procalcitonin",
                  "lab_ck",
                  "lab_fibrinogen",
                  "lab_appt",
                  "lab_inr",
                  "lab_creatinine", "lab_nk", "lab_n_l_ratio")

identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex","ancestry","height","weight","smoking")
italy_lab <- a %>% select(identification, any_of(blood_values))
italy_lab[,-1] <- apply(italy_lab[,-1], 2, function(x){as.numeric(x)})

write.table(italy_lab, file="Italy.lab_20201008TN.tsv", sep="\t", col.names = T, row.names = F, quote = F)
