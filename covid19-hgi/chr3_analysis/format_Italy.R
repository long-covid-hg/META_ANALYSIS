setwd("/home/tomoko/01.data.QC")
library(xlsx)
library(dplyr)
library(tidyr)
a <- read.xlsx("../covid19-hgi-clinical-values/hgi_italy/Clinical_characteristics_1001_COVID19_patients_rev.xlsx", sheetName = "Foglio1")#1002
a <- a[!is.na(a$Sample.ID),]#1001
b <- read.xlsx("../covid19-hgi-clinical-values/hgi_italy/Clinical_characteristics_37_COVID19_patients_to substitute.xlsx", sheetName = "Foglio1")#1002

for(i in b$Sample.ID){
  a[a$Sample.ID == i,] <- b[b$Sample.ID == i,]
}

a <- a %>% rename(anonymized_patient_id = Sample.ID, age_at_diagnosis = Age, sex = Gender, ancestry = Ethnicity)
a1 <- a %>% mutate(height = -1, weight = -1, smoking = -1)

a2 <- a1 %>% mutate(hospitalization = ifelse(Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 0, 0, 1))
a2 <- a2 %>% mutate(hospitalization_start = -1, hospitalization_end = -1)

a3 <- a2 %>% mutate(death = ifelse(Status..A..Alive..D.Dead. == "D", 1, 0))
a3 <- a3 %>% mutate(highest_respiratory_support = case_when(Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 4 ~ 1,
                                                            Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 3 ~ 2,
                                                            TRUE ~ 0))

a3 <- a3 %>% mutate(cause_of_death = -1, date_of_death = -1, icu_admit = -1, icu_duration = -1)

a3 <- a3 %>% mutate(days_ventilator = -1, hosp_dvt = -1, hosp_thrombo = -1, hosp_stroke = -1, hosp_bleeding = -1)

#According to their definition, they used T-Troponin >15 ng/L as a definition of heart involvement. But this is not always true..
a3 <- a3 %>% mutate(hosp_infarction = ifelse(grepl("T", Heart.involvement..T.T.Troponin..15.ng.L..B..NT.proBNP.M..88..F..153.pg.ml..A.arrhythmia.), 1, 0))
#According to their definition, Kidney involvement (creatinine M > 0.7-1.20; F > 0.5-1.10 mg/dl)
a3 <- a3 %>% mutate(hosp_renal = ifelse(Kidney.involvement..creatinine.M...0.7.1.20..F...0.5.1.10.mg.dl. == "Yes", 1, 0))
#According to their definition, Liver involvement (both ALT and AST > 40 U/L)
a3 <- a3 %>% mutate(hosp_hepatic = ifelse(Liver.involvement..both.ALT.and.AST...40.U.L. == "Yes", 1, 0))
                      
##COVID test
a3 <- a3 %>% mutate(covid19_test = 1, covid19_test_date = -1, covid19_test_type = -1, covid19_first_symptoms_date = -1)

##medication
a4 <- a3 %>% mutate(steroids = ifelse(Steroids..M.Methylprednisolone..D.Dexamethasone. == "No", 0, 1))
a4 <- a4 %>% mutate(biologics = ifelse(Immunotherapy..T.Tocilizumab..R.Ruxolitinib..B.Baricitinib. == "No", 0, 1))
a4 <- a4 %>% mutate(lmwh = ifelse(Enoxaparin..Yes.No. == "No", 0, 1))
a4 <- a4 %>% mutate(hcq = ifelse(grepl("HCQ", Other.anti.COVID19.drugs..HCQ.Hydroxicloroquine..A.Azithromycine..REM.Remdesevir..L.R.Lopinavir.Ritonavir..C.Cloroquine..D.R..Darunavir.Ritonavir..), 1, 0))
a4 <- a4 %>% mutate(remdesivir = ifelse(grepl("REM", Other.anti.COVID19.drugs..HCQ.Hydroxicloroquine..A.Azithromycine..REM.Remdesevir..L.R.Lopinavir.Ritonavir..C.Cloroquine..D.R..Darunavir.Ritonavir..), 1, 0))

##highest who score
a5 <- a4 %>% mutate(highest_who_score = case_when(Status..A..Alive..D.Dead. == "D" ~ 10,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 4 & Respiratory.Severity.P.F.score..worst.value. < 150 ~ 8,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 4 & !(Respiratory.Severity.P.F.score..worst.value. < 150) ~ 7,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 3 ~ 6,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 2 ~ 5,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 1 ~ 4,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 0 & grepl("Yes", Taste.Smell.system.involvement..B..Both..T..A.Hypogeusia..S..A.Hyposmia) ~ 2,
                                                  Clinical.Category..4.intubation..3.CPAP.BiPAP..2.Oxygen.therapy..1.no.Oxygen.therapy..0.not.hospitalized. == 0 & !(grepl("Yes", Taste.Smell.system.involvement..B..Both..T..A.Hypogeusia..S..A.Hyposmia)) ~ 1,
                                                  TRUE ~ 0))

identification <- c("anonymized_patient_id")
demographics <- c("age_at_diagnosis",
                  "sex",
                  "ancestry",
                  "height",
                  "weight",
                  "smoking")
hosps <- c("hospitalization",
                     "hospitalization_start",
                     "hospitalization_end",
                     "death",
                     "cause_of_death",
                     "date_of_death",
                     "icu_admit",
                     "icu_duration",
                     "highest_who_score",
                     "highest_respiratory_support",
                     "days_ventilator",
                     "hosp_dvt",
                     "hosp_thrombo",
                     "hosp_stroke",
                     "hosp_infarction",
                     "hosp_renal",
                     "hosp_bleeding",
                     "hosp_hepatic")

covid <- c("covid19_test",
           "covid19_test_date",
           "covid19_test_type",
           "covid19_first_symptoms_date")

medications <- c("steroids",
                 "biologics",
                 "lmwh",
                 "hcq",
                 "remdesivir")

#without_com <- a5 %>% select(identification, demographics, hosps, covid, medications)
#without_com$sex <- ifelse(without_com$sex == "F", 1, 0)
#without_com$ancestry <- case_when(grepl("White", without_com$ancestry) ~ 0,
                                  #grepl("Black", without_com$ancestry) ~ 1,
                                  #grepl("Hispanic", without_com$ancestry) ~ 2,
                                  #TRUE ~ -1)

#write.table(without_com, "Italy.withoutCom.tsv", sep="\t", col.names = T, row.names = F, quote = F)

#############to be continued.

a5 <- a5 %>% mutate(sex = ifelse(sex == "F", 1, 0), ancestry = case_when(grepl("White", ancestry) ~ 0,
                                                                         grepl("Black", ancestry) ~ 1,
                                                                         grepl("Hispanic", ancestry) ~ 2,
                                                                         TRUE ~ -1))
##comorbidities
COM <- a$Clinical.known.co.morbidities
COM1 <- unique(unlist(strsplit(COM, ", ")))
a5 <- a5 %>% mutate(com_hiv = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                        grepl("HIV",Clinical.known.co.morbidities) ~ 1, TRUE ~ 0))
a5 <- a5 %>% mutate(com_immunocomp = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                        grepl("Immunodeficiency",Clinical.known.co.morbidities) ~ 1, TRUE ~ 0))
a5 <- a5 %>% mutate(com_transplant = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                               grepl("transplant",Clinical.known.co.morbidities, ignore.case = T) ~ 1, TRUE ~ 0))
a5 <- a5 %>% mutate(com_autoimm_rheum = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                  grepl("rheu",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Inflammatory disease",Clinical.known.co.morbidities, ignore.case = T) ~ 1, 
                                                  grepl("RCU",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("raynaud's disease",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("psoriatic arthritis",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Autoimmune",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Sjogren",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("myositys",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Wegener",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Spondiloentesoartrite",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("thyroiditis",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("spondy",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  grepl("Crohn",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  TRUE ~ 0))
a5 <- a5 %>% mutate(com_type_ii_diabetes = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             grepl("DM",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             grepl("Diabetes",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             grepl("diabetic",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             TRUE ~ 0))
a5 <- a5 %>% mutate(com_type_i_diabetes = -1)
a5 <- a5 %>% mutate(com_asthma = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                           grepl("asthma",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           grepl("Astma",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             TRUE ~ 0))
a5 <- a5 %>% mutate(com_chronic_pulm = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                 grepl("BPCO",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("COLD",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Chronic bronchitis",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("COPD",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Emphysema",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Emphisema",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           TRUE ~ 0))
a5 <- a5 %>% mutate(com_sleep_apnea = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                grepl("OSAS",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 TRUE ~ 0))
a5 <- a5 %>% mutate(com_liver = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                          grepl("Hypertransaminasemia",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("Liver",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("HCV",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("hepatitis",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("hepatic",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("HBV",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("transaminases",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("Gilbert",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                          grepl("Hepatitis B",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           TRUE ~ 0))
a5 <- a5 %>% mutate(com_gallbl = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                           grepl("Gallbladder stones",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           grepl("cholecystectomy for calculus",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                TRUE ~ 0))
a5 <- a5 %>% mutate(com_pancreas = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             grepl("pancreas",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             TRUE ~ 0))
a5 <- a5 %>% mutate(com_chronic_kidney = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                   grepl("Kidney",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                   grepl("Polycystic Kidney",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                   grepl("glomerulonephritis",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                   grepl("CKD",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 TRUE ~ 0))
a5 <- a5 %>% mutate(com_dialysis = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             com_chronic_kidney == 1 ~ -1,
                                             TRUE ~ 0))
a5 <- a5 %>% mutate(com_heart_failure = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                  grepl("Congestive Heart Failure",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 TRUE ~ 0))
a5 <- a5 %>% mutate(com_hypertension = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                 grepl("pulmonary hypertension",Clinical.known.co.morbidities, ignore.case = T) ~ 0,
                                                 grepl("Hypertension",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("high blood pressure",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("hypertensive",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 TRUE ~ 0))
a5 <- a5 %>% mutate(com_infarction = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                               grepl("Ischemic Cardiopathy",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("AMI",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("Ischemic heart disease",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("myocardial infar heart attack",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("Acute myocardial infarction",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("coronary artery disease",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                               grepl("IMA",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                  TRUE ~ 0))
a5 <- a5 %>% mutate(com_vascular = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             grepl("PTA",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             grepl("Periferical Vascular diseae",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             TRUE ~ 0))

a5 <- a5 %>% mutate(com_stroke = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                           grepl("Ictus",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           grepl("subarachnoid hemorrhage",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           grepl("cerebral",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           grepl("stroke",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                    TRUE ~ 0))
a5 <- a5 %>% mutate(com_afib = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                         grepl("FA",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                         grepl("AF",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                         grepl("Fibrillation",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                           TRUE ~ 0))
a5 <- a5 %>% mutate(com_dementia = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                 grepl("Cognitive",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Dementia",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Alzheimer",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                TRUE ~ 0))



a5 <- a5 %>% mutate(com_neurological = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                 grepl("Parkinson",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Schizofrenia",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Neurological disease",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Mood disorder",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Psichiatric",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("migraine",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Bipolar",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("anxiety",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("epilepsy",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Emicrania",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("Epilessia",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 grepl("SAD",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                         TRUE ~ 0))

a5 <- a5 %>% mutate(com_leukemia = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             grepl("Leukemia",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                    TRUE ~ 0))
a5 <- a5 %>% mutate(com_lymphoma = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                             grepl("lynphoma",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             grepl("lymphoma",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             grepl("DLB-CL",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                             TRUE ~ 0))

a5 <- a5 %>% mutate(com_malignant_solid = case_when(grepl("unknown",Clinical.known.co.morbidities, ignore.case = T) ~ -1,
                                                    grepl("cancer",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                    grepl("tumor",Clinical.known.co.morbidities, ignore.case = T) ~ 1,
                                                 TRUE ~ 0))


com <- c("com_hiv",
         "com_immunocomp",
         "com_transplant",
         "com_autoimm_rheum",
         "com_type_ii_diabetes",
         "com_type_i_diabetes",
         "com_asthma",
         "com_chronic_pulm",
         "com_sleep_apnea",
         "com_liver",
         "com_gallbl",
         "com_pancreas",
         "com_chronic_kidney",
         "com_dialysis",
         "com_heart_failure",
         "com_hypertension",
         "com_infarction",
         "com_vascular",
         "com_stroke",
         "com_afib",
         "com_dementia",
         "com_neurological",
         "com_leukemia",
         "com_lymphoma",
         "com_malignant_solid")

italy <- a5 %>% select(identification, demographics, hosps, covid, medications, com)
head(italy)

write.table(italy, "Italy.withCom_20201008TN.tsv", sep="\t", col.names = T, row.names = F, quote = F)
