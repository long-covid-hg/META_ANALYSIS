setwd("/home/tomoko/01.data.QC/hgi_spain_bujanda/")
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)

path <- "../../covid19-hgi-clinical-values/"

## READ IN SPANISH BUJANDA DATA
spain_bujanda_manifest <- read.xlsx(paste0(path,"hgi_spain_bujanda/international_database_Spanish_Bujanda_v5_editedTN.xlsx"),sheetName = "Spanish patients_One-visit",stringsAsFactors=FALSE,colClasses="character")
spain_bujanda_lab <- read.xlsx(paste0(path,"hgi_spain_bujanda/international_database_Spanish_Bujanda_v5_editedTN.xlsx"),sheetName = "Laboratory data_time",stringsAsFactors=FALSE)

## Variable name dictionary
identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex")

hospital_info <- 
  c("highest_respiratory_support",
    "highest_who_score")

covid19 <-
  c("covid19_test_date")

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
                  "lab_trop_i",
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
                  "lab_aptt",
                  "lab_inr",
                  "lab_creatinine",
                  "lab_nk", "lab_n_l_ratio","lab_na","lab_k", "lab_mcv")

### HARMONIZING SPANISH BUJANDA DATA
#remove non-lab-confirmed cases
spain_bujanda_manifest <- spain_bujanda_manifest[spain_bujanda_manifest$rt.PCR.date != "1899-12-29" & !is.na(spain_bujanda_manifest$anonymized_patient_id),]

spain_bujanda_lab <- spain_bujanda_lab %>% filter(anonymized_patient_id %in% spain_bujanda_manifest$anonymized_patient_id)

spain_bujanda_lab <- spain_bujanda_lab %>% mutate(lab_wbc = lab_wbc..K.uL., lab_lymphocytes = lab_lympho..K.uL.,
                                                  lab_neutrophils = lab_neutrophils..K.uL., lab_platelets = lab_platelets..K.uL.,
                                                  lab_crp = lab_crp..mg.dL., lab_trop_t = lab_trop_t..pg.mL., lab_ast = lab_ast..IU.L.,
                                                  lab_alt = lab_alt..IU.L., lab_bilirubin = lab_bilirubin..mg.dl., lab_ldh = lab_ldh..U.L.,
                                                  lab_ggt = lab_ggt..U.L., lab_alp = lab_alp..U.L., lab_il_6 = lab_il_6..pg.mL.,
                                                  lab_serum_ferritin = lab_serum_ferritin..ng.mL., lab_procalcitonin = lab_procalcitonin..ng.mL.,
                                                  lab_ck = lab_ck..IU.L., lab_fibrinogen = lab_fibrinogen..g.L., lab_creatinine = lab_creatinine..mg.dL.,
                                                  lab_d_dimer = lab_d_dimer..mg.L.)

spain_bujanda_lab_rev <- merge(spain_bujanda_lab, spain_bujanda_manifest, by="anonymized_patient_id", all=T)

spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - as.Date(rt.PCR.date)))

spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% mutate(study = "spain_bujanda")
spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% rename(lab_aptt = lab_appt, lab_inr = lab_inr.)
spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% rename(covid19_test_date = rt.PCR.date) %>% mutate_at(vars(date_of_death), function(x){format(as.Date(x, origin = "1899-12-30"))})
spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% mutate(lab_fibrinogen = lab_fibrinogen/100)

c <- spain_bujanda_lab_rev %>% select(any_of(c(identification, demographics,blood_values,hospital_info, "study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(vars(colnames(spain_bujanda_lab_rev)[colnames(spain_bujanda_lab_rev) %in% blood_values[-1]], "highest_respiratory_support"),as.numeric)
c <- c %>% mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_who_score", "highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))

#the number of measurements per biomarker
df <- c %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(anonymized_patient_id) 
ggplot(df,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#the mean of the number of measurements per biomarker
df <- c %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))
ggplot(df,aes(y=mean,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df <- c %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>% 
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,any_of(blood_values))) %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(any_of(blood_values),-lab_result_date),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(any_of(blood_values)) %>% mutate(value = ifelse(!is.na(value) & !is.infinite(value), value, NA))

df1 <- df %>% filter(name == blood_values[2])
ggplot(df1,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw()  +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
for(i in seq(3,30)){
  df1 <- df %>% filter(name == blood_values[i])
  ggplot(df1,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}

saveRDS(c, file="spain_bujanda_lab.rds")
