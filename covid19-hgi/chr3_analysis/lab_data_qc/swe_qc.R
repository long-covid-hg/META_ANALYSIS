setwd("/home/tomoko/01.data.QC/hgi_sweden/")
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)

path <- "../../covid19-hgi-clinical-values/"

## READ IN SWEDISH DATA
swe_manifest <- read.xlsx(paste0(path,"hgi_swe/PronMed genetic fenotypes.xlsx"),sheetName = "Data",stringsAsFactors=FALSE)
swe_lab <- read.xlsx(paste0(path,"hgi_swe/PronMed blood chemistry v2_editedAG.xlsx"),sheetName = "values",stringsAsFactors=FALSE)
swe_lab_plus <- read.xlsx(paste0(path,"hgi_swe/PronMedBloodChemistry_v2_editedTN.xlsx"),sheetIndex = "New Labs",stringsAsFactors=FALSE)

swe_lab_rev <- merge(swe_lab, swe_lab_plus, by=c("anonymized_patient_id", "lab_result_date.relative.ICU.admission"), all = T)

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

### HARMONIZING SWEDISH DATA
swe_lab <- swe_lab_rev %>% rename(lab_wbc = lab_leukocytes) %>% mutate(study = "sweden")

#impute by .manual if the value is missing
swe_lab$lab_lymphocytes[is.na(swe_lab$lab_lymphocytes)] <- swe_lab$lab_lymphocytes..manual.[is.na(swe_lab$lab_lymphocytes)]
swe_lab$lab_neutrophils[is.na(swe_lab$lab_neutrophils)] <- swe_lab$lab_neutrophils.manual[is.na(swe_lab$lab_neutrophils)]
swe_lab$lab_monocytes[is.na(swe_lab$lab_monocytes)] <- swe_lab$lab_monocytes.manual[is.na(swe_lab$lab_monocytes)]
swe_lab$lab_eosinophils[is.na(swe_lab$lab_eosinophils)] <- swe_lab$lab_eosinophils.manual[is.na(swe_lab$lab_eosinophils)]
swe_lab$lab_basophils[is.na(swe_lab$lab_basophils)] <- swe_lab$lab_basophils.manual[is.na(swe_lab$lab_basophils)]
#CRP mg/L > mg/dL
swe_lab$lab_crp <- swe_lab$lab_crp/10
#trop i to trop T??
swe_lab <- swe_lab %>% rename(lab_trop_i = lab_trop_t)

#conversion from ukat/L to U/L
swe_lab$lab_ast <- swe_lab$lab_ast/0.0166
swe_lab$lab_alt <- swe_lab$lab_alt/0.0166
swe_lab$lab_ggt <- swe_lab$lab_ggt/0.0166
#convert from umol/L to mg/dL
swe_lab$lab_bilirubin <- swe_lab$lab_bilirubin/17.1
#convert D-Dimer mg/L to FEU ug/L  ## double check with them!
swe_lab$lab_d_dimer <- swe_lab$lab_d_dime
#creatinine convert from umol/L to mg/dL
swe_lab$lab_creatinine <- swe_lab$lab_creatinine * 0.0113
##!!## at this moment use relative date to ICU admission but need COVID test date. 
swe_lab <- swe_lab %>% rename(lab_result_date = lab_result_date.relative.ICU.admission)

# add new lab value
swe_lab <- swe_lab %>% mutate(lab_alp = 59.99*P.ALP, lab_aptt = P.APTT, lab_fibrinogen = P.Fibrinogen, 
                              lab_ck = 59.99*P.CK, lab_ldh = 59.99*P.LD, lab_inr = P.Protrombinkomplex)

swe_lab_rev <- merge(swe_lab, swe_manifest, by="anonymized_patient_id", all=T)

swe_lab <- swe_lab_rev %>% select(any_of(c(identification, blood_values, demographics,hospital_info,covid19,"study")))

#write.table(swe_lab, "swe_lab_cleanedTN20200925.tsv", col.names = T, row.names = F, sep="\t", quote = F)

b <- swe_lab %>% select(any_of(c(identification, demographics,blood_values,hospital_info,"study"))) %>% mutate_at(vars("anonymized_patient_id","sex"),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))

#the number of measurements per biomarker
df <- b %>% mutate_at(.vars = vars(any_of(-anonymized_patient_id)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support)) %>%
  group_by(anonymized_patient_id) 
ggplot(df,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#the mean of the number of measurements per biomarker
df <- b %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support)) %>%
  group_by(name) %>% summarize(mean=mean(value))
ggplot(df,aes(y=mean,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df <- b %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>% 
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,any_of(blood_values))) %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(any_of(blood_values),-lab_result_date),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(any_of(blood_values)) %>% mutate(value = ifelse(!is.na(value) & !is.infinite(value), value, NA))

df1 <- df %>% filter(name == blood_values[2])
ggplot(df1,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
for(i in seq(3,30)){
  df1 <- df %>% filter(name == blood_values[i])
  ggplot(df1,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
}

saveRDS(b, file="swe_lab.rds")
