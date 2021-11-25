setwd("/home/tomoko/01.data.QC/hgi_germany_schulte/")
library(readxlsb)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)

path <- "../../covid19-hgi-clinical-values/"


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

## READ IN GERMANY SCHULTE
germany_schulte_manifest <- read_xlsb(paste0(path,"hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 2)
germany_schulte_lab <- read_xlsb(paste0(path,"hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 3)

scale <- read_xlsb(paste0(path,"../covid19-hgi-clinical-values/hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 4)
germany_schulte_labrev <- germany_schulte_lab %>% mutate_at(vars(-c(anonymized_patient_id ,lab_result_date)),function(x){as.numeric(x)}) %>% 
  mutate(lab_lymphocytes = lab_lymphocytes*lab_wbc, lab_neutrophils = lab_neutrophils*lab_wbc,
         lab_monocytes = lab_monocytes*lab_wbc, lab_eosinophils = lab_eosinophils**lab_wbc, 
         lab_basophils = lab_basophils*lab_wbc)
hosp2 <- c("TUM_058", "TUM_059", "TUM_060", "TUM_061")
germany_schulte_labrev <- germany_schulte_labrev %>% 
  mutate(lab_d_dimer = case_when(anonymized_patient_id %in% hosp2 ~ lab_d_dime,
                                 TRUE ~ lab_d_dime/1000),
         lab_crp = case_when(anonymized_patient_id %in% hosp2 ~ lab_crp/10,
                             TRUE ~ lab_crp),
         lab_trop_t = case_when(anonymized_patient_id %in% hosp2 ~ -1,
                                TRUE ~ lab_trop_t*1000),
         lab_trop_i = case_when(anonymized_patient_id %in% hosp2 ~ lab_trop_t,
                                TRUE ~ -1)
  )
germany_schulte_labrev  <- germany_schulte_labrev  %>%
  mutate(lab_fibrinogen = lab_fibrinogen/100)

tmp <- merge(germany_schulte_labrev, germany_schulte_manifest, by="anonymized_patient_id", all=T)
h <- tmp %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date, origin="1900-01-01") - as.Date(covid19_test_date))) %>% 
  mutate(study = "germany_schulte") %>%
  select(any_of(c(identification, demographics,hospital_info,blood_values,"study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support","highest_who_score")), .funs = funs(ifelse(. == -1,NA,.)))

#error correction
h$lab_wbc[h$anonymized_patient_id == "TUM_049"][1] <- 7.52
h$lab_na[h$anonymized_patient_id == "TUM_020" & h$lab_na > 200] <- 149
h$lab_mcv[h$anonymized_patient_id == "TUM_018" & h$lab_mcv > 300] <- 98

## START PROCESSING HOSPITAL INFORMATION

#the number of measurements per biomarker
df <- h %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(anonymized_patient_id) 
ggplot(df,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#the mean of the number of measurements per biomarker
df <- h %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))
ggplot(df,aes(y=mean,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df <- h %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))

df <- h %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>% 
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

saveRDS(h, file="germany_schulte_lab.rds")
