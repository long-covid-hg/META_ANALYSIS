setwd("/home/tomoko/01.data.QC/hgi_canada/")
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)

## READ IN CANADA DATA
canada_manifest <- fread("../../covid19-hgi-clinical-values/hgi_canada/bqc19_one_time_V2.tsv")
canada_lab <- fread("../../covid19-hgi-clinical-values/hgi_canada/bqc19_visit_time_V2.tsv")
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

tmp <- merge(canada_lab, canada_manifest, by="anonymized_patient_id", all.x=T, all.y=T)
tmp <- tmp %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)))
g <- tmp %>% mutate(study="canada") %>% select(any_of(c(identification, demographics,hospital_info,blood_values, "study"))) %>% mutate_at("anonymized_patient_id",as.character) %>% 
  mutate(sex = ifelse(sex=="F", as.character(1), as.character(0)))

g <- g %>% mutate(lab_d_dimer = lab_d_dimer/1000) %>% mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_who_score", "highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))

#the number of measurements per biomarker
df <- g %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(anonymized_patient_id) 
ggplot(df,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#the mean of the number of measurements per biomarker
df <- g %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))
ggplot(df,aes(y=mean,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df <- g %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>% 
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
### monocyte check ID 760 Na, K, mcv too!!
saveRDS(g, file="canada_lab.rds")
