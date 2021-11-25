setwd("/home/tomoko/01.data.QC/hgi_belgium")
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)

path <- "../../covid19-hgi-clinical-values/"

## READ IN BELGIUM DATA
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-111020_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)
bel_manifest <- bel_manifest[-dim(bel_manifest)[1],]
bel_lab <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-111020_TN.xls"),sheetName = "visit_time",stringsAsFactors=FALSE)
bel_lab <- bel_lab[!is.na(bel_lab$anonymized_patient_id),]
bel_lab <- bel_lab[-dim(bel_lab)[1],]

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

## START PROCESSING HOSPITAL INFORMATION

### HARMONIZING BELGIUM DATA
a <- bel_manifest %>% select(any_of(c(identification,demographics,hospital_info,covid19))) %>% 
  mutate(study="belgium", covid19_test_date=as.Date(covid19_test_date,"%d/%m/%Y")) 

a <- a %>% mutate(highest_who_score = case_when(highest_who_score == "8/9" ~ 8,
                                                TRUE ~ as.numeric(highest_who_score)))
a <- a %>% mutate(highest_respiratory_support = as.numeric(highest_respiratory_support))

bel_lab_rev <- merge(bel_lab, a, by="anonymized_patient_id", all.x=T, all.y=T)
bel_lab_rev <- bel_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date,"%d/%m/%Y") - as.Date(covid19_test_date))) %>% 
  rename(lab_aptt = APTT, lab_inr = INR)

set_numeric <- function(x){as.numeric(x)}
bel_lab_rev <- bel_lab_rev %>% mutate_at(vars(c(any_of(blood_values),-lab_result_date)), set_numeric)
bel_lab_rev <- bel_lab_rev %>% mutate(lab_lymphocytes = lab_wbc * lab_lymphocytes/100, lab_neutrophils = lab_wbc * lab_neutrophils/100,
                                      lab_monocytes = lab_wbc * lab_monocytes/100, lab_eosinophils = lab_wbc * lab_eosinophils/100,
                                      lab_basophils = lab_wbc * lab_basophils/100, lab_d_dimer = lab_d_dimer/1000)
bel_lab_rev <- bel_lab_rev %>% mutate(lab_fibrinogen = lab_fibrinogen/100)
a <- bel_lab_rev %>% select(any_of(c(identification, demographics,blood_values,hospital_info, "study"))) %>% mutate_at("anonymized_patient_id",as.character)%>%
  mutate_at(vars("age_at_diagnosis","highest_who_score", "highest_respiratory_support"),as.numeric) %>% mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_who_score", "highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))
 
#the number of measurements per biomarker
df <- a %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(anonymized_patient_id) 
ggplot(df,aes(y=value,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#the mean of the number of measurements per biomarker
df <- a %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))
ggplot(df,aes(y=mean,x=name)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + theme_bw() + ylab("Average number of non-missing measurement per individual") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

df <- a %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>%
  filter(lab_result_date <= 30 & lab_result_date >= -2) %>% group_by(anonymized_patient_id) %>% 
  summarize_all(~sum(!is.na(.))) %>% 
  pivot_longer(!c(anonymized_patient_id,study,age_at_diagnosis,sex,lab_result_date,highest_respiratory_support,highest_who_score)) %>%
  group_by(name) %>% summarize(mean=mean(value))

df <- a %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. < 0,NA,.))) %>% 
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

saveRDS(a, file="belgium_lab.rds")
