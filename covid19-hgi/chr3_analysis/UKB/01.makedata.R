setwd("/project/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/06.clinical_values")

data <- readRDS("../../../scratch/data/01.UKBB/UKBall20201023.457941.rds")
data <- data %>% filter(status == 1)

identification <- c("anonymized_patient_id")
demographics <- c("age_at_diagnosis", "sex", "ancestry", "height", "weight", "smoking")
comorbidities <- c("com_hiv",
                   "com_immunocomp",
                   "com_transplant",
                   "com_autoimm_rheum",
                   "com_type_i_diabetes",
                   "com_type_ii_diabetes",
                   "com_diabetes",
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

hospital <- c("hospitalization",
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

covid19 <- c("covid19_test",
             "covid19_test_date",
             "covid19_test_type",
             "covid19_first_symptoms_date")

rx <- c("steroids",
        "biologics",
        "lmwh",
        "hcq",
        "remdesivir")

library(dplyr)
library(tidyr)

data <- data %>% rename(anonymized_patient_id = ID, age_at_diagnosis = AGE) %>% mutate(sex = ifelse(SEX == 0, 1, 0))

data <- data %>% mutate(death = ifelse(!is.na(date_of_death) & as.Date(date_of_death) < as.Date(COVIDdate) + 30, 1, 0))
data$date_of_death[data$death == 0] <- NA

d <- fread("../../../scratch/data/01.UKBB/death_20201021.txt.gz", sep="\t")
unique(d$date_of_death[grepl("09\\/2020$", d$date_of_death)])
d1 <- fread("../../../scratch/data/01.UKBB/death_cause_20201021.txt.gz", sep="\t")
D <- merge(d, d1, by=c("eid", "ins_index"))
colnames(D)[1] <- "ID"
D1 <- D[as.Date(D$date_of_death, "%d/%m/%Y") >= "2020/03/16",]
D1 <- D1 %>% rename(anonymized_patient_id = ID)

data1 <- data %>% merge(D1[,c("anonymized_patient_id", "cause_icd10")], by="anonymized_patient_id", all.x=T)

data2 <- data1 %>% group_by(anonymized_patient_id) %>% mutate(cause_of_death = case_when(any(cause_icd10 == "U071") ~ 1,
                                                                                         death == 1 ~ 2,
                                                                                         TRUE ~ 0)) %>%
  ungroup() %>% distinct(anonymized_patient_id, .keep_all=TRUE) %>% select(c(-cause_icd10,-SEX,-CENTRE)) %>% mutate(date_of_death = as.Date(date_of_death, "%d/%m/%Y", origin="01/01/1970"))

data <- data2 %>% rename(covid19_test_date = COVIDdate, covid19_test = status)

#critical care
c <- fread("../../../scratch/data/01.UKBB/hesin_critical_20200920.txt.gz", sep="\t")
c <- c[c$eid %in% data$anonymized_patient_id]
c <- c %>% rename(anonymized_patient_id = eid)
data1 <- data %>% merge(c, by="anonymized_patient_id", all.x=T)
data1 <- data1 %>% filter(as.Date(ccstartdate, "%d/%m/%Y") < 30 + as.Date(covid19_test_date))
data1 <- data1 %>% filter(as.Date(ccstartdate, "%d/%m/%Y") > as.Date(covid19_test_date) - 14)
data1 <- data1 %>% select(c("anonymized_patient_id", "bressupdays", "aressupdays", "acardsupdays", "rensupdays"))

data2 <- data %>% merge(data1, by="anonymized_patient_id", all.x=T) %>% group_by(anonymized_patient_id) %>%
  mutate_at(vars(c("bressupdays", "aressupdays", "acardsupdays", "rensupdays")), max, na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE)

#hesin 
hesin_diag <- fread("../../../scratch/data/01.UKBB/hesin_diag_20200910.txt.gz", sep="\t")
hesin <- fread("../../../scratch/data/01.UKBB/hesin_20200910.txt.gz", sep="\t")
hesin1 <- hesin[grepl("/2020$", hesin$epiend),]
#hesin2 <- hesin1[as.Date(hesin1$admidate, "%d/%m/%Y") >= as.Date("2020/03/16"),]
tmp <- hesin1[hesin1$eid == unique(data2$anonymized_patient_id)[1],]
TMP <- tmp
for(i in seq(2,length(unique(data2$anonymized_patient_id)))){
  tmp <- hesin1[hesin1$eid == unique(data2$anonymized_patient_id)[i],]
  TMP <- rbind(TMP,tmp)
}
hesin_new <- merge(TMP, hesin_diag, by=c("eid", "ins_index"))
hesin_oper <- fread("../../../scratch/data//01.UKBB/hesin_oper_20200910.txt.gz", sep="\t")
hesin_new1 <- merge(TMP, hesin_oper, by=c("eid", "ins_index"))

data2$resp_severe <- 0
##resp_severe
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 ,]
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "J80" %in% tmp1$diag_icd10] <- 1
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "J9600" %in% tmp1$diag_icd10] <- 1
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "J9609" %in% tmp1$diag_icd10] <- 1
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "Z991" %in% tmp1$diag_icd10] <- 1
}
for(i in seq(1,length(unique(hesin_new1$eid)))){
  tmp <- hesin_new1[hesin_new1$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 ,]
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new1$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "E851" %in% tmp1$oper4] <- 1
  data2$resp_severe[data2$anonymized_patient_id == unique(hesin_new1$eid)[i] & as.Date(data2$covid19_test_date) <= max(as.Date(tmp1$epiend, "%d/%m/%Y")) & "E852" %in% tmp1$oper4] <- 1
}
sum(data2$resp_severe)#143

data2 <- data2 %>% mutate(resp_severe = case_when(aressupdays > 0  & !is.infinite(aressupdays) ~ 1,
                                                  TRUE ~ resp_severe))
sum(data2$resp_severe)#156

##resp_mild
data2$resp_mild <- 0
data2 <- data2 %>% mutate(resp_mild = case_when(bressupdays > 0  & !is.infinite(bressupdays) ~ 1,
                                                resp_severe == 1 ~ 1,
                                                TRUE ~ 0))
sum(data2$resp_mild)#182

##hosp_dvt I81 I82* 
data2$hosp_dvt <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  data2$hosp_dvt[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & "I81" %in% tmp1$diag_icd10] <- 1
  data2$hosp_dvt[data2$anonymized_patient_id == unique(hesin_new$eid)[i] & any(grepl("^I82", tmp1$diag_icd10))] <- 1
}
sum(data2$hosp_dvt)#3

##hosp_thrombo I26*
data2$hosp_thrombo <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("^I26", tmp1$diag_icd10))){
    data2$hosp_thrombo[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_thrombo)#23

##hosp_stroke I63*, 
data2$hosp_stroke <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("^I63", tmp1$diag_icd10))){
    data2$hosp_stroke[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_stroke)#16

#hosp_infarction I21*
data2$hosp_infarction <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("^I21", tmp1$diag_icd10))){
    data2$hosp_infarction[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_infarction)#12

#hosp_renal N17*
data2$hosp_renal <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("^N17", tmp1$diag_icd10))){
    data2$hosp_renal[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_renal)#189
data2 <- data2 %>% mutate(hosp_renal = case_when(rensupdays > 0  & !is.infinite(rensupdays) ~ 1,
                                                TRUE ~ hosp_renal))
sum(data2$hosp_renal)#194
#hosp_bleeding I60*, I61*, I62*, K250 K252 K260 K262 K270 K272 K280 K282 K625 K922 I850 
data2$hosp_bleeding <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("^I60", tmp1$diag_icd10) | grepl("^I61", tmp1$diag_icd10) | grepl("^I62", tmp1$diag_icd10) |
         grepl("K250", tmp1$diag_icd10) | grepl("K252", tmp1$diag_icd10) | grepl("K260", tmp1$diag_icd10) |
         grepl("K262", tmp1$diag_icd10) | grepl("K272", tmp1$diag_icd10) | grepl("K270", tmp1$diag_icd10) |
         grepl("K280", tmp1$diag_icd10) | grepl("K282", tmp1$diag_icd10) | grepl("K625", tmp1$diag_icd10) | 
         grepl("K922", tmp1$diag_icd10) | grepl("I850", tmp1$diag_icd10))){
    data2$hosp_bleeding[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_bleeding)#20

#hosp_hepatic K720
data2$hosp_hepatic <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new1$eid)[i],]
  TMP <- data2[data2$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp[(as.Date(tmp$epistart, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14 | as.Date(tmp$admidate, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)) - 14) & as.Date(tmp$epiend, "%d/%m/%Y") >= min(as.Date(TMP$covid19_test_date)),]
  if(any(grepl("K720", tmp1$diag_icd10))){
    data2$hosp_hepatic[data2$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
}
sum(data2$hosp_hepatic)#1

final <- data2 %>% select(-c(bressupdays,aressupdays,acardsupdays,rensupdays))

saveRDS(final, file="../../../scratch/data/01.UKBB/clinical_value_20201024.rds")
