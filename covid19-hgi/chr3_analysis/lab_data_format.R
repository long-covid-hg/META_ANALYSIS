library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)
library(readr)

parse_date_covid19 <- function(date,order_type)
{
  date_out <- as_date(rep(NA,length(date)))
  for (i in 1:length(date))
  {
    i <- as.character(date[i])
    if(i=="-1" | is.na(i) | i=="1899-12-29" | i=="1899-12-30") {
      o <- NA}
    else if (nchar(i)==5) {
      o <- as_date(as.numeric(i),origin="1899-12-30") 
    }
    else if (nchar(i)>5) {
      o <- as_date(parse_date_time(i,order_type))
    }
    else {o <- NA}
    date_out[which(i==date)] <- o
  }
  return(date_out)
}


path <- "../covid19-hgi-clinical-values/"

## READ IN BELGIUM DATA
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-250920_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)
bel_manifest <- bel_manifest[-dim(bel_manifest)[1],]
bel_lab <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-250920_TN.xls"),sheetName = "visit_time",stringsAsFactors=FALSE)
bel_lab <- bel_lab[!is.na(bel_lab$anonymized_patient_id),]
bel_lab <- bel_lab[-dim(bel_lab)[1],]

## READ IN SWEDISH DATA
swe_manifest <- read.xlsx(paste0(path,"hgi_swe/PronMed genetic fenotypes.xlsx"),sheetName = "Data",stringsAsFactors=FALSE)
swe_lab <- read.xlsx(paste0(path,"hgi_swe/PronMed blood chemistry v2_editedAG.xlsx"),sheetName = "values",stringsAsFactors=FALSE)
swe_lab_plus <- read.xlsx(paste0(path,"hgi_swe/PronMedBloodChemistry_v2_editedTN.xlsx"),sheetIndex = "New Labs",stringsAsFactors=FALSE)

swe_lab_rev <- merge(swe_lab, swe_lab_plus, by=c("anonymized_patient_id", "lab_result_date.relative.ICU.admission"), all = T)

## READ IN SPANISH BUJANDA DATA
spain_bujanda_manifest <- read.xlsx(paste0(path,"hgi_spain_bujanda/international database_Spanish_Bujanda_Bañales_Nafría_v5_editedTN.xlsx"),sheetName = "Spanish patients_One-visit",stringsAsFactors=FALSE,colClasses="character")
spain_bujanda_lab <- read.xlsx(paste0(path,"hgi_spain_bujanda/international database_Spanish_Bujanda_Bañales_Nafría_v5_editedTN.xlsx"),sheetName = "Laboratory data_time",stringsAsFactors=FALSE)

## READ IN SPANISH BUTTI DATA
spain_butti_manifest <- read.xlsx(paste0(path,"hgi_spain_butti/Data dictionary HGI Covid-19 08092020 finished_cleaned_AG.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
spain_butti_lab <- read.xlsx(paste0(path,"hgi_spain_butti/Data dictionary HGI Covid-19 08092020 finished_cleaned_AG.xlsx"),sheetName = "Example_visit_time",stringsAsFactors=FALSE)

## READ IN ITALIAN DATA
italy_manifest <- fread("../covid19-hgi-clinical-values/hgi_italy/Italy.withCom_20201008TN.tsv")
italy_lab <- fread("../covid19-hgi-clinical-values/hgi_italy/Italy.lab_20201008TN.tsv")

## READ IN BRAZIL DATA
brazil_manifest <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
brazil_lab <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_visit_time",stringsAsFactors=FALSE,colClasses="character")

## READ IN CANADA DATA
canada_manifest <- fread("../covid19-hgi-clinical-values/hgi_canada/bqc19_one_time_V2.tsv")
canada_lab <- fread("../covid19-hgi-clinical-values/hgi_canada/bqc19_visit_time_V2.tsv")

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
                  "lab_nk", "lab_n_l_ratio")

## START PROCESSING HOSPITAL INFORMATION

### HARMONIZING BELGIUM DATA
a <- bel_manifest %>% select(any_of(c(identification,demographics,hospital_info,covid19))) %>% 
  mutate(study="belgium", covid19_test_date=as.Date(covid19_test_date,"%d/%m/%Y")) 

bel_lab_rev <- merge(bel_lab, a, by="anonymized_patient_id", all.x=T)
bel_lab_rev <- bel_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date,"%d/%m/%Y") - as.Date(covid19_test_date))) %>% 
  rename(lab_aptt = APTT, lab_inr = INR)

set_numeric <- function(x){as.numeric(x)}
bel_lab_rev <- bel_lab_rev %>% mutate_at(vars(c(any_of(blood_values),-lab_result_date)), set_numeric)
bel_lab_rev <- bel_lab_rev %>% mutate(lab_lymphocytes = lab_wbc * lab_lymphocytes/100, lab_neutrophils = lab_wbc * lab_neutrophils/100,
                                      lab_monocytes = lab_wbc * lab_monocytes/100, lab_eosinophils = lab_wbc * lab_eosinophils/100,
                                      lab_basophils = lab_wbc * lab_basophils/100, lab_d_dimer = lab_d_dimer/1000)
a <- bel_lab_rev %>% select(any_of(c(identification, demographics,blood_values,hospital_info, "study"))) %>% mutate_at("anonymized_patient_id",as.character)%>%
  mutate_at(vars("age_at_diagnosis","highest_who_score", "highest_respiratory_support"),as.numeric)

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

swe_lab_rev <- merge(swe_lab, swe_manifest, by="anonymized_patient_id", all.x=T)

swe_lab <- swe_lab_rev %>% select(any_of(c(identification, blood_values, demographics,hospital_info,covid19,"study")))

#write.table(swe_lab, "swe_lab_cleanedTN20200925.tsv", col.names = T, row.names = F, sep="\t", quote = F)

b <- swe_lab %>% select(any_of(c(identification, demographics,blood_values,hospital_info,"study"))) %>% mutate_at(vars("anonymized_patient_id","sex"),as.character)

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
                                                  lab_ck = lab_ck..IU.L., lab_fibrinogen = lab_fibrinogen..g.L., lab_creatinine = lab_creatinine..mg.dL.)

spain_bujanda_lab_rev <- merge(spain_bujanda_lab, spain_bujanda_manifest, by="anonymized_patient_id", all.x=T)

spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - as.Date(rt.PCR.date)))

spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% mutate(study = "spain_bujanda")
spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% rename(lab_aptt = lab_appt, lab_inr = lab_inr.)
spain_bujanda_lab_rev <- spain_bujanda_lab_rev %>% rename(covid19_test_date = rt.PCR.date) %>% mutate_at(vars(date_of_death), function(x){format(as.Date(x, origin = "1899-12-30"))})


c <- spain_bujanda_lab_rev %>% select(any_of(c(identification, demographics,blood_values,hospital_info, "study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(vars(colnames(spain_bujanda_lab_rev)[colnames(spain_bujanda_lab_rev) %in% blood_values[-1]], "highest_respiratory_support"),as.numeric)


#write.table(unique(spain_bujanda_lab_rev$anonymized_patient_id[spain_bujanda_lab_rev$lab_result_date < -2]), file="spain_bujanda_confirm_covid_date_id.txt", quote = F, col.names = F, row.names = F)
## spain_butti_manifest
spain_butti_manifest_rev <- spain_butti_manifest %>% select(any_of(c(identification,hospital_info,"hospitalization_start","covid19_test_date")))
spain_butti_lab_rev <- merge(spain_butti_lab, spain_butti_manifest_rev, by="anonymized_patient_id", all.x=T)
spain_butti_lab_rev$hospitalization_start <- parse_date_covid19(spain_butti_lab_rev$hospitalization_start,"dmy")
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(l.b_result_d.te) - as.Date(hospitalization_start) + covid19_test_date))
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(study="spain_butti")

spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_wbc = lab_leukocytes/1000)
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_neutrophils = as.numeric(lab_neutrophils)*lab_wbc, lab_lymphocytes = as.numeric(lab_lymphocytes)*lab_wbc,
                                                      lab_monocytes = as.numeric(lab_monocytes)*lab_wbc, lab_eosinophils=as.numeric(lab_eosinophils)*lab_wbc,
                                                      lab_basophils = as.numeric(lab_basophils)*lab_wbc)


d <- spain_butti_lab_rev %>% mutate(study="spain_butti") %>% select(any_of(c(identification, demographics,blood_values,hospital_info,"study"))) %>% mutate_at("anonymized_patient_id",as.character) %>% 
  mutate_at(vars(colnames(spain_butti_lab_rev)[colnames(spain_butti_lab_rev) %in% blood_values[-1]], "highest_respiratory_support"),as.numeric)

### HARMONIZING ITALIAN DATA
italy_lab_rev <- merge(italy_lab, italy_manifest, by="anonymized_patient_id", all.x=T)

e <- italy_lab_rev %>% mutate(study="italy") %>% select(any_of(c(identification, demographics,hospital_info,blood_values, "study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character)

### HARMONIZING BRAZIL DATA
brazil_lab <- brazil_lab %>% mutate(lab_wbc = lab_leukocytes/1000, lab_platelets = lab_platelets/1000)
brazil_lab <- brazil_lab %>% mutate(lab_lymphocytes = lab_lymphocytes*lab_wbc, lab_neutrophils = lab_neutrophils*lab_wbc,
                                    lab_monocytes = lab_monocytes*lab_wbc, lab_eosinophils = lab_eosinophils**lab_wbc, 
                                    lab_basophils = lab_basophils*lab_wbc)

brazil_lab_rev <- merge(brazil_lab, brazil_manifest, by="anonymized_patient_id", all.x=T)
f <- brazil_lab_rev %>% mutate(study="brazil") %>% select(any_of(c(identification, demographics,hospital_info,blood_values,"study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character)
f <- f %>% group_by(anonymized_patient_id) %>% mutate(first_date=min(lab_result_date)) %>% ungroup() %>%
  mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - first_date))

### HARMONIZING CANADA DATA
tmp <- merge(canada_lab, canada_manifest, by="anonymized_patient_id", all.x=T)
tmp <- tmp %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)))
g <- tmp %>% mutate(study="canada") %>% select(any_of(c(identification, demographics,hospital_info,blood_values, "study"))) %>% mutate_at("anonymized_patient_id",as.character) %>% 
  mutate(sex = ifelse(sex=="F", as.character(1), as.character(0)))

g <- g %>% mutate(lab_d_dimer = lab_d_dimer/1000)

lab_data <- bind_rows(list(a,b,c,d,e,f,g))

write.table(lab_data, file="All.cohort.lab.tsv", col.names = T, row.names = F, sep="\t", quote = F)
