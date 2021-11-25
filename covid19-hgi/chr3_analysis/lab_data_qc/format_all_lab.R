setwd("/home/tomoko/02.associations")

a <- readRDS("../01.data.QC/hgi_belgium/belgium_lab.rds")
b <- readRDS("../01.data.QC/hgi_sweden/swe_lab.rds")
c <- readRDS("../01.data.QC/hgi_spain_bujanda/spain_bujanda_lab.rds")
d <- readRDS("../01.data.QC/hgi_canada/canada_lab.rds")

path <- "../covid19-hgi-clinical-values/"

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

## READ IN SPANISH BUTTI DATA
spain_butti_manifest <- read.xlsx(paste0(path,"hgi_spain_butti/Data dictionary HGI Covid-19 08092020 finished_cleaned_AG.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
spain_butti_lab <- read.xlsx(paste0(path,"hgi_spain_butti/Data dictionary HGI Covid-19 08092020 finished_cleaned_AG.xlsx"),sheetName = "Example_visit_time",stringsAsFactors=FALSE)

spain_butti_manifest_rev <- spain_butti_manifest %>% select(any_of(c(identification,hospital_info,"hospitalization_start","covid19_test_date")))
spain_butti_lab_rev <- merge(spain_butti_lab, spain_butti_manifest_rev, by="anonymized_patient_id", all=T)
spain_butti_lab_rev$hospitalization_start <- parse_date_covid19(spain_butti_lab_rev$hospitalization_start,"dmy")
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_result_date = as.numeric(as.Date(l.b_result_d.te) - as.Date(hospitalization_start) + covid19_test_date))
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(study="spain_butti")

spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_wbc = lab_leukocytes/1000)
spain_butti_lab_rev <- spain_butti_lab_rev %>% mutate(lab_neutrophils = as.numeric(lab_neutrophils)*lab_wbc, lab_lymphocytes = as.numeric(lab_lymphocytes)*lab_wbc,
                                                      lab_monocytes = as.numeric(lab_monocytes)*lab_wbc, lab_eosinophils=as.numeric(lab_eosinophils)*lab_wbc,
                                                      lab_basophils = as.numeric(lab_basophils)*lab_wbc)

e <- spain_butti_lab_rev %>% mutate(study="spain_butti") %>% select(any_of(c(identification, demographics,blood_values,hospital_info,"study"))) %>% mutate_at("anonymized_patient_id",as.character) %>% 
  mutate_at(vars(colnames(spain_butti_lab_rev)[colnames(spain_butti_lab_rev) %in% blood_values[-1]], "highest_respiratory_support"),as.numeric)
e <- e %>% mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))


## READ IN ITALIAN DATA
ita_renieri_manifest <- fread("../covid19-hgi-clinical-values/hgi_italy/Italy.withCom_20201008TN.tsv")
ita_renieri_lab <- fread("../covid19-hgi-clinical-values/hgi_italy/Italy.lab_20201008TN.tsv")

ita_renieri_lab_rev <- merge(ita_renieri_lab, ita_renieri_manifest, by="anonymized_patient_id", all=T)

f <- ita_renieri_lab_rev %>% mutate(study="italy") %>% select(any_of(c(identification, demographics,hospital_info,blood_values, "study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support", "highest_who_score")), .funs = funs(ifelse(. == -1,NA,.)))

## READ IN BRAZIL DATA
bra_manifest <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
bra_lab <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_visit_time",stringsAsFactors=FALSE,colClasses="character")

brazil_lab <- bra_lab %>% mutate(lab_wbc = lab_leukocytes/1000, lab_platelets = lab_platelets/1000)
brazil_lab <- brazil_lab %>% mutate(lab_lymphocytes = lab_lymphocytes*lab_wbc, lab_neutrophils = lab_neutrophils*lab_wbc,
                                    lab_monocytes = lab_monocytes*lab_wbc, lab_eosinophils = lab_eosinophils**lab_wbc, 
                                    lab_basophils = lab_basophils*lab_wbc)

brazil_lab_rev <- merge(brazil_lab, bra_manifest, by="anonymized_patient_id", all=T)
g <- brazil_lab_rev %>% mutate(study="brazil") %>% select(any_of(c(identification, demographics,hospital_info,blood_values,"study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character)
g <- g %>% group_by(anonymized_patient_id) %>% mutate(first_date=min(lab_result_date)) %>% ungroup() %>%
  mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - first_date))

g <- g %>% mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support")), .funs = funs(ifelse(. == -1,NA,.)))

## READ IN GERMANY SCHULTE
library(readxlsb)
germany_schulte_manifest <- read_xlsb(paste0("../covid19-hgi-clinical-values/hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 2)
germany_schulte_lab <- read_xlsb(paste0("../covid19-hgi-clinical-values/hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 3)

scale <- read_xlsb(paste0("../covid19-hgi-clinical-values/hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 4)
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
  
tmp <- merge(germany_schulte_labrev, germany_schulte_manifest, by="anonymized_patient_id", all=T)
h <- tmp %>% mutate(lab_result_date = as.numeric(as.Date(lab_result_date, origin="1900-01-01") - as.Date(covid19_test_date))) %>% 
  select(any_of(c(identification, demographics,hospital_info,blood_values,"study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(c("highest_respiratory_support","highest_who_score")), .funs = funs(ifelse(. == -1,NA,.)))

## READ IN ITALY VALENTI
ita_valenti_manifest <- read.xlsx(paste0(path,"hgi_italy_valenti/Clinical_Data_V2_Milan.xlsx"),sheetName = "One_Time",stringsAsFactors=FALSE,colClasses="character")
ita_valenti_lab <- read.xlsx(paste0(path,"hgi_italy_valenti/Clinical_Data_V2_Milan.xlsx"),sheetIndex = 1, stringsAsFactors=FALSE,colClasses="character")

tmp <- merge(ita_valenti_lab, ita_valenti_manifest, by="anonymized_patient_id", all=T)
i <- tmp %>% group_by(anonymized_patient_id) %>% mutate(first_date=min(lab_result_date), study="italy_valenti") %>% ungroup() %>%
  mutate(lab_result_date = as.numeric(as.Date(lab_result_date) - as.Date(first_date))) %>% 
  select(any_of(c(identification, demographics,hospital_info,blood_values,"study"))) %>% mutate_at(vars(c("anonymized_patient_id","sex")),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values)), .funs = funs(ifelse(. < 0,NA,.))) %>%
  mutate_at(.vars = vars(any_of(c("highest_respiratory_support","highest_who_score"))), .funs = funs(ifelse(. == -1,NA,.)))

lab_data <- bind_rows(list(a,b,c,d,e,f,g,h,i))

saveRDS(lab_data, file="all_lab_data.rds")
