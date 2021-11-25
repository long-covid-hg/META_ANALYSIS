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



path <- "/home/aganna/"

## READ IN BELGIUM DATA
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-111020_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)


## READ IN SWEDISH DATA
swe_manifest <- read.xlsx(paste0(path,"hgi_swe/PronMed genetic fenotypes with WHO.xlsx"),sheetName = "Data",stringsAsFactors=FALSE)


## READ IN SPANISH BUJANDA DATA
spain_bujanda_manifest <- read.xlsx(paste0(path,"hgi_spain_bujanda/international_database_Spanish_Bujanda_v5_editedTN.xlsx"),sheetName = "Spanish patients_One-visit",stringsAsFactors=FALSE,colClasses="character")
spain_bujanda_manifest <- spain_bujanda_manifest %>% filter(!is.na(anonymized_patient_id))


## READ IN SPANISH BUTTI DATA
spain_butti_manifest <- read.xlsx(paste0(path,"hgi_spain_butti/Data_dictionary_HGI_Covid-19_V2_editedTN_Corrections.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")


## READ IN ITALIAN DATA
ita_renieri_manifest <- read_tsv(paste0(path,"hgi_italy/Italy.withCom_20201008TN.tsv"))

## READ BRASIL
bra_manifest <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")

## READ CANADA
can_manifest <- read_tsv(paste0(path,"hgi_canada/bqc19_one_time_V2.tsv"))

## READ ITALY VALENTI
italy_valenti_manifest <- read.xlsx(paste0(path,"hgi_italy_valenti/Form_GWAS_COVID_DataDictConv_Ultima_13OCT20_mod_FM_09102020.xlsx"),sheetName = "One_Time",stringsAsFactors=FALSE,colClasses="character")

## READ IN GERMAN_SCHULTE
germany_schulte_manifest <- read.xlsx(paste0(path,"hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xls"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")



## Variable name dictionary
identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex","ancestry","height","weight","smoking")

comorbidities <- 
c("com_hiv",
"com_immunocomp",
"com_transplant",
"com_autoimm_rheum",
"com_type_i_diabetes",
"com_type_ii_diabetes",
"com_asthma",
"com_chronic_pulm",
"com_sleep_apnea",
"com_liver",
"com_gallbl",
"com_pancreas",
"com_chronic_kidney",
"com_heart_failure",
"com_hypertension",
"com_infarction",
"com_vascular",
"com_stroke",
"com_dementia",
"com_neurological",
"com_leukemia",
"com_lymphoma",
"com_malignant_solid",
"com_dialysis",
"com_afib",
"com_dialysis")


hospital_info <- 
c("hospitalization",
"hospitalization_start",
"hospitalization_end",
"hospitalization_end_cause",
"icu_admit",
"icu_duration",
"highest_respiratory_support",
"days_ventilator",
"hosp_dvt",
"hosp_thrombo",
"hosp_stroke",
"hosp_infarction",
"hosp_renal",
"hosp_hepatic",
"hosp_bleeding",
"death",
"cause_of_death",
"date_of_death",
"highest_who_score")


covid19 <-
c("covid19_test",
"covid19_test_date",
"covid19_test_type",
"covid19_first_symptoms_date")

medications <- 
c("steroids",
"biologics",
"lmwh",
"hcq",
"remdesivir")

## START PROCESSING HOSPITAL INFORMATION

## Demographics ##
a <- bel_manifest %>% select(any_of(c(identification,demographics)))  %>% mutate_at(vars(-("anonymized_patient_id")),as.numeric) %>% mutate(study="belgium") 
b <- swe_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(study="sweden") %>% mutate_at("anonymized_patient_id",as.character)
c <- spain_bujanda_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(weight=as.numeric(weight),study="spain_bujanda") %>% mutate_at("anonymized_patient_id",as.character)
d <- spain_butti_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(sex=as.numeric(sex),age_at_diagnosis=as.numeric(age_at_diagnosis),study="spain_butti") %>% mutate_at("anonymized_patient_id",as.character)
e <-  ita_renieri_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(study="italy_renieri") 
f <- bra_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(study="brasil") %>% mutate_at("anonymized_patient_id",as.character)
g <- can_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(study="canada",sex=ifelse(sex=="M",0,1)) %>% mutate_at("anonymized_patient_id",as.character)
h <- italy_valenti_manifest %>% select(any_of(c(identification,demographics))) %>% mutate(study="italy_valenti") %>% mutate_at(c("sex","ancestry"),as.numeric)

# Combine
demog_data <- bind_rows(list(a,b,c,d,e,f,g,h))

# Recode missing and fix small issues (height==1 in the swedish data and height or weight=0 in the spain_butti data)
demog_data <- demog_data %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. == -1,NA,.))) %>% mutate(height=ifelse(height==1 | height==0 | height==16,NA,height),weight=ifelse(weight==0,NA,height))

# Plotting - missing %
df <- demog_data %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name, ncol=3) + theme_bw() + ylab("% missing") + ylim(0,100) + theme(axis.text.x=element_text(angle=45, hjust=1)) 

# Plotting - continuos
df <- demog_data %>% select(anonymized_patient_id,age_at_diagnosis,height,weight,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name,scales = "free_y") + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) 

# Plotting - dicotom
df <- demog_data %>% 
mutate(ancestry2=ifelse(ancestry==0,0,1)) %>%
select(anonymized_patient_id,sex,ancestry2,study) %>%
pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name,scales = "free_y") + theme_bw() + ylab("%") + theme(axis.text.x=element_text(angle=45, hjust=1)) 

# Export
demog_data %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"curated_clinical/demog_data.tsv"),col_names = TRUE)


## Comorbidities ##
a <- bel_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate_at(vars(-("anonymized_patient_id")),as.numeric) %>% mutate(study="belgium")
b <- swe_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="sweden") %>% mutate_at("anonymized_patient_id",as.character)
c <- spain_bujanda_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="spain_bujanda") %>% mutate_at("anonymized_patient_id",as.character)
d <- spain_butti_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="spain_butti") %>% mutate_at("anonymized_patient_id",as.character)
e <- ita_renieri_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="italy_renieri")
f <- bra_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="brasil") %>% mutate_at("anonymized_patient_id",as.character)
g <- can_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="canada") %>% mutate_at("anonymized_patient_id",as.character)
h <- italy_valenti_manifest %>% select(any_of(c(identification,comorbidities))) %>% mutate(study="italy_valenti")

# Combine
comorb_data <- bind_rows(list(a,b,c,d,e,f,g,h))

# Recode missing and fix small issues
comorb_data <- comorb_data %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. == -1,NA,.)))


# Plotting - missing %
df <- comorb_data %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Exclude variables with too many missing
comorb_data_miss_cleaned <- comorb_data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=80)

# Plotting - dicotom
df <- comorb_data_miss_cleaned %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100)) %>% filter(value!=0)
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% with the comorbidity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
comorb_data %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"curated_clinical/comorb_data.tsv"),col_names = TRUE)




## COVID19 info ##
a <- bel_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="belgium", covid19_test_date=parse_date_covid19(covid19_test_date,"dmy"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"dmy"))

b <- swe_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="sweden", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

c <- spain_bujanda_manifest %>% rename("covid19_test_date"="rt.PCR.date") %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="spain_bujanda", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

d <- spain_butti_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="spain_butti", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

e <- ita_renieri_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="italy_renieri", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))

f <- bra_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="brasil", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

g <- can_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="canada", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

h <- italy_valenti_manifest %>% select(any_of(c(identification,covid19))) %>% 
  mutate(study="italy_valenti", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))


# Combine
covid_data <- bind_rows(list(a,b,c,d,e,f,g,h))

# Recode missing and fix small issues - 
# exclude dates in the future (spain_butti) 
covid_data <- covid_data %>% 
  mutate_at(.vars = vars(-c(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date)), .funs = funs(ifelse(. == -1,NA,.))) #%>% filter(covid19_first_symptoms_date < "2020-08-31" | is.na(covid19_first_symptoms_date))


# Plotting - missing %
df <- covid_data %>% select(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date,study) %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting - dates
df <- covid_data %>% select(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
covid_data %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"curated_clinical/covid_data.tsv"),col_names = TRUE)




## hospital_info ##
a <- bel_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
mutate(study="belgium", hospitalization_start=parse_date_covid19(hospitalization_start,"dmy"),hospitalization_end=parse_date_covid19(hospitalization_end,"dmy"), date_of_death=parse_date_covid19(date_of_death,"dmy"), highest_who_score=as.numeric(highest_who_score)) 

b <- swe_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
mutate(study="sweden", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

c <- spain_bujanda_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
mutate(study="spain_bujanda", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"),date_of_death=parse_date_covid19(date_of_death,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

d <- spain_butti_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="spain_butti", hospitalization_start=parse_date_covid19(hospitalization_start,"dmy"),hospitalization_end=parse_date_covid19(hospitalization_end,"dmy"),highest_respiratory_support=as.numeric(highest_respiratory_support)) %>% mutate_at("anonymized_patient_id",as.character)

e <- ita_renieri_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="italy_renieri", hospitalization_start=parse_date_covid19(hospitalization_start,"dmy"),hospitalization_end=parse_date_covid19(hospitalization_end,"dmy"),date_of_death=parse_date_covid19(date_of_death,"dmy"))

f <- bra_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="brasil", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

g <- can_manifest %>% select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="canada", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"),date_of_death=parse_date_covid19(date_of_death,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

h <- italy_valenti_manifest %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="italy_valenti", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd")) %>% mutate_at("anonymized_patient_id",as.character) %>% mutate_at("icu_admit",as.numeric)

# Combine
hosp_data <- bind_rows(list(a,b,c,d,e,f,g,h))

# Recode missing and fix small issues - 
# add hospitalization duration - 
# remove negative hospitalization duration
# create consistency in the highest respirator support -- still to do
hosp_data <- hosp_data %>% 
mutate_at(.vars = vars(-c(anonymized_patient_id,hospitalization_start,hospitalization_end,date_of_death)), .funs = funs(ifelse(. == -1,NA,.))) %>%
mutate(hospitalization_duration=as.numeric(hospitalization_end-hospitalization_start)) #%>% 
#filter(hospitalization_duration>=0 | is.na(hospitalization_duration))

#hosp_data %>%  select(highest_respiratory_support,days_ventilator)
#hosp_data[hosp_data$hospitalization_duration < 0 | hosp_data$hospitalization_end > today(),]
# Check those with ICU duration if they all have icu admission
#table(is.na(hosp_data$icu_duration), is.na(hosp_data$icu_admit))

# Plotting - missing %
df <- hosp_data  %>% mutate(hospitalization_start=ifelse(hospitalization_start==-1,-1,0),hospitalization_end=ifelse(hospitalization_end==-1,-1,0)) %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Exclude variables with too many missing
#hosp_data_data_miss_cleaned <- hosp_data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=80)

# Plotting - dates
df <- hosp_data %>% select(anonymized_patient_id,hospitalization_start,hospitalization_end,date_of_death,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plotting - durations
df <- hosp_data %>% select(anonymized_patient_id,hospitalization_duration,icu_duration,days_ventilator,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name,scales = "free_y") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plotting - dicotom
df <- hosp_data %>% select(anonymized_patient_id,matches("hosp_"),hospitalization,icu_admit,study) %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100)) %>% filter(value!=0)
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("%") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
hosp_data %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"curated_clinical/hosp_data.tsv"),col_names = TRUE)




## Medications ##
a <- bel_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="belgium")
b <- swe_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="sweden") %>% mutate_at("anonymized_patient_id",as.character)
c <- spain_bujanda_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="spain_bujanda") %>% mutate_at("anonymized_patient_id",as.character)
d <- spain_butti_manifest %>% select(any_of(c(identification,medications))) %>% mutate(remdesivir=as.numeric(remdesivir),study="spain_butti") %>% mutate_at("anonymized_patient_id",as.character)
e <- ita_renieri_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="italy_renieri")
f <- bra_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="brasil") %>% mutate_at("anonymized_patient_id",as.character)
g <- can_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="canada") %>% mutate_at("anonymized_patient_id",as.character)
h <- italy_valenti_manifest %>% select(any_of(c(identification,medications))) %>% mutate(study="italy_valenti")


# Combine
medication_data <- bind_rows(list(a,b,c,d,e,f,g,h))

# Recode missing and fix small issues
medication_data <- medication_data %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. == -1,NA,.)))


# Plotting - missing %
df <- medication_data %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name, ncol=5) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting -  distribution
df <- medication_data %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100)) %>% filter(value!=0)
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name,ncol=5) + theme_bw() + ylab("% takinf the medication") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
medication_data %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"curated_clinical/medication_data.tsv"),col_names = TRUE)

