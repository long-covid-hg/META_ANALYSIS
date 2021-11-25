### Association between genetic variants and mortality

path <- "/home/aganna/"

# Read in demograhics 
demog <- read_tsv(paste0(path,"curated_clinical/demog_data.tsv"))

# Read in genetics
genes <- read_tsv(paste0(path,"all_variant_rs35081325_rs11385942.tsv")) 

# Read in mortality
hosp <- read_tsv(paste0(path,"curated_clinical/hosp_data.tsv"))

# Combine
hosp <- hosp %>% select("anonymized_patient_id","death","hospitalization_end_cause","hospitalization","highest_who_score","hosp_thrombo","hosp_stroke","hosp_infarction","hosp_renal","hosp_hepatic","hosp_bleeding","hospitalization_duration","highest_respiratory_support")

all <- demog %>% select("anonymized_patient_id","age_at_diagnosis","sex") %>% left_join(genes,by="anonymized_patient_id") %>% inner_join(hosp,by="anonymized_patient_id") %>% mutate(death=ifelse(!is.na(death),death,ifelse(hospitalization_end_cause==1,1,0)), ventilation=ifelse(highest_respiratory_support %in% c(1,2),1,0)) 



### Association with ventilation
summary(glm(ventialtion~age_at_diagnosis+sex+study+I(2-rs35081325)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=all[all$pop=="EUR",], family=binomial()))
summary(glm(ventilation~age_at_diagnosis+sex+study+I(2-rs11385942)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=all[all$pop=="EUR",], family=binomial()))

sum(all$ventilation==1 & all$pop=="EUR",na.rm=T)
sum(all$ventilation==0 & all$pop=="EUR",na.rm=T)

### Association with death
summary(glm(death~age_at_diagnosis+sex+study+I(2-rs35081325)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=all[all$pop=="EUR",], family=binomial()))
summary(glm(death~age_at_diagnosis+sex+study+I(2-rs11385942)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=all[all$pop=="EUR",], family=binomial()))

sum(all$death==1 & all$pop=="EUR",na.rm=T)
sum(all$death==0 & all$pop=="EUR",na.rm=T)

### test different complications

summary(glm(hosp_hepatic~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_hepatic==1 & all$hospitalization==1,na.rm=T)

summary(glm(hosp_hepatic~age_at_diagnosis+I(2-rs11385942),data=all[all$hospitalization==1 & all$study.x=="swe",], family = binomial()))
summary(glm(hosp_hepatic~age_at_diagnosis+sex+I(2-rs11385942),data=all[all$hospitalization==1 & all$study.x=="belgium",], family = binomial()))
summary(glm(hosp_hepatic~age_at_diagnosis+sex+I(2-rs11385942),data=all[all$hospitalization==1 & all$study.x=="italy",], family = binomial()))
summary(glm(hosp_hepatic~age_at_diagnosis+sex+I(2-rs11385942),data=all[all$hospitalization==1 & all$study.x=="spain_bujanda",], family = binomial()))


summary(glm(hosp_thrombo~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_thrombo==1 & all$hospitalization==1,na.rm=T)

summary(glm(hosp_stroke~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_stroke==1 & all$hospitalization==1,na.rm=T)

summary(glm(hosp_infarction~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_infarction==1 & all$hospitalization==1,na.rm=T)

summary(glm(hosp_renal~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_renal==1 & all$hospitalization==1,na.rm=T)

summary(glm(hosp_bleeding~age_at_diagnosis+sex+study+I(2-rs11385942),data=all[all$hospitalization==1,], family = binomial()))
sum(all$hosp_bleeding==1 & all$hospitalization==1,na.rm=T)


summary(glm(hosp_hepatic~age_at_diagnosis+sex+study.x+I(2-rs17078348)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=all[all$pop=="EUR" & all$hospitalization==1,], family = binomial()))
sum(all$hosp_hepatic==1 & all$hospitalization==1 & all$pop=="EUR",na.rm=T)




