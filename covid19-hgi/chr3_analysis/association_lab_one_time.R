setwd("/home/tomoko/02.associations")
lab_data <- readRDS("all_lab_data.rds")

identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex")

hospital_info <- 
  c("highest_respiratory_support",
    "highest_who_score","resp_mild", "resp_severe")

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


lab_data <- lab_data %>% mutate(lab_n_l_ratio = case_when(is.na(lab_n_l_ratio) ~ lab_neutrophils/lab_lymphocytes,
                                                          TRUE ~ lab_n_l_ratio))

lab_data$resp_mild <- ifelse(lab_data$highest_respiratory_support >= 0, 1, 0)
lab_data$resp_severe <- ifelse(lab_data$highest_respiratory_support >= 1, 1, 0)

lab_data <- lab_data %>% mutate(anonymized_patient_id = case_when(study == "sweden" ~ paste0("SW_",anonymized_patient_id),
                                                                  study == "belgium" ~ paste0("BB_",anonymized_patient_id),
                                                                  study == "ita_renieri" ~ paste0("ITR_",anonymized_patient_id),
                                                                  study == "spain_bujanda" ~ paste0("SB_",anonymized_patient_id),
                                                                  study == "spain_butti" ~ paste0("SB_",anonymized_patient_id),
                                                                  study == "canada" ~ paste0("CA_",anonymized_patient_id),
                                                                  study == "brazil" ~ paste0("BR_",anonymized_patient_id),
                                                                  study == "germany_schulte" ~ paste0("GS_",anonymized_patient_id),
                                                                  study == "ita_valenti" ~ paste0("ITV_",anonymized_patient_id)
))

geno <- fread("../covid19-hgi-clinical-values/all_variant_rs35081325_rs11385942.tsv")

data <- merge(lab_data, geno[,-"study"], by="anonymized_patient_id", all.x=T)

##EUR
data_EUR <- data %>% filter(pop == "EUR")

##max
df <- data %>% filter(study != "ita_valenti") %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(c(any_of(blood_values), -lab_result_date)),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(c(any_of(blood_values[-1])))

df1 <- df %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

out <- data.frame(matrix(0,31,5))
colnames(out) <- c("biomarker", "OR", "OR_se", "pvalue","N")
for(i in c(2:4,6:12,14:28,30:33)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis + sex + study")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(13,29)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis + sex")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(5)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

out[order(out$pvalue),]
#max 
data_EUR <- data_EUR %>% filter(study != "ita_valenti")

df <- data_EUR %>% filter(study != "ita_valenti") %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(c(any_of(blood_values), -lab_result_date)),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(c(any_of(blood_values[-1])))

df1 <- df %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

out <- data.frame(matrix(0,31,5))
colnames(out) <- c("biomarker", "beta", "se", "pvalue","N")
for(i in c(2:4,6:12,14:28,30:33)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + sex + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(5)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(13,29)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

out[order(out$pvalue),]

#first_value
get_first <- function(x){
  tmp <- x[!is.na(x)]
  tmp[1]
}
df <- data %>% filter(study != "ita_renieri") %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(c(any_of(blood_values), -lab_result_date)),get_first) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(c(any_of(blood_values[-1])))

df1 <- df %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

out <- data.frame(matrix(0,31,5))
colnames(out) <- c("biomarker", "OR", "OR_se", "pvalue","N")
for(i in c(2:3,6:12,14:28,30:33)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis + sex + study")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(13)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis + sex")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(4,5)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("resp_mild ~ value + age_at_diagnosis")), dat=df2,family = "binomial")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- exp(summary(LM)$coefficient[2,1])
  out[i-1,3] <- exp(summary(LM)$coefficient[2,2])
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(29)){
  df2 <- df1 %>% filter(name == blood_values[i])
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- NA
  out[i-1,3] <- NA
  out[i-1,4] <- NA
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

out[order(out$pvalue),]
#max 
data_EUR <- data %>% filter(pop == "EUR")

df <- data_EUR %>% filter(study != "ita_renieri") %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(c(any_of(blood_values), -lab_result_date)),get_first) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(c(any_of(blood_values[-1])))

df1 <- df %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

out <- data.frame(matrix(0,31,5))
colnames(out) <- c("biomarker", "beta", "se", "pvalue","N")
for(i in c(2:3,6:12,14:28,30:33)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + sex + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(5)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(13)){
  df2 <- df1 %>% filter(name == blood_values[i])
  LM <- glm(as.formula(paste0("value ~ I(2 - rs35081325) + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- summary(LM)$coefficient[2,1]
  out[i-1,3] <- summary(LM)$coefficient[2,2]
  out[i-1,4] <- summary(LM)$coefficient[2,4]
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

for(i in c(4,29)){
  df2 <- df1 %>% filter(name == blood_values[i])
  out[i-1,1] <- blood_values[i]
  out[i-1,2] <- NA
  out[i-1,3] <- NA
  out[i-1,4] <- NA
  out[i-1,5] <- dim(df2[!is.na(df2$value),])[1]
}

out[order(out$pvalue),]
