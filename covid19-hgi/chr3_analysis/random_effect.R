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

#lab_data$resp_mild <- ifelse(lab_data$highest_respiratory_support >= 0, 1, 0)
lab_data$resp_severe <- ifelse(lab_data$highest_respiratory_support >= 1, 1, 0)

lab_data <- lab_data %>% mutate(anonymized_patient_id = case_when(study == "sweden" ~ paste0("SW_",anonymized_patient_id),
                                                                  study == "belgium" ~ paste0("BB_",anonymized_patient_id),
                                                                  study == "ita_renieri" ~ paste0("ITR_",anonymized_patient_id),
                                                                  study == "spain_bujanda" ~ paste0("SBuj_",anonymized_patient_id),
                                                                  study == "spain_butti" ~ paste0("SBut_",anonymized_patient_id),
                                                                  study == "canada" ~ paste0("CA_",anonymized_patient_id),
                                                                  study == "brazil" ~ paste0("BR_",anonymized_patient_id),
                                                                  study == "germany_schulte" ~ paste0("GS_",anonymized_patient_id),
                                                                  study == "ita_valenti" ~ paste0("ITV_",anonymized_patient_id)
))

geno <- readRDS("../covid19-hgi-clinical-values/bgen/all_ans_pca_sex_raw.rds")

data <- merge(lab_data, geno[,-"study"], by="anonymized_patient_id", all.x=T)

##longitudinal
df <- data %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
 filter(pop == "EUR") %>% group_by(anonymized_patient_id) %>%
  pivot_longer(c(any_of(blood_values[-1])))
library(lme4)
blood_values_rev <- blood_values[c(2:3,6:11,13:17,19:20,22:23,26:28,30:33)]

for(j in seq(1,10)){
  out <- data.frame(matrix(0,24,6))
  colnames(out) <- c("biomarker", "beta", "se","pvalue","N","SNP")
  for(i in c(1:24)){
    df2 <- df %>% filter(name == blood_values_rev[i]) %>% filter(!is.na(value))
    df2 <- df2 %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
      mutate(value = scale(value)) %>% ungroup()
    if(length(unique(df2$study)) > 1 & length(unique(df2$sex)) > 1){
      LM <- glmer(as.formula(paste0("value ~ `",snps[j],"` + lab_result_date + age_at_diagnosis + (1 + lab_result_date | anonymized_patient_id) + sex + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=df2,family = "gaussian")
      out[i,1] <- blood_values_rev[i]
      out[i,2] <- summary(LM)$coefficients[2,1]
      out[i,3] <- summary(LM)$coefficients[2,2]
      out[i,4] <- sjstats::p_value(LM)[2,2]
      out[i,5] <- dim(df2[!is.na(df2$value),])[1]
      out[i,6] <- snps[j]
    }
    if(length(unique(df2$study)) == 1){
      LM <- glmer(as.formula(paste0("value ~ `",snps[j],"` + lab_result_date + age_at_diagnosis + (1 + lab_result_date | anonymized_patient_id) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=df2,family = "gaussian")
      out[i,1] <- blood_values_rev[i]
      out[i,2] <- summary(LM)$coefficients[2,1]
      out[i,3] <- summary(LM)$coefficients[2,2]
      out[i,4] <- sjstats::p_value(LM)[2,2]
      out[i,5] <- dim(df2[!is.na(df2$value),])[1]
      out[i,6] <- snps[j]
    }
  }
  #print(out[order(out$pvalue),])
  write.table(out, file="longitudinal_EUR_biomarker.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}

blood_values_rev <- blood_values[c(2:3,6:11,13:17,19:20,22:23,26:28,30:33)]

for(j in seq(1,10)){
  out <- data.frame(matrix(0,24,6))
  colnames(out) <- c("biomarker", "beta", "se","pvalue","N","SNP")
  for(i in c(1:24)){
    df2 <- df %>% filter(name == blood_values_rev[i]) %>% filter(!is.na(value))
    df2 <- df2 %>% group_by(name) %>% mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
      mutate(value = scale(value)) %>% ungroup()
    if(length(unique(df2$study)) > 1 & length(unique(df2$sex)) > 1){
      LM <- glmer(as.formula(paste0("value ~ `",snps[j],"` + lab_result_date + `",snps[j],"`*lab_result_date + age_at_diagnosis + (1 + lab_result_date | anonymized_patient_id) + sex + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=df2,family = "gaussian")
      out[i,1] <- blood_values_rev[i]
      out[i,2] <- summary(LM)$coefficients[dim(summary(LM)$coefficients)[1],1]
      out[i,3] <- summary(LM)$coefficients[dim(summary(LM)$coefficients)[1],2]
      out[i,4] <- sjstats::p_value(LM)[dim(summary(LM)$coefficients)[1],2]
      out[i,5] <- dim(df2[!is.na(df2$value),])[1]
      out[i,6] <- snps[j]
    }
    if(length(unique(df2$study)) == 1){
      LM <- glmer(as.formula(paste0("value ~ `",snps[j],"` + lab_result_date + `",snps[j],"`*lab_result_date + age_at_diagnosis + (1 + lab_result_date | anonymized_patient_id) + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=df2,family = "gaussian")
      out[i,1] <- blood_values_rev[i]
      out[i,2] <- summary(LM)$coefficients[dim(summary(LM)$coefficients)[1],1]
      out[i,3] <- summary(LM)$coefficients[dim(summary(LM)$coefficients)[1],2]
      out[i,4] <- sjstats::p_value(LM)[dim(summary(LM)$coefficients)[1],2]
      out[i,5] <- dim(df2[!is.na(df2$value),])[1]
      out[i,6] <- snps[j]
    }
  }
  #print(out[order(out$pvalue),])
  write.table(out, file="longitudinal_slope_EUR_biomarker.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}
