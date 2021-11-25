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

##longitudinal
df <- data_EUR %>% filter(study != "ita_renieri" & study != "ita_valenti") %>% filter(lab_result_date <= 30 & lab_result_date >= -2)  %>%
  group_by(anonymized_patient_id) %>%
  pivot_longer(c(any_of(blood_values[-1])))
library(lme4)
i <- 17
df2 <- df %>% filter(name == blood_values[i]) %>% filter(!is.na(value)) %>% 
  mutate(value = log10(value)) %>% mutate(value = ifelse(is.infinite(value), min(value[!is.infinite(value) & !is.na(value)])/2, value)) %>%
  mutate(value = (value - mean(value, na.rm = T))/sd(value, na.rm = T))

LM <- glmer(as.formula(paste0("value ~ I(2 - rs35081325) + lab_result_date + age_at_diagnosis + (1 + lab_result_date | anonymized_patient_id) + sex + study + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), dat=df2,family = "gaussian")
summary(LM)
sjstats::p_value(LM)

sum(!is.na(df2$value))/length(unique(df2$anonymized_patient_id))
