setwd("/home/tomoko/02.lab_data")

lab_data <- fread("All.cohort.lab.tsv")
#lab_data <- bind_rows(list(a,b,c,d,e,f,g))

demographics <- c("age_at_diagnosis","sex")
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

lab_data <- lab_data %>% 
  mutate_at(.vars = vars(-c(anonymized_patient_id)), .funs = funs(ifelse(. < 0,NA,.)))

lab_data <- lab_data %>% mutate(lab_n_l_ratio = case_when(is.na(lab_n_l_ratio) ~ lab_neutrophils/lab_lymphocytes,
                                                          TRUE ~ lab_n_l_ratio))

lab_data$resp <- ifelse(lab_data$highest_respiratory_support >= 1, 1, 0)

lab_data <- lab_data %>% mutate(anonymized_patient_id = case_when(study == "sweden" ~ paste0("swe_",anonymized_patient_id),
                                                      study == "belgium" ~ paste0("BB_",anonymized_patient_id),
                                                      study == "italy" ~ paste0("IT_",anonymized_patient_id),
                                                      study == "spain_bujanda" ~ paste0("SB_",anonymized_patient_id),
                                                      study == "spain_butti" ~ paste0("SB_",anonymized_patient_id),
                                                      study == "canada" ~ paste0("Ca_",anonymized_patient_id),
                                                      study == "brazil" ~ paste0("Bra_",anonymized_patient_id)))


#max value distribution
library(ggridges)
png("lab_distribution_raw.png", width=3000, height=2000)
df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,blood_values[2])) %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(-anonymized_patient_id,-study),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(blood_values[2])
df <- df %>% mutate(value = ifelse(!is.na(value) & !is.infinite(value), value, NA))

a2 <- ggplot(df,aes(x=value,y=study,fill = study)) + geom_density_ridges(scale = 5) + ggtitle(blood_values[2])

for(i in seq(3,30)){
  df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,blood_values[i])) %>%
    group_by(anonymized_patient_id) %>%
    mutate_at(vars(-anonymized_patient_id,-study),max,na.rm=T) %>% 
    ungroup()  %>% 
    distinct(anonymized_patient_id, .keep_all=TRUE) %>%
    pivot_longer(blood_values[i])
  df <- df %>% mutate(value = ifelse(!is.na(value) & !is.infinite(value), value, NA))
  a2 <- a2 + ggplot(df,aes(x=value,y=study,fill = study)) + geom_density_ridges(scale = 5) + ggtitle(blood_values[i])
}
a2
dev.off()

i <- 2
df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,blood_values[i])) %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(-anonymized_patient_id,-study),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>%
  pivot_longer(blood_values[i])
df <- df[!is.na(df$value) & !is.infinite(df$value),]
df$value <- log10(df$value)
df$value[is.infinite(df$value)] <- min(df$value[!is.infinite(df$value)])/2
df$value <- (df$value - mean(df$value, na.rm = T))/sd(df$value, na.rm = T)
a2 <- ggplot(df,aes(x=value,y=study,fill = study)) + geom_density_ridges(scale = 5) + ggtitle(blood_values[i])

for(i in seq(3,30)){
  df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= -2) %>% select(c(identification,demographics,study,blood_values[i])) %>%
    group_by(anonymized_patient_id) %>%
    mutate_at(vars(-anonymized_patient_id,-study),max,na.rm=T) %>% 
    ungroup()  %>% 
    distinct(anonymized_patient_id, .keep_all=TRUE) %>%
    pivot_longer(blood_values[i])
  df <- df[!is.na(df$value) & !is.infinite(df$value),]
  df$value <- log10(df$value)
  df$value[is.infinite(df$value)] <- min(df$value[!is.infinite(df$value)])/2
  df$value <- (df$value - mean(df$value, na.rm = T))/sd(df$value, na.rm = T)
  a2 <- a2 + ggplot(df,aes(x=value,y=study,fill = study)) + geom_density_ridges(scale = 5) + ggtitle(blood_values[i])
}
png("lab_distribution_log.png", width=3000, height=2000)
a2
dev.off()
