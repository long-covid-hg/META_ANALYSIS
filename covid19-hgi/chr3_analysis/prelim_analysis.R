library(lme4)

a <- bel_manifest %>% select(any_of(c(identification, "age_at_diagnosis", "sex","highest_respiratory_support"))) %>% mutate(study ="belgium")
b <- swe_manifest %>% select(any_of(c(identification, "age_at_diagnosis", "sex","highest_respiratory_support"))) %>% mutate(study ="sweden")
c <- spain_bujanda_manifest %>% select(any_of(c(identification, "age_at_diagnosis", "sex","highest_respiratory_support"))) %>% mutate(study ="spain_bujanda")
d <- spain_butti_manifest %>% select(any_of(c(identification, "age_at_diagnosis", "sex","highest_respiratory_support"))) %>% mutate(study ="spain_butti")

b$anonymized_patient_id <- as.character(b$anonymized_patient_id)
c$anonymized_patient_id <- as.character(c$anonymized_patient_id)
d$anonymized_patient_id <- as.character(d$anonymized_patient_id)
d$age_at_diagnosis <- as.numeric(d$age_at_diagnosis)
d$sex <- as.numeric(d$sex)

a <- a %>% mutate(highest_respiratory_support = case_when(highest_respiratory_support == 2 | highest_respiratory_support == 1 ~ 1,
                                                          TRUE ~ 0))
b <- b %>% mutate(highest_respiratory_support = case_when(highest_respiratory_support == 2 | highest_respiratory_support == 1 ~ 1,
                                                          TRUE ~ 0))
c <- c %>% mutate(highest_respiratory_support = case_when(highest_respiratory_support == 2 | highest_respiratory_support == 1 ~ 1,
                                                          TRUE ~ 0))
d <- d %>% mutate(highest_respiratory_support = case_when(highest_respiratory_support == 2 | highest_respiratory_support == 1 ~ 1,
                                                          TRUE ~ 0))
outcome <- bind_rows(list(a,b,c,d))
final <- merge(lab_data, outcome, by=c("anonymized_patient_id", "study"), all.x=T)
final <- final %>% filter(lab_result_date <= 30 & lab_result_date >= 0)

LM <- glmer(as.formula(paste0("highest_respiratory_support ~ lab_wbc + lab_wbc*lab_result_date + (1 + lab_result_date | anonymized_patient_id) + age_at_diagnosis + sex + study")), dat=final,family = "binomial")
summary(LM)

LM <- glm(as.formula(paste0("highest_respiratory_support ~ lab_wbc + lab_wbc*lab_result_date + age_at_diagnosis + sex + study")), dat=final,family = "binomial")
summary(LM)

LM <- glmer(as.formula(paste0("highest_respiratory_support ~ lab_il_6 + lab_il_6*lab_result_date + (1 + lab_result_date | anonymized_patient_id) + age_at_diagnosis + sex + study")), dat=final,family = "binomial")
summary(LM)



df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= 0) %>% 
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(-anonymized_patient_id,-study),max,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>% 
  pivot_longer(!c(anonymized_patient_id,study))

final <- merge(df, outcome, by=c("anonymized_patient_id", "study"), all.x=T)

final1 <- final[final$name == "lab_il_6" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ value + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)

final1 <- final[final$name == "lab_procalcitonin" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ value + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)

final1 <- final[final$name == "lab_serum_ferritin" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ value + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)

final1 <- final[final$name == "lab_wbc" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ value + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)

df <- lab_data %>% filter(lab_result_date <= 30 & lab_result_date >= 0) %>% 
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(-anonymized_patient_id,-study),min,na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>% 
  pivot_longer(!c(anonymized_patient_id,study))
final <- merge(df, outcome, by=c("anonymized_patient_id", "study"), all.x=T)

final1 <- final[final$name == "lab_wbc" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ value + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)


#variant
genes <- fread(paste0("../covid19-hgi-clinical-values/hgi_belgium/bel_variant_rs17078348.tsv"))
genes$study <- "belgium"
genes1 <- fread(paste0("../covid19-hgi-clinical-values/hgi_spain_bujanda/spain_bujanda_variant_rs17078348.tsv"))
genes1$study <- "spain_bujanda"
genes2 <- fread(paste0("../covid19-hgi-clinical-values/hgi_swe/swe_variant_rs17078348.tsv"))
genes2$study <- "sweden"
gene <- rbind(genes, genes1, genes2)

final <- merge(final, gene, by=c("anonymized_patient_id", "study"), all.x=T)
final <- final[!is.na(final$rs17078348),]

final1 <- final[final$name == "lab_il_6" & !is.infinite(final$value),]
final1$value <- log10(final1$value)
LM <- glm(as.formula(paste0("value ~ I(2-rs17078348) + age_at_diagnosis + sex + study")), dat=final1,family = "gaussian")
summary(LM)

final1 <- final[final$name == "lab_procalcitonin" & !is.infinite(final$value),]
final1$value <- log10(final1$value)
hist(final1$value)
LM <- glm(as.formula(paste0("value ~ I(2-rs17078348) + age_at_diagnosis + sex + study")), dat=final1,family = "gaussian")
summary(LM)

final1 <- final[final$name == "lab_wbc" & !is.infinite(final$value),]
final1$value <- log10(final1$value)
hist(final1$value)
LM <- glm(as.formula(paste0("value ~ I(2-rs17078348) + age_at_diagnosis + sex + study")), dat=final1,family = "gaussian")
summary(LM)

final1 <- final[final$name == "lab_serum_ferritin" & !is.infinite(final$value),]
final1$value <- log10(final1$value)
hist(final1$value)
LM <- glm(as.formula(paste0("value ~ I(2-rs17078348) + age_at_diagnosis + sex + study")), dat=final1,family = "gaussian")
summary(LM)

final1 <- final[final$name == "lab_wbc" & !is.infinite(final$value),]
LM <- glm(as.formula(paste0("highest_respiratory_support ~ I(2-rs17078348) + age_at_diagnosis + sex + study")), dat=final1,family = "binomial")
summary(LM)
