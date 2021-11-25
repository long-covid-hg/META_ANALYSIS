setwd("/project/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/06.clinical_values")
path <- "/project/richards/tomoko.nakanishi/09.COVID19/scratch/01.UKBB/06.clinical_values/"
data <- readRDS("../../../scratch/data/01.UKBB/clinical_value_20201024.rds")
geno <- readRDS(paste0(path,"allå_chr.rds"))

library(tidyr)
library(dplyr)
final <- data %>% inner_join(geno, by=c("anonymized_patient_id" = "FID"))

eur <- fread("../../../scratch/data/01.UKBB/EMC-White-British.20pc.txt")
colnames(eur) <- c("V1", paste0("PC",1:20))
final <- final %>% inner_join(eur, by=c("anonymized_patient_id" = "V1"))

final <- final %>% mutate(a1 = ifelse(death == 1 | resp_severe == 1, 1, 0),
                          resp_severe = ifelse(death == 1 & resp_severe == 0, -1, resp_severe),
                          resp_mild = ifelse(death == 1 & resp_mild == 0, -1, resp_mild),
                          hosp_vte = case_when(hosp_thrombo == 1 | hosp_dvt == 1 ~ 1,
                                               death == 1 & hosp_thrombo == 0 & hosp_dvt == 0 ~ -1,
                                               TRUE ~ 0),
                          hosp_vte = case_when(hosp_thrombo == 1 | hosp_dvt == 1 ~ 1,
                                               death == 1 & hosp_thrombo == 0 & hosp_dvt == 0 ~ -1,
                                               TRUE ~ 0),
                          hosp_cvd = case_when(hosp_stroke == 1 | hosp_infarction == 1 ~ 1,
                                               death == 1 & hosp_stroke == 0 & hosp_infarction == 0 ~ -1,
                                               TRUE ~ 0),
                          hosp_renal = ifelse(death == 1 & hosp_renal == 0, -1, hosp_renal),
                          hosp_bleeding = ifelse(death == 1 & hosp_bleeding == 0, -1, hosp_bleeding)
                          )

final <- final %>% mutate_at(.vars = vars(pheno), .funs = funs(ifelse(. < 0,NA,.)))

snps <- colnames(final)[grepl("rs", colnames(final)) & !grepl("HET", colnames(final))]
pheno <- c("death","a1","resp_severe", "resp_mild" , colnames(final)[grepl("hosp", colnames(final))])
pheno <- pheno[-5]
out <- data.table(matrix(0,10,7))
colnames(out) <- c("SNP", "outcome", "beta", "se", "pval", "Ncase", "Ncontrol")
for(i in seq(1,12)){
  out$outcome <- rep(pheno[i], 10)
  for(j in seq(1,10)){
    out$SNP[j] <- snps[j]
    LM <- glm(as.formula(paste0(pheno[i],"~ `",snps[j],"` + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=final,family = "binomial")
    out$beta[j] <- summary(LM)$coefficients[2,1]
    out$se[j] <- summary(LM)$coefficients[2,2]
    out$pval[j] <- summary(LM)$coefficients[2,4]
    out$Ncase[j] <- sum(final[,pheno[i]] == 1,na.rm=T)
    out$Ncontrol[j] <- sum(!is.na(final[,pheno[i]])) - out$Ncase[j]
  }
  write.table(out, file="/project/richards/tomoko.nakanishi/09.COVID19/results/01.UKBB/06.clinical_values/UKB.outcome.result.tsv", quote=F, col.names = F, row.names = F, append=T)
}

out <- fread("/project/richards/tomoko.nakanishi/09.COVID19/results/01.UKBB/06.clinical_values/UKB.outcome.result.tsv")
colnames(out) <- c("SNP", "outcome", "beta", "se", "pval", "Ncase", "Ncontrol")

out1 <- out %>% filter(pval < 0.05) %>% filter(pval != 0)
out1[order(out1$pval),]

out1 <- out %>% filter(SNP == "rs13050728_C")
out1

cor(final$resp_severe, final$hosp_renal, use = "complete.obs")
cor(final$resp_severe, final$hosp_vte, use = "complete.obs")

LM <- glm(as.formula(paste0(pheno[9],"~ `",snps[1],"` + ",pheno[3]," + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=final,family = "binomial")
LM <- glm(as.formula(paste0(pheno[11],"~ `",snps[1],"` + ",pheno[3]," + age_at_diagnosis + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")), dat=final,family = "binomial")
