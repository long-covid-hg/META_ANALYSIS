rm(list=ls())

library(readxl)   
library(data.table)
library(dplyr)
library(stringr)

setwd('~/release_15042021/')


main_pheno <- fread('phenotypes_main_2021_05_12_wo_mismatch.tsv')
pcs <- fread("cov_15042021.qc.merged_controls_scores.tsv")
fam <- fread("cov_15042021.qc.merged_controls.fam")

# join with PCs and keep only samples in the fam/PCs file

pheno <- main_pheno %>% 
  inner_join(pcs, by = c("ID" = "s")) %>% 
  # change ID to match the ID in the imputed data
  mutate(ID = paste0("0_", ID)) %>% 
  select(ID, Cohort, PC1:PC20, Age, Sex, Age2, AgeSex, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2) %>%
  distinct()

# check some duplicated IDs
dups <- pheno$ID[duplicated(pheno$ID)]
# add row index to pheno file to understand when the dup was added
pheno$row_id <- 1:nrow(pheno)
pheno_dups <- pheno %>% 
  filter(ID %in% dups) %>%
  arrange(ID) %>%
  select(ID, Cohort, Age, Sex, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2, row_id)

# 10 IDs are duplicated, cohort italy, duplicated with different ages and phenotype definitions. 
# Keep rows where sex is consistent, only as C2 cases (we don't have age as covariate when running C2 for italy) ?
# Seem to be added at very different stages, we don't actually know if the samples is the same or the ID has been duplicated. Discard
fwrite(pheno_dups, "ita_dup.tsv", quote = F, na = "NA", sep = "\t")

pheno <- pheno %>%
  filter(!ID %in% dups) %>%
  select(-row_id)

# Add stratified analyses columns:
# - Analysis_X_M: males
# - Analysis_X_F: females
# - Analysis_X_LT60: <=60
# - Analysis_X_GT60: >60
head(pheno)

pheno <- pheno %>% 
  mutate(Analysis_A2_M = ifelse(Sex == 1, Analysis_A2, NA),
         Analysis_A2_F = ifelse(Sex == 2, Analysis_A2, NA),
         Analysis_A2_LT60 = ifelse(Age <= 60, Analysis_A2, NA),
         Analysis_A2_GT60 = ifelse(Age > 60, Analysis_A2, NA),
         Analysis_B1_M = ifelse(Sex == 1, Analysis_B1, NA),
         Analysis_B1_F = ifelse(Sex == 2, Analysis_B1, NA),
         Analysis_B1_LT60 = ifelse(Age <= 60, Analysis_B1, NA),
         Analysis_B1_GT60 = ifelse(Age > 60, Analysis_B1, NA),
         Analysis_B2_M = ifelse(Sex == 1, Analysis_B2, NA),
         Analysis_B2_F = ifelse(Sex == 2, Analysis_B2, NA),
         Analysis_B2_LT60 = ifelse(Age <= 60, Analysis_B2, NA),
         Analysis_B2_GT60 = ifelse(Age > 60, Analysis_B2, NA),
         Analysis_C2_M = ifelse(Sex == 1, Analysis_C2, NA),
         Analysis_C2_F = ifelse(Sex == 2, Analysis_C2, NA),
         Analysis_C2_LT60 = ifelse(Age <= 60, Analysis_C2, NA),
         Analysis_C2_GT60 = ifelse(Age > 60, Analysis_C2, NA))


# double check we match all IDs in the genotype file
system("gsutil cp gs://dsge-covid19-data/15042021/conf/bgen_samples.txt .")
bgen_samples <- fread("bgen_samples.txt", header = F)

no_match <- setdiff(bgen_samples$V1, pheno$ID)
length(intersect(bgen_samples$V1, pheno$ID))

pheno_analyzed <- pheno %>% 
  filter(ID %in% bgen_samples$V1)

# separate cohorts and correct stratified analyses (e.g. for cohorts with no age info for pop controls)
# write ID males samples for each cohort
belgium <- pheno %>%
  filter(Cohort == "belgium") %>% 
  mutate(Analysis_A2_LT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_LT60),
         Analysis_A2_GT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_GT60),
         Analysis_B2_LT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_LT60),
         Analysis_B2_GT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_GT60),
         Analysis_C2_LT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_LT60),
         Analysis_C2_GT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_GT60))

# check sex mismatch
sex_mism <- fread('sex_mism_belgium', header = FALSE)

sex_mism$V1 <- paste0("0_", sex_mism$V1)

sex_mism_belgium <- belgium %>% 
  inner_join(sex_mism, by = c("ID" = "V1")) %>% 
  filter(Sex != V2)


write.table(belgium, file = "pheno_gwas/pheno_cov_belgium_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

belgium %>% 
  filter(Sex == 1) %>% 
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_belgium_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)

egypt <- pheno %>%
  filter(Cohort == "egypt")
summary(egypt$Age) # No NA's

write.table(egypt, file = "pheno_gwas/pheno_cov_egypt_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

egypt %>% 
  filter(Sex == 1) %>% 
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_egypt_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)


germany <- pheno %>%
  filter(Cohort == "germany")
summary(germany$Age) # No NA's

write.table(germany, file = "pheno_gwas/pheno_cov_germany_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

germany %>% 
  filter(Sex == 1) %>% 
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_germany_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)


# # # IRAN
iran <- pheno %>%
  filter(Cohort == "iran")
summary(iran$Age) # No NA's

write.table(iran, file = "pheno_gwas/pheno_cov_iran_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

iran %>% 
  filter(Sex == 1) %>%
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_iran_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)


italy <- pheno %>%
  filter(Cohort == "italy") %>% 
  mutate(Analysis_A2_LT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_LT60),
         Analysis_A2_GT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_GT60),
         Analysis_B2_LT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_LT60),
         Analysis_B2_GT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_GT60),
         Analysis_C2_LT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_LT60),
         Analysis_C2_GT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_GT60))

write.table(italy, file = "pheno_gwas/pheno_cov_italy_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

italy %>% 
  filter(Sex == 1) %>% 
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_italy_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)


sweden <- pheno %>%
  filter(Cohort == "sweden")
summary(sweden$Age) # 3748 NA's 
summary(sweden$Age[sweden$Analysis_A2 == 0]) # all controls have NA age

sweden <- pheno %>%
  filter(Cohort == "sweden") %>% 
  mutate(Analysis_A2_LT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_LT60),
         Analysis_A2_GT60 = ifelse(is.na(Age) & Analysis_A2 == 0, 0, Analysis_A2_GT60),
         Analysis_B2_LT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_LT60),
         Analysis_B2_GT60 = ifelse(is.na(Age) & Analysis_B2 == 0, 0, Analysis_B2_GT60),
         Analysis_C2_LT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_LT60),
         Analysis_C2_GT60 = ifelse(is.na(Age) & Analysis_C2 == 0, 0, Analysis_C2_GT60))

write.table(sweden, file = "pheno_gwas/pheno_cov_sweden_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

sweden %>% 
  filter(Sex == 1) %>% 
  select(ID) %>% 
  write.table(file = "pheno_gwas/males_sweden_15042021.tsv", append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = FALSE, quote = FALSE)


# Check analysis to run for each cohort
ana_cols <- colnames(pheno)[grepl("^Analysis_", colnames(pheno))]

# belgium
res <- NULL
for (ana in ana_cols) {
  t <- table(belgium[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}

res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis

# egypt
res <- NULL
for (ana in ana_cols) {
  t <- table(egypt[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}

res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis

# germany
res <- NULL
for (ana in ana_cols) {
  t <- table(germany[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}

res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis


# iran
res <- NULL
for (ana in ana_cols) {
  t <- table(iran[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}

res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis


# italy
res <- NULL
for (ana in ana_cols) {
  t <- table(italy[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}
res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis


# sweden
res <- NULL
for (ana in ana_cols) {
  t <- table(sweden[[ana]])
  r <- c(ana, t["0"], t["1"])
  res <- rbind(res, r)
}
res <- data.frame(res, stringsAsFactors = F)
colnames(res) <- c("Analysis", "N.controls", "N.cases")
res <- res %>% 
  filter(as.numeric(N.cases) >= 100)
res
res$Analysis


# # # 
control_to_keep <- pheno %>% 
  filter(Analysis_C2 == 0) %>% 
  pull(ID)

fwrite(list(control_to_keep), "controls.to.keep.txt")
