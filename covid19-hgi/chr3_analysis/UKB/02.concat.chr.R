setwd("/project/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/06.clinical_values")
library(data.table)
library(tidyr)
library(dplyr)
path <- "/project/richards/tomoko.nakanishi/09.COVID19/scratch/01.UKBB/06.clinical_values/"
chr3 <- fread(paste0(path,"3.raw"))
chr6 <- fread(paste0(path,"6.raw")) %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
chr9 <- fread(paste0(path,"9.raw")) %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
chr12 <- fread(paste0(path,"12.raw")) %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
chr19 <- fread(paste0(path,"19.raw")) %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
chr21 <- fread(paste0(path,"21.raw")) %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))

data <- chr3 %>% merge(chr6, by="FID") %>% merge(chr9, by="FID") %>% merge(chr12, by="FID") %>%
  merge(chr19, by="FID") %>% merge(chr21, by="FID")

saveRDS(data, file=paste0(path,"all_chr.rds"))
