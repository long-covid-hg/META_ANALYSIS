library(readxl)   
library(data.table)
library(dplyr)
library(stringr)

main_pheno <- fread('phenotypes_main_2021_04_06.tsv')
# Update IDs for the first Tanta Egyptian batch we received. Add "COV." for both cases and controls
main_pheno <- main_pheno %>% mutate(ID = ifelse(Cohort == "egypt", paste0("COV.", ID), ID))

########### SWEDEN ###########
sweden <- as.data.frame(read_excel('age_and_sex_karolinska.xlsx'))
sweden_edited <- sweden %>% dplyr::mutate(Cohort="sweden",
                                          ID=anonymized_patient_id,
                                          Age=as.numeric(age_at_diagnosis),
                                          Sex=as.numeric(sex) + 1,
                                          Age2=Age**2,
                                          AgeSex=Age*Sex,
                                          Analysis_A1=1,
                                          Analysis_A2=1,
                                          Analysis_B1=1,
                                          Analysis_B2=1,
                                          Analysis_C2=1) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)

########### Egypt Mansoura Controls ###########
egypt_mansoura_controls <- as.data.frame(read_excel('control_mansoura_batch_1.xlsx', sheet = "Sheet1"))
nrow(egypt_mansoura_controls)

# Two samples, MANC45 and MANC96, were excluded, so we have to remove them
egypt_man_con_edited <- egypt_mansoura_controls %>% filter(!(age == 'excluded'))
egypt_man_con_edited <- egypt_man_con_edited %>% dplyr::mutate(Cohort="egypt",
                                                               ID=...7,
                                                               Age=as.numeric(age),
                                                               Sex=as.numeric(gender),
                                                               Age2=Age**2,
                                                               AgeSex=Age*Sex,
                                                               Analysis_A1=NA,
                                                               Analysis_A2=0,
                                                               Analysis_B1=NA,
                                                               Analysis_B2=0,
                                                               Analysis_C2=0) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit ID so it matches with FAM
str_sub(egypt_man_con_edited$ID, 4, 3) <- "_"
egypt_man_con_edited <- egypt_man_con_edited %>% dplyr::mutate(ID=paste0("COV.", str_to_upper(ID)))


########### Egypt Mansoura Cases ###########
egypt_mansoura_patients <- as.data.frame(read_excel('patients_mansoura_batch_1.xlsx', sheet = 'FREEZE_2_data_dictionary'))
egypt_mansoura_patient_transposed <- transpose(egypt_mansoura_patients)
# V4 = ID, V5 = age_at_diagnosis, V6 = Sex, V52 = hospitilization, V64 = highest_respiratory_support, V66 = icu_admit
egypt_mansoura_patient_transposed <- egypt_mansoura_patient_transposed %>% dplyr::select(V4, V5, V6, V52, V64, V66)
names(egypt_mansoura_patient_transposed) <- c("ID", "Age", "Sex", "hospitilization", "highest_respiratory_support", "icu_admit")
# remove lines 1-4 since they will have have previous column names
egypt_mansoura_patient_transposed <- egypt_mansoura_patient_transposed[-c(1, 2, 3, 4), ] 

egypt_man_cas_edited <- egypt_mansoura_patient_transposed %>% dplyr::mutate(Cohort="egypt",
                                                                     ID=ID,
                                                                     Age=as.numeric(Age),
                                                                     Sex=as.numeric(Sex) + 1,
                                                                     Age2=Age**2,
                                                                     AgeSex=Age*Sex,
                                                                     Analysis_A1=NA,
                                                                     Analysis_A2=ifelse(icu_admit == 1 | highest_respiratory_support >=0, 1, NA),
                                                                     Analysis_B1=ifelse(hospitilization == 1, 1, 0),
                                                                     Analysis_B2=ifelse(hospitilization == 1, 1, NA),
                                                                     Analysis_C2=1) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit ID so it matches with FAM
str_sub(egypt_man_cas_edited$ID, 4, 3) <- "_"
egypt_man_cas_edited <- egypt_man_cas_edited %>% dplyr::mutate(ID=paste0("COV.", str_to_upper(ID)))


########### Germany cases ###########
german_cas <- as.data.frame(read_excel('FIMM_MetaData_COVID-19_GWAS_part2_20210329.xlsx'))
# ! REMOVE CoMRI-DNA-141, CoMRI-DNA-165, CoMRI-DNA-177, and CoMRI-DNA-197 from genotype data (ger, withdraw consent)
german_cas <- german_cas %>%
  dplyr::filter(ID != c('CoMRI-DNA-141','CoMRI-DNA-165','CoMRI-DNA-177', 'CoMRI-DNA-197')) 


german_cas_edited <- german_cas %>% dplyr::mutate(Cohort="germany",
                                                  ID=ID,
                                                  Age=Age,
                                                  Sex=ifelse(Sex == "1=male", 1, 2),
                                                  Age2=Age**2,
                                                  AgeSex=Age*Sex,
                                                  Analysis_A1=NA,
                                                  Analysis_A2=ifelse(Phenotype == "A2", 1, NA),
                                                  Analysis_B1=ifelse(Phenotype == "A2" | Phenotype == "B1" | Phenotype == "B2", 1, NA),
                                                  Analysis_B2=ifelse(Phenotype == "A2" | Phenotype == "B1" | Phenotype == "B2", 1, NA),
                                                  Analysis_C2=ifelse(Phenotype == "A2" | Phenotype == "B1" | Phenotype == "B2" | Phenotype == "C2", 1, NA)) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit the IDs
german_cas_edited <- german_cas_edited %>% dplyr::mutate(ID=paste0("COV.", str_to_upper(str_replace_all(ID, "-", "_"))))
head(german_cas_edited)


########### Italian cases ###########
italy_cas_one <- as.data.frame(read_excel('samples-sheet_clinical_cathegory_19-01_2021_shipment.xlsx'))
names(italy_cas_one)[1:5] <- c("ID", "family_relationship", "gender", "age", "clinical_category")
italy_cas_one <- italy_cas_one %>%
  mutate(age = replace(age,age=="NA",0))

italy_cas_one_edited <- italy_cas_one %>% dplyr::mutate(Cohort="italy",
                                                        ID=ID,
                                                        Age=ifelse(!is.na(age), as.numeric(age), NA),
                                                        Sex=gender,
                                                        Age2=ifelse(!is.na(Age), Age**2, NA),
                                                        AgeSex=Age*Sex,
                                                        Analysis_A1=NA,
                                                        Analysis_A2=ifelse(clinical_category >= 2, 1, NA),
                                                        Analysis_B1=ifelse(clinical_category > 0, 1, 0),
                                                        Analysis_B2=ifelse(clinical_category >= 1, 1, NA),
                                                        Analysis_C2=1) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit the IDs
italy_cas_one_edited <- italy_cas_one_edited %>% dplyr::mutate(ID=paste0("COV.", str_replace_all(ID, "-", "_")))


italy_cas_two <- as.data.frame(read_excel('Missing_data_samples_March_2021.xlsx'))
names(italy_cas_two)[1:5] <- c("ID", "family_relationship", "gender", "age", "clinical_category")
italy_cas_two <- italy_cas_two %>%
  mutate(age = replace(age,age=="NA",0))

italy_cas_two <- italy_cas_two %>%
  dplyr::filter(!(ID %in% c('COV2735-1159','COV3408-1416','COV725-372', 'COV1616-683', 'COV1114-529', 'COV1857-789'))) 

italy_cas_two_edited <- italy_cas_two %>% dplyr::mutate(Cohort="italy",
                                                        ID=ID,
                                                        Age=ifelse(!is.na(age), as.numeric(age), NA),
                                                        Sex=as.numeric(gender),
                                                        Age2=ifelse(!is.na(Age), Age**2, NA),
                                                        AgeSex=Age*Sex,
                                                        Analysis_A1=NA,
                                                        Analysis_A2=ifelse(clinical_category >= 2, 1, NA),
                                                        Analysis_B1=ifelse(clinical_category > 0, 1, 0),
                                                        Analysis_B2=ifelse(clinical_category >= 1, 1, NA),
                                                        Analysis_C2=1) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit the IDs
italy_cas_two_edited <- italy_cas_two_edited %>% dplyr::mutate(ID=paste0("COV.", str_replace_all(ID, "-", "_")))


iran <- as.data.frame(read_excel('DDRI-FIMM-data.xlsx', sheet = 'DDRI-FIMM'))
head(iran$hospitalization)
iran_edited <- iran %>% dplyr::mutate(Cohort="iran",
                                      ID=SAMPLE_IDENTIFIER,
                                      Age=as.numeric(age_at_diagnosis...5),
                                      Sex=as.numeric(sex...6) + 1,
                                      Age2=Age**2,
                                      AgeSex=Age*Sex,
                                      Analysis_A1=NA,
                                      Analysis_A2=ifelse(icu_admit == 1 | highest_respiratory_support >=0, 1, NA),
                                      Analysis_B1=ifelse(hospitalization == 1, 1, 0),
                                      Analysis_B2=ifelse(hospitalization == 1, 1, NA),
                                      Analysis_C2=1) %>%
  select(Cohort, ID, Age, Sex, Age2, AgeSex, Analysis_A1, Analysis_A2, Analysis_B1, Analysis_B2, Analysis_C2)
# Edit the IDs
iran_edited <- iran_edited %>% dplyr::mutate(ID=paste0("COV.", ID))

# files to combine
sweden_edited # For this one, we edit the .fam file since the IDs there are not "normal"
egypt_man_con_edited # DONE EDITING IDs
egypt_man_cas_edited # DONE EDITING IDs
german_cas_edited # DONE EDITING IDs
italy_cas_one_edited # DONE EDITING IDs
italy_cas_two_edited # DONE EDITING IDs
iran_edited # DONE EDITING IDs


merged_data <- rbind(main_pheno, sweden_edited, egypt_man_con_edited, egypt_man_cas_edited,
                     german_cas_edited, italy_cas_one_edited, italy_cas_two_edited, iran_edited)
write.table(merged_data, file = "phenotypes_main_2021_04_09.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



# Update IIDs in fam file for Swedish samples
sweden_fam <- fread("TC-2500_210205.0_removed_GRCh38.fam")
names(sweden_fam) <- c("FID", "IID", "FatherID", "MotherID", "Sex", "Phenotype")
sweden_fam$id <- 1:nrow(sweden_fam) # add id col to the df so we can use it to re-order the df to its original state at the end
sweden_fam <- sweden_fam %>% dplyr::mutate(IID = ifelse(IID < 10, paste0("GB00", IID), paste0("GB0", IID)))

# revert the df to its original order (how the initial fam fail was ordered)
ordered_sweden_fam <- sweden_fam[order(sweden_fam$id), ]

final_sweden_fam <- ordered_sweden_fam %>% select(-id) # remove the id col now that everything is back to its original order

# write the updated fam information into a file
write.table(final_sweden_fam, file = paste("TC-2500_210205.0_removed_GRCh38_updated_ID", ".fam", sep = ""), append = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)



fam_1 <- fread('PLUS_COV_ILLUMINA_19032021.fam')
pheno_1 <- fread('phenotypes_main_2021_04_09.tsv')
setdiff(fam_1$V2, pheno_1$ID)


# 11.04.2021 - Update pheno file to convert age == 0 and sex == 9 to NA
setwd('~/release_15042021/')
main_pheno <- fread('phenotypes_main_2021_04_10.tsv')

main_pheno <- main_pheno %>% 
  mutate(Sex = ifelse(Sex == 9, NA, Sex),
         Age = ifelse(Age == 0, NA, Age))

write.table(main_pheno, file = "phenotypes_main_2021_04_11.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# 11.04.2021 - Update phebo file to remove covid- samples for iran, and set to NA sex mismatches
setwd('~/release_15042021/')
main_pheno <- fread('phenotypes_main_2021_04_11.tsv')

iran_neg <- as.data.frame(read_excel('DDRI-FIMM-data.xlsx', sheet = 'DDRI-FIMM')) %>%
  dplyr::mutate(ID=paste0("COV.", SAMPLE_IDENTIFIER)) %>%
  filter(covid19_test == 0)

main_pheno$Analysis_A1[main_pheno$ID %in% iran_neg$ID] <- NA
main_pheno$Analysis_A2[main_pheno$ID %in% iran_neg$ID] <- 0
main_pheno$Analysis_B1[main_pheno$ID %in% iran_neg$ID] <- NA
main_pheno$Analysis_B2[main_pheno$ID %in% iran_neg$ID] <- 0
main_pheno$Analysis_C2[main_pheno$ID %in% iran_neg$ID] <- 0

write.table(main_pheno, file = "phenotypes_main_2021_04_11_w_mismatch.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

sex_mismatch <- fread('imputed_sex_discordant.tsv') %>%
  mutate(ID = s)

# temporarily set all pheno to NA for samples with sex mismatches
main_pheno$Analysis_A1[main_pheno$ID %in% sex_mismatch$ID] <- NA
main_pheno$Analysis_A2[main_pheno$ID %in% sex_mismatch$ID] <- NA
main_pheno$Analysis_B1[main_pheno$ID %in% sex_mismatch$ID] <- NA
main_pheno$Analysis_B2[main_pheno$ID %in% sex_mismatch$ID] <- NA
main_pheno$Analysis_C2[main_pheno$ID %in% sex_mismatch$ID] <- NA

write.table(main_pheno, file = "phenotypes_main_2021_04_11_wo_mismatch.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# # # UPDATE IRAN

# read in pheno file, put ID as first column, delete "sex" column
main_pheno <- fread("phenotypes_main_2021_04_11_w_mismatch.tsv") %>% 
  select(ID, Cohort, Age:Analysis_C2)


# 16.06.2021: NOT NECESSARY, IDs have been already corrected in the plink file
# # jsut for check
# # main_pheno$ID_corr <- NA
# 
# iran <- main_pheno %>% 
#   filter(Cohort == "iran")
# 
# # correct info for each ID on plate 5, due to sample swapping
# id_map <- fread("iran_new_list_anu.tsv") %>% 
#   mutate(id.new = paste0("COV.", `Corrected ID`)) %>% 
#   select(id.ori = `Sample ID`, id.new)
# 
# for (i in 1:nrow(id_map)) {
#   main_pheno[main_pheno$ID == id_map[i,id.ori], 2:ncol(main_pheno)] <- iran[iran$ID == id_map[i,id.new], 2:ncol(main_pheno)]
#   # main_pheno$ID_corr[main_pheno$ID == id_map[i,id.ori]] <- id_map[i,id.new]
# }


# # # # # # # # # # # # # # # # # # # # # 
# # # 16.06.2021 UPDATE PHENO FILE  # # #
# # # # # # # # # # # # # # # # # # # # #

main_pheno <- fread("phenotypes_main_2021_04_11_w_mismatch.tsv") %>% 
  select(ID, Cohort, Age:Analysis_C2)

# Correct sex for sample A221

main_pheno$Sex[main_pheno$ID == "COV.A221"] <- 2

write.table(main_pheno, file = "phenotypes_main_2021_06_16.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

system('gsutil cp phenotypes_main_2021_06_16.tsv gs://dsge-covid19-data/')


# # temporarily set all pheno to NA for samples with sex mismatches
# sex_mismatch <- fread('imputed_sex_discordant.tsv') %>%
#   mutate(ID = s) %>% 
#   filter(!startsWith(x = ID, prefix = "COV.A"))
# 
# main_pheno$Analysis_A1[main_pheno$ID %in% sex_mismatch$ID] <- NA
# main_pheno$Analysis_A2[main_pheno$ID %in% sex_mismatch$ID] <- NA
# main_pheno$Analysis_B1[main_pheno$ID %in% sex_mismatch$ID] <- NA
# main_pheno$Analysis_B2[main_pheno$ID %in% sex_mismatch$ID] <- NA
# main_pheno$Analysis_C2[main_pheno$ID %in% sex_mismatch$ID] <- NA
# 
# write.table(main_pheno, file = "phenotypes_main_2021_06_16_wo_mismatch.tsv", append = FALSE, sep = "\t", dec = ".",
#             row.names = FALSE, col.names = TRUE, quote = FALSE)
# system('gsutil cp phenotypes_main_2021_06_16_wo_mismatch.tsv gs://dsge-covid19-data/')
