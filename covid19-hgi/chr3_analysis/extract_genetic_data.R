### Read in genetics data and extract SNP

### Extract SNPs

#--- NEW APPROACH ----
## Spain_bujanda

for chr in {1..22}; do
/home/aganna/plink2 --vcf /home/aganna/hgi_spain_bujanda/imputed_data/${chr}.vcf.gz  \
--make-bed --out /home/aganna/hgi_spain_bujanda/temp${chr}
done

rm /home/aganna/hgi_spain_bujanda/mergelist.txt

for chr in {1..22}; do
echo /home/aganna/hgi_spain_bujanda/temp${chr} >> /home/aganna/hgi_spain_bujanda/mergelist.txt
done

# Notice this is in plink 1.9 format
/home/aganna/plink --merge-list /home/aganna/hgi_spain_bujanda/mergelist.txt --make-bed --out /home/aganna/hgi_spain_bujanda/spain_bujanda_imputed_all

/home/aganna/plink2 --bfile /home/aganna/hgi_spain_bujanda/spain_bujanda_imputed_all  \
--extract bed1 /home/aganna/list_snp_to_extract  --export A --out /home/aganna/hgi_spain_bujanda/spain_bujanda_all_snps

rm /home/aganna/hgi_spain_bujanda/temp*
  

#---- OLD APPROACH -----

## Brasil
/home/aganna/plink2 --bfile /home/aganna/hgi_brasil/chr3_brazil.dose \
--snps "chr3:45848429:A:T,chr3:45834967:G:GA,chr9:136132908:T:TC,chr3:101434448:T:G,chr21:33242905:T:C,chr19:10355447:C:T,chr19:4723658:C:A,chr6:41535823:G:A,chr6:31153649:G:A,chr12:112944151:C:T" --export A --out /home/aganna/hgi_brasil/chr3_brazil_all_snps

## Italy-Belgium-Sweden
/home/aganna/plink2 --bfile /home/aganna/chr3_ita_be_swe.dose \
--snps "chr3:45848429:A:T,chr3:45834967:G:GA,chr9:136132908:T:TC,chr3:101434448:T:G,chr21:33242905:T:C,chr19:10355447:C:T,chr19:4723658:C:A,chr6:41535823:G:A,chr6:31153649:G:A,chr12:112944151:C:T" --export A --out /home/aganna/chr3_ita_be_swe_rs35081325

## Spain_alarcon
/home/aganna/plink2 --bfile /home/aganna/hgi_spain_alarcon/SPGRX.Alarcon.chr3.20mks.GRCh38 \
--snps "3:45848429:A:T" --export A --out /home/aganna/hgi_spain_alarcon/hgi_spain_alarcon_rs35081325_part1

/home/aganna/plink2 --bfile /home/aganna/hgi_spain_alarcon/SPGRX.Alarcon_CHR3region_80samples/chr3reg20MKsGRCh38 \
--snps "3:45848429:A:T" --export A --out /home/aganna/hgi_spain_alarcon/hgi_spain_alarcon_rs35081325_part2


## Canada -- Does not exists
/home/aganna/plink2 --bfile /home/aganna/hgi_canada/bqc19 \
--snps "rs35081325" --export A --out /home/aganna/hgi_canada/hgi_canada_rs35081325

## Italy_valenti
/home/aganna/plink2 --vcf /home/aganna/hgi_italy_valenti/imputed_data/3.vcf.gz \
--snps "chr3:45848429:A:T,chr3:45834967:G:GA,chr9:136132908:T:TC,chr3:101434448:T:G,chr21:33242905:T:C,chr19:10355447:C:T,chr19:4723658:C:A,chr6:41535823:G:A,chr6:31153649:G:A,chr12:112944151:C:T" --export A --out /home/aganna/hgi_italy_valenti/italy_valenti_rs35081325

## Spain_butti
/home/aganna/plink2 --vcf /home/aganna/hgi_spain_butti/imputed_data/3.vcf.gz \
--snps "chr3:45848429:A:T,chr3:45834967:G:GA,chr9:136132908:T:TC,chr3:101434448:T:G,chr21:33242905:T:C,chr19:10355447:C:T,chr19:4723658:C:A,chr6:41535823:G:A,chr6:31153649:G:A,chr12:112944151:C:T" --export A --out /home/aganna/hgi_spain_butti/spain_butti_all_snps


library(data.table)
library(xlsx)
library(ggplot2)
library(readr)
library(stringr)

af <- function(geno)
{
  ## calc_n
  n0 <- sum(geno==0,na.rm=T)
  n1 <- sum(geno==1,na.rm=T)
  n2 <- sum(geno==2,na.rm=T)
  n <- n0 + n1 + n2
  ## calculate allele frequencies
  p <- ((2*n0)+n1)/(2*n)
  return(p)
}
  

path <- "/home/aganna/"

# Read manifests with sample ID
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-250920_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)
swe_manifest <- read.xlsx(paste0(path,"hgi_swe/PronMed genetic fenotypes.xlsx"),sheetName = "Data",stringsAsFactors=FALSE)
spain_bujanda_manifest <- read.xlsx(paste0(path,"hgi_spain_bujanda/international_database_Spanish_Bujanda_v5_editedTN.xlsx"),sheetName = "Spanish patients_One-visit",stringsAsFactors=FALSE,colClasses="character")
spain_bujanda_manifest <- spain_bujanda_manifest %>% filter(!is.na(anonymized_patient_id))
spain_bujanda_key <- fread(paste0(path,"hgi_spain_bujanda/key_geno_pheno.csv"))
ita_renieri_manifest <- read_tsv(paste0(path,"hgi_italy/Italy.withCom_20200925TN.tsv"))
bra_manifest <- read.xlsx(paste0(path,"hgi_brasil/Data_HGI_Chr3_manuscript.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
can_manifest <- read_tsv(paste0(path,"hgi_canada/bqc19_one_time_V2.tsv"))
ita_valenti_manifest <- read.xlsx(paste0(path,"hgi_italy_valenti/Clinical_Data_V2_Milan.xlsx"),sheetName = "One_Time",stringsAsFactors=FALSE,colClasses="character")


## Read ethnicity assignment and PCs
pcs_eth_swe_ita_bel <- fread(paste0(path,"cases_controls_pca_sup_pops_0.5_probs_sweitalbel.txt"))
pcs_eth_brasil <- fread(paste0(path,"hgi_brasil/pca/cases_controls_pca_sup_pops_0.8_probs.txt"))
pcs_eth_spain_bujanda <- fread(paste0(path,"hgi_spain_bujanda/pca/cases_controls_pca_sup_pops_0.8_probs.txt"))
pcs_eth_canada <- fread(paste0(path,"hgi_canada/pca/cases_controls_pca_sup_pops_0.8_probs.txt"))


# Read genetic data
swe_ita_bel_rs35081325 <- fread(paste0(path,"chr3_ita_be_swe_rs35081325.raw"))
colnames(swe_ita_bel_rs35081325) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
swe_ita_bel_rs11385942 <- fread(paste0(path,"chr3_ita_be_swe_rs11385942.raw")) %>% select("chr3:45834967:G:GA_G")
colnames(swe_ita_bel_rs11385942) <- c("rs11385942")
swe_ita_bel <- bind_cols(swe_ita_bel_rs35081325,swe_ita_bel_rs11385942)

brasil_gene_rs35081325 <- fread(paste0(path,"hgi_brasil/chr3_brazil_rs35081325.raw"))
colnames(brasil_gene_rs35081325) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
brasil_gene_rs11385942 <- fread(paste0(path,"hgi_brasil/chr3_brazil_rs11385942.raw")) %>% select("chr3:45834967:G:GA_G")
colnames(brasil_gene_rs11385942) <- c("rs11385942")
brasil_gene <- bind_cols(brasil_gene_rs35081325,brasil_gene_rs11385942)


spain_bujanda_gene_rs35081325 <- fread(paste0(path,"hgi_spain_bujanda/spain_bujanda_rs35081325.raw"))
colnames(spain_bujanda_gene_rs35081325) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
spain_bujanda_gene_rs11385942 <- fread(paste0(path,"hgi_spain_bujanda/spain_bujanda_rs11385942.raw")) %>% select("chr3:45834967:G:GA_G")
colnames(spain_bujanda_gene_rs11385942) <- c("rs11385942")
spain_bujanda_gene <- bind_cols(spain_bujanda_gene_rs35081325,spain_bujanda_gene_rs11385942)


spain_alarcon_gene_part1 <- fread(paste0(path,"hgi_spain_alarcon/hgi_spain_alarcon_rs35081325_part1.raw"))
colnames(spain_alarcon_gene_part1) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
spain_alarcon_gene_part2 <- fread(paste0(path,"hgi_spain_alarcon/hgi_spain_alarcon_rs35081325_part2.raw"))
colnames(spain_alarcon_gene_part2) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
spain_alarcon_gene_rs35081325 <- rbind(spain_alarcon_gene_part1,spain_alarcon_gene_part2)

spain_alarcon_gene_part1 <- fread(paste0(path,"hgi_spain_alarcon/hgi_spain_alarcon_rs11385942_part1.raw")) %>% select("3:45834967:G:GA_G")
colnames(spain_alarcon_gene_part1) <- c("rs11385942")
spain_alarcon_gene_part2 <- fread(paste0(path,"hgi_spain_alarcon/hgi_spain_alarcon_rs11385942_part2.raw")) %>% select("3:45834967:G:GA_G")
colnames(spain_alarcon_gene_part2) <- c("rs11385942")
spain_alarcon_gene_rs11385942 <- rbind(spain_alarcon_gene_part1,spain_alarcon_gene_part2)

spain_alarcon_gene <- bind_cols(spain_alarcon_gene_rs35081325,spain_alarcon_gene_rs11385942)


canada_gene_rs35081325 <- fread(paste0(path,"hgi_canada/rs35081325.raw"))
colnames(canada_gene_rs35081325) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
canada_gene_rs11385942 <- fread(paste0(path,"hgi_canada/rs11385942.raw")) %>% select("chr3:45834967:G:GA_GA")
colnames(canada_gene_rs11385942) <- c("rs11385942")
canada_gene <- bind_cols(canada_gene_rs35081325,canada_gene_rs11385942) %>% mutate(rs35081325=2-round(canada_gene$rs35081325,0),rs11385942=2-round(canada_gene$rs35081325,0)) # INVERT BECAUSE THEY USE THE OPPOSITE ALLELE AND USE DOSAGES


italy_valenti_gene_rs35081325 <- fread(paste0(path,"hgi_italy_valenti/italy_valenti_rs35081325.raw"))
colnames(italy_valenti_gene_rs35081325) <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE","rs35081325")
italy_valenti_gene_rs11385942 <- fread(paste0(path,"hgi_italy_valenti/italy_valenti_rs11385942.raw")) %>% select("chr3:45834967:G:GA_G")
colnames(italy_valenti_gene_rs11385942) <- c("rs11385942")
italy_valenti_gene <- bind_cols(italy_valenti_gene_rs35081325,italy_valenti_gene_rs11385942)




# Select Swedish samples
temp <- pcs_eth_swe_ita_bel %>% mutate(anonymized_patient_id=as.character(as.numeric(gsub("COV.COV-","",s))))
swe_gene_pp <- swe_ita_bel %>% mutate(anonymized_patient_id=as.character(as.numeric(gsub("0_COV.COV-","",IID)))) %>% filter(anonymized_patient_id %in% as.character(swe_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="sweden") %>% left_join(temp,by=c("anonymized_patient_id"="anonymized_patient_id"))  %>% select(-s)

# Select Belgium samples
temp <- pcs_eth_swe_ita_bel %>% mutate(anonymized_patient_id=as.character(gsub("COV.","",s)))
bel_gene_pp <- swe_ita_bel %>% mutate(anonymized_patient_id=as.character(gsub("0_COV.","",IID))) %>% filter(anonymized_patient_id %in% as.character(bel_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="belgium") %>% left_join(temp,by=c("anonymized_patient_id"="anonymized_patient_id"))  %>% select(-s)

# Select Italian renieri samples
temp <- pcs_eth_swe_ita_bel %>% mutate(anonymized_patient_id=gsub("_","-",as.character(gsub("^COV.","",pcs_eth_swe_ita_bel$s))))
ita_renieri_gene_pp <- swe_ita_bel %>% mutate(anonymized_patient_id=gsub("_","-",as.character(gsub("^0_COV.","",IID)))) %>% filter(anonymized_patient_id %in% as.character(ita_renieri_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="italy_renieri") %>% left_join(temp,by=c("anonymized_patient_id"="anonymized_patient_id"))  %>% select(-s)

# Select Spanish_bujanda samples
temp <- pcs_eth_spain_bujanda %>% mutate(anonymized_patient_id=as.character(s)) %>% select(-s)
spain_bujanda_gene_pp <- spain_bujanda_gene %>% left_join(pcs_eth_spain_bujanda,by=c("IID"="s")) %>% inner_join(spain_bujanda_key,by=c("IID"="Tube number")) %>% mutate(anonymized_patient_id=as.character(CIC)) %>% filter(anonymized_patient_id %in% as.character(spain_bujanda_manifest$anonymized_patient_id)) %>% mutate(study="spain_bujanda") %>% select(-c("FID","IID","PAT","MAT","SEX","PHENOTYPE","CIC"))


# Select canada samples
temp <- pcs_eth_canada %>% mutate(anonymized_patient_id=as.character(s)) %>% select(-s)
canada_gene_pp <- canada_gene %>% mutate(anonymized_patient_id=as.character(as.numeric(str_sub(IID, start=1, end=4)))) %>% filter(anonymized_patient_id %in% as.character(can_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="canada") %>% left_join(temp,by=c("anonymized_patient_id"="anonymized_patient_id"))


# Select brasil samples
temp <- pcs_eth_brasil %>% mutate(anonymized_patient_id=str_split(s, "_",simplify=TRUE)[,1]) %>% select(-s)
brasil_gene_pp <- brasil_gene %>% mutate(anonymized_patient_id=str_split(IID, "_",simplify=TRUE)[,2]) %>% filter(anonymized_patient_id %in% as.character( brasil_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="brasil") %>% left_join(temp,by=c("anonymized_patient_id"="anonymized_patient_id"))


# Select Italian valenti samples -- POPULATION ASSIGNMENT IN TEMPORARY
ita_valenti_gene_pp <- italy_valenti_gene %>% mutate(anonymized_patient_id=IID) %>% filter(anonymized_patient_id %in% as.character(italy_valenti_manifest$anonymized_patient_id)) %>% select("anonymized_patient_id","rs35081325","rs11385942","SEX_GENE"="SEX") %>% mutate(study="italy_valenti") %>% mutate(pop="EUR")


# Combine
gene_all <- bind_rows(list(swe_gene_pp,bel_gene_pp,ita_renieri_gene_pp,spain_bujanda_gene_pp,canada_gene_pp,brasil_gene_pp,ita_valenti_gene_pp))


# Plot allele frequencies
df <- gene_all %>% filter(pop=="EUR") %>% select(rs35081325,anonymized_patient_id,study) %>% pivot_longer(!c(anonymized_patient_id,study))
df <- df %>% group_by(study,name) %>% summarise(freq=af(value))
ggplot(df,aes(x=study,y=freq*100)) + geom_bar(stat="identity") + facet_wrap(~name, ncol=7) + theme_bw() + ylab("Variant Frequency") + theme(axis.text.x=element_text(angle=45, hjust=1)) 


df <- gene_all %>% filter(pop=="EUR") %>% select(rs11385942,anonymized_patient_id,study) %>% pivot_longer(!c(anonymized_patient_id,study))
df <- df %>% group_by(study,name) %>% summarise(freq=af(value))
ggplot(df,aes(x=study,y=freq*100)) + geom_bar(stat="identity") + facet_wrap(~name, ncol=7) + theme_bw() + ylab("Variant Frequency") + theme(axis.text.x=element_text(angle=45, hjust=1)) 

## Correlation between the two variables
gene_all %>% group_by(study) %>% summarize(cor(rs11385942, rs35081325, use="complete.obs"))

## Write out combined genetic data
gene_all  %>% mutate(temp=recode(study, belgium="BB",sweden="SW",spain_bujanda="SB",spain_butti="SU",italy_renieri="ITR",brasil="BR",canada="CA",italy_valenti="ITV"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% select(-temp) %>% write_tsv(paste0(path,"all_variant_rs35081325_rs11385942.tsv"),col_names = TRUE)
