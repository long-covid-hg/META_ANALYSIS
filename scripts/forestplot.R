setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/07.LongCOVID_meta/DF2")
library(meta)
library(grid)
library(meta)
library(grid)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

#data loading  IT is faster if you pre-select the variant of your interest.
data <- fread(paste0("W2.3_meta_out.chr6.LD.tsv"))
data <- data %>% dplyr::select(-colnames(data)[grepl("leave", colnames(data))])

for(i in c(1:dim(data)[1])){
  tmp <- data[i,] %>% dplyr::select("#CHR","POS", "REF","ALT","SNP",
                                colnames(data)[grepl("_beta|_sebeta|_pval|_af_alt|_INFO|_Nsample", colnames(data))]) %>%
    dplyr::select(-any_of(colnames(data)[grepl("^all_|^lmso_", colnames(data))]))
  
  tmp_long <- reshape(tmp, direction='long', 
                      varying=colnames(tmp)[6:dim(tmp)[2]], 
                      timevar='cohort',
                      times=unique(paste0(str_split(colnames(tmp)[6:dim(tmp)[2]], pattern="\\_", simplify = T)[,1],"_",str_split(colnames(tmp)[6:dim(tmp)[2]], pattern="\\_", simplify = T)[,2])),
                      v.names=c('INFO', 'Nsample',  'af_alt', 'beta', 'pval','sebeta' ),
                      idvar='SNP')
  
  colnames(tmp_long)[7:12] <- c( 'beta', 'sebeta', 'pval','af_alt', 'INFO', 'Nsample')
  tmp_long <- tmp_long %>% arrange(cohort)
  m1 <- metagen(beta,
                sebeta,
                data=tmp_long,
                studlab=paste(cohort),
                fixed = TRUE,
                random = TRUE,
                prediction=FALSE,
                sm="OR", control=list(maxiter=1000))
  png(paste0("chr6/forestplot_W2.3_",tmp$SNP,"_",tmp$rsid,".png"), width=500, height=500)
  forest(m1, smlab="",leftcols=c("studlab"),
         rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
         leftlabs = c("Cohorts"),zero.pval = TRUE)
  grid.text(paste0(tmp$SNP," ",tmp$rsid), .5, 0.95, gp=gpar(cex=1.2))
  dev.off()
}




