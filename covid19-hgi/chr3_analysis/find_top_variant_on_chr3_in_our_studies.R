library(data.table)
library(meta)


d <- fread("B2_ALL_inv_var_meta")

d_top <- d[d$'#CHR'==3 & (d$POS<(45834967+100000) & d$POS>(45834967-100000)),]

RES <- NULL
for(i in 1:nrow(d_top))
{
  mod <- metagen(c(d_top$HOSTAGE_EUR_beta[i],d_top$GENCOVID_EUR_beta[i],d_top$BelCovid_EUR_beta[i],d_top$SweCovid_EUR_beta[i]),c(d_top$HOSTAGE_EUR_sebeta[i],d_top$GENCOVID_EUR_sebeta[i],d_top$BelCovid_EUR_sebeta[i],d_top$SweCovid_EUR_sebeta[i]))
  
  RES <- rbind(RES,c(d_top$SNP[i],mod$TE.fixed,mod$pval.fixed))
  print(i)
}

head(RES[order(as.numeric(RES[,3])),])
