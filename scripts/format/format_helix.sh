#!/usr/bin/env bash

m=$(zcat $1 | awk 'NR >1 {if(m < ($14+$18)) m=($14+$18)}END{print m}')
gunzip -c $1 | awk '
BEGIN {FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
NR >1 {print $a["Chr"],$a["Pos"],$a["Ref"],$a["Alt"],$a["AAF"],($a["Num_Cases"]+$a["Num_Controls"])/'$m',log($a["Effect"]),(log($a["UCI_Effect"])-log($a["LCI_Effect"]))/3.92,$a["Pval"],($a["Num_Cases"]+$a["Num_Controls"])}' | \
awk '$9 != "NA"' | \
bgzip -@4 > $1.formatted.gz
