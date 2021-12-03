#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
NR >1 {print $a["CHR"],$a["BP"],$a["Allele1"],$a["Allele2"],$a["AF_Allele2"],$a["imputationInfo"],$a["BETA"],$a["SE"],$a["P"],$a["N"]}' | \
bgzip -@4 > $1.formatted.gz
