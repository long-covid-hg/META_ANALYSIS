#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
NR >1 {print $a["CHROM"],$a["GENPOS"],$a["ALLELE0"],$a["ALLELE1"],$a["A1FREQ"],$a["INFO"],$a["BETA"],$a["SE"],10^-$a["LOG10P"],$a["N"]}' | \
bgzip -@4 > $1.formatted.gz
