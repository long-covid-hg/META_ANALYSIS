#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
NR >1 {b=$a[""]print $a["Chr"],$a["Pos"],$a["Ref"],$a["Alt"],$a["AAF"],$a["Info"],$a["BETA"],$a["SE"],10^-$a["LOG10P"],$a["N"]}' | \
bgzip -@4 > $1.formatted.gz
