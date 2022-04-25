#!/bin/bash
#

# store input argument as task ID
#taskid=$1

# get list of munged files - will loop over all task IDs inputted
# if a sumstats file is listed more than once (i.e. has been munged multiple times)
# the awk command will pick the last occurrence
{
   for taskid in "$@"
   do
      gsutil ls gs://long-covid-hg-cromwell/format/$taskid/call-formatting/** | grep ".\.gz$"
   done
} | awk '{split($1,a,"/");c[a[length(a)]]=$1}END{for(filen in c){print c[filen]}}' > formatted_sumstats.txt

# generate ethnicity column - defaults to "oth" if not "EUR", "FIN", "NFE", "SAS", "EAS", "AFR" and "AMR"
awk '{split($1,a,"/");split(a[length(a)],b,"_");eth="oth"}b[7]=="EUR"||b[7]=="NFE"{eth="nfe"}tolower(b[1])~/finngen/||b[7]=="FIN"{eth="fin"}b[7]=="SAS"||b[7]=="EAS"||b[7]=="AFR"||b[7]=="AMR"{eth=tolower(b[7])}{print $1"\t"eth}' formatted_sumstats.txt > data/DF2/step2_munge.txt

# delete intermediate file
rm formatted_sumstats.txt

