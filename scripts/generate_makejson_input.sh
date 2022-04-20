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
      gsutil ls gs://long-covid-hg-cromwell/munge/$taskid/call-lift_postprocess/shard-*/** | grep ".\.gz$"
   done
} | awk '{split($1,a,"/");c[a[length(a)]]=$1}END{for(filen in c){print c[filen]}}' > munged_sumstats.txt

# create input file for makejson
{
   # header
   echo -e "pheno\tname\tfile\tn_cases\tn_controls\tchr\tpos\tref\talt\taf_alt\tinfo\teffect\teffect_type\tpval\tse\tNsamples"
   # NQs and NEs
   awk 'BEGIN{OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}pheno~/^N/{print pheno,name_eth,$1,ncases,ncontrols,"#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta\tN"}' munged_sumstats.txt | sort -k1
   # WQ1.3 (including NQ1.3 if no WQ1.3)
   awk 'BEGIN{OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}pheno=="WQ1.3"{line[name_eth]="WQ1.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}pheno=="NQ1.3"&&(!(name_eth in line)){line[name_eth]="WQ1.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
   # WQ2.3 (including NQ2.3 if no WQ2.3)
   awk 'BEGIN{OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}pheno=="WQ2.3"{line[name_eth]="WQ2.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}pheno=="NQ2.3"&&(!(name_eth in line)){line[name_eth]="WQ2.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
   # N1.3 - prioritise NE1.3 over NQ1.3
   awk 'BEGIN{OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}pheno=="NE1.3"{line[name_eth]="N1.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}pheno=="NQ1.3"&&(!(name_eth in line)){line[name_eth]="N1.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
   # N2.3 - prioritise NE2.3 over NQ2.3
   awk 'BEGIN{OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}pheno=="NE2.3"{line[name_eth]="N2.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}pheno=="NQ2.3"&&(!(name_eth in line)){line[name_eth]="N2.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
   # W1.3 prioritise "E" over "Q" and "W" over "N"
   awk 'BEGIN{p["WE1.3"]=4;p["WQ1.3"]=3;p["NE1.3"]=2;p["NQ1.3"]=1;OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}{ln="W1.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}!(pheno in p){p[pheno]=5}!(name_eth in pp){pp[name_eth]=5}p[pheno]<pp[name_eth]{line[name_eth]=ln}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
   # W2.3 prioritise "E" over "Q" and "W" over "N"
   awk 'BEGIN{p["WE2.3"]=4;p["WQ2.3"]=3;p["NE2.3"]=2;p["NQ2.3"]=1;OFS="\t"}{split($1,a,"/");split(a[length(a)],b,"_")}{pheno=b[3];name=b[1];eth=b[7];name_eth=b[1]"_"b[7];ncases=b[8];ncontrols=b[9]}{ln="W2.3\t"name_eth"\t"$1"\t"ncases"\t"ncontrols"\t#CHR\tPOS\tREF\tALT\taf_alt\timputationInfo\tbeta\tbeta\tpval\tsebeta"}!(pheno in p){p[pheno]=5}!(name_eth in pp){pp[name_eth]=5}p[pheno]<pp[name_eth]{line[name_eth]=ln}END{for(ne in line){print line[ne]}}' munged_sumstats.txt
} > data/DF2/config_meta_F2.tsv

# clean up
rm munged_sumstats.txt
