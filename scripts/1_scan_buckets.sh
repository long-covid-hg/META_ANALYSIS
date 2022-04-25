#!/bin/bash
#
# Samuel E. Jones on behalf of the LongCOVID HGI
# March 2022
#
# This script takes an input list of buckets, builds an inventory
# and generates a report of the data uploaded to each
#

# get script name
scriptname="$(basename "$(test -L "$0" && readlink "$0" || echo "$0")")"
blank=${scriptname//?/ }

# store number of input args
nargs=$#

# echo empty line
echo

#############
### USAGE ###
#############
usage()
{
   # Display Help
   echo -e "\tThis script takes an input list of buckets, builds an inventory of each bucket and"
   echo -e "\tgenerates a report of which studies have uploaded data for which build."
   echo
   echo -e "\tUsage: $scriptname -b BUCKETLIST [-h]"
   echo
   echo -e "\tOptions:"
   echo
   echo -e "\t-b BUCKETLIST   Input text file containing a list of all Long COVID HGI buckets to be"
   echo -e "\t                scanned diagnosis codes. File should be plain-text and bucket addresses"
   echo -e "\t                should start with \"gs://\"."
   echo
   echo -e "\t-h              [OPTIONAL] Prints the script usage and exits. Ignores all other options"
   echo -e "\t                and does not accept arguments."
   echo
   exit $1
}

###################################
### LOGGING AND PRINT FUNCTIONS ###
###################################

# set spaces to pad out subsequent lines printed by "print_" commands
sp11="           "

__log_init__() {
    if [[ -t 1 ]]; then
        # colors for logging in interactive mode
        [[ $COLOR_BOLD ]]   || COLOR_BOLD="\033[1m"
        [[ $COLOR_RED ]]    || COLOR_RED="\033[0;31m"
        [[ $COLOR_GREEN ]]  || COLOR_GREEN="\033[0;34m"
        [[ $COLOR_YELLOW ]] || COLOR_YELLOW="\033[0;33m"
        [[ $COLOR_BLUE ]]   || COLOR_BLUE="\033[0;32m"
        [[ $COLOR_OFF ]]    || COLOR_OFF="\033[0m"
    else
        # no colors to be used if non-interactive
        COLOR_RED= COLOR_GREEN= COLOR_YELLOW= COLOR_BLUE= COLOR_OFF=
    fi
    readonly COLOR_RED COLOR_GREEN COLOR_YELLOW COLOR_BLUE COLOR_OFF

    #
    # map log level strings (FATAL, ERROR, etc.) to numeric values
    #
    # Note the '-g' option passed to declare - it is essential
    #
    unset _log_levels _loggers_level_map
    declare -gA _log_levels _loggers_level_map
    _log_levels=([FATAL]=0 [ERROR]=1 [WARN]=2 [INFO]=3 [DEBUG]=4 [VERBOSE]=5)

}

print_error() {
    {
        printf "${COLOR_RED}ERROR ::   "
        printf '%s\n' "$1"
	[ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
        printf "$COLOR_OFF\n"
    } >&2
    exit 1
}

print_warn() {
    printf "${COLOR_YELLOW}WARNING :: "
    printf '%s\n' "$1"
    [ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
    printf "$COLOR_OFF"
}

print_info() {
    printf "${COLOR_BLUE}INFO ::    "
    printf '%s\n' "$1"
    [ "$#" -gt "1" ] && { shift; printf "${sp11}%s\n" "$@"; }
    printf "$COLOR_OFF"
}

# print only if output is going to terminal
print_tty() {
    if [[ -t 1 ]]; then
        printf "${sp11}%s\n" "$@"
    fi
}

################################
### INPUT CHECKING FUNCTIONS ###
################################

# define function to catch badly parsed arguments
fregex="^-[A-Za-z]$"
catch_badargs() {
   [[ "$1" =~ $fregex ]] && print_error "Missing argument for \"-$2\" option."
}

# function to check if file exists and has non-zero size
check_file_exists() {
   [ ! -s "$1" ] && print_error "Input file \"$1\" is missing or empty."
}

# check file is the right shape
check_file_columns() {
   badline=`awk -v nf=$2 'NF!=nf{print NR;exit}' $1`
   if [ ! -z "$badline" ]
   then
      badlinefields=`awk -v bl=$badline 'NR==bl{print NF}' $1`
      [ "$2" -gt "1" ] && print_error "File \"$1\" is malformed on line $badline. Expected $2 fields but found $badlinefields." || print_error "File \"$1\" is malformed on line $badline. Expected $2 field but found $badlinefields."
   fi
}


#######################
### OPTIONS PARSING ###
#######################

# if no arguments passed, display usage and exit with code 1
if [[ ! $@ =~ ^\-.+ ]]
then
   print_error "No options or arguments. Run \"./$scriptname -h\" to see usage."
fi

# parse arguments
while getopts ":hb:" option; do
   case "$option" in
      h) # display usage
         usage 0
         ;;
      b) # input code list
         bucketlist="$OPTARG"
         catch_badargs $OPTARG b
         ;;
      :) # missing argument for option that requires an argument
         print_error "Missing argument for \"-$OPTARG\" option."
         ;;
     \?) # unrecognised option
         print_error "Invalid option chosen. Run \"./$scriptname -h\" to see usage."
         ;;
   esac
done


########################
### CHECK INPUT FILE ###
########################

# check specified bucket list exists
check_file_exists $bucketlist

# check only one column exists for file
check_file_columns $bucketlist 1

# check all lines start with gs://
badline=`awk '$0!~/^gs:\/\//{print NR;exit}' $bucketlist`
[ ! -z "$badline" ] && print_error "Bucket list file \"$bucketlist\" malformed one line $badline."


####################
### SCAN BUCKETS ###
####################

# print info message as this step can take a minute or so
print_info "Scanning buckets listed in \"$bucketlist\" - this may take a few minutes"

# loop over all buckets and dump inventory into file
buckets=`awk '{print $1}' $bucketlist`
bcount=0
{
   for bucket in $buckets
   do
      # run ls -lh command on bucket and reformat to have filename then datestamp then size
      gsutil ls -lh $bucket | awk '$1~/^TOTAL/{next}NF==1{str=$1;gsub(/^[ \t]+/,"",str);print str;next}{print $4"\t"$3"\t"$1$2}'
      bcount=$((bcount+1))
   done
} > .full_inventory_all_buckets.txt
# get count of non-empty buckets
nnonemptybuckets=`awk '{split($1,a,"/");print a[3]}' .full_inventory_all_buckets.txt | sort | uniq | wc -l`
print_info "Bucket scans show $nnonemptybuckets buckets containing data out of $bcount buckets listed"


############################################
### CHECK FOR REQUESTED FILES IN BUCKETS ###
############################################

### check files matching the naming format
# specify regex segments
dsregex="[[:alnum:].-]+"            #[dataset]
snregex="[[:alpha:].-]+"             #[surname]
pnregex="[NW][QE][12][.][369]"        #[phenotype]
dfregex="F[1-9][0-9]*"                #[freeze_number]
aregex="(ALL|LE_60|GT_60)"            #[age]
sregex="(ALL|FEMALE|MALE)"            #[sex]
eregex="[[:alpha:]]{3}"               #[ancestry]
caregex="[0-9]+"                      #[n_cases]
coregex="[0-9]+"                      #[n_controls]
gwregex="(REGENIE|SAIGE)"             #[gwas software]
dregex="202[1-9][0-1][0-9][0-3][0-9]" #[YYYYMMDD]

# define regex for GWAS results file
# requested format: [dataset]_[surname]_[phenotype]_[freeze_number]_[age]_[sex]_[ancestry]_[n_cases]_[n_controls]_[gwas software]_[YYYYMMDD].txt.gz
nregex="^${dsregex}_${snregex}_${pnregex}_${dfregex}_${aregex}_${sregex}_${eregex}_${caregex}_${coregex}_${gwregex}_${dregex}[.]txt(.gz)?$"

# define regex for sample summary file
# requested format: [dataset]_[surname]_[freeze_number].zip
mregex="^${dsregex}_${snregex}_F[1-9][0-9]*[.]zip$"

# clear temporary files for next section (if they already exist)
rm -f .gwas_results_all_buckets.txt .gwas_results_all_buckets_filtered.txt .meta_files_all_buckets.txt .other_files_all_buckets.txt

# extract filenames that match the regex string
for line in `awk '{print $0}' .full_inventory_all_buckets.txt`
do
   # break line into fields and get field count
   fullfilen=`echo $line | awk '{print $1}'`
   datestamp=`echo $line | awk '{print $2}'`
   filesize=`echo $line | awk '{print $3}'`
   numfields=`echo $line | awk '{print NF}'`
   # get name of file without bucket
   filen=`echo $fullfilen | awk '{split($1,a,"/");print a[length(a)]}'`
   if [[ "$filen" =~ $nregex ]]
   then
      # write to gwas results file if matching GWAS results regex
      echo -e "$fullfilen\t$datestamp\t$filesize" >> .gwas_results_all_buckets.txt
   elif [[ "$filen" =~ $mregex ]]
   then
      # write to sample summary file if matching sample summary regex
      echo -e "$fullfilen\t$datestamp\t$filesize" >> .meta_files_all_buckets.txt
   else
      # else write to other file list
      [ "${numfields}" -eq "1" ] && echo $fullfilen >> .other_files_all_buckets.txt || echo -e "$fullfilen\t$datestamp\t$filesize" >> .other_files_all_buckets.txt
   fi
done

# now generate csv showing which studies provide results for which phenotypes
# (build list of studies from recognised GWAS results and sample summary files)
{
   echo "STUDY,BUCKET"
   awk '{split($1,a,"/");split(a[4],b,"_");print b[1]","a[3]}' .gwas_results_all_buckets.txt .meta_files_all_buckets.txt | sort -k1 | uniq
} > .study_contribution_table.csv

for fno in {1..20} # assume 20 max data freezes?
do
   # check if any results uploaded for freeze $fno
   if grep -q "_F${fno}_" .gwas_results_all_buckets.txt
   then
      # if so, check for the sample summary file
      awk -v fno=$fno 'BEGIN{FS=","}NR==FNR{a[$1]++;next}FNR==1{print $0",F"fno"_sampleinfo";next}{cc=""}($1 in a){cc="YES"}{print $0","cc}' <(grep "_F${fno}.zip" .meta_files_all_buckets.txt | awk 'BEGIN{print "DUMMY"}{split($1,a,"/");split(a[4],b,"_");print b[1]}') .study_contribution_table.csv > .tmp.csv && mv .tmp.csv .study_contribution_table.csv
      # then check all phenotypes
      for phe in WQ1.3 WQ2.3 NQ1.3 NQ2.3 WE1.3 WE2.3 NE1.3 NE2.3
      do
	 # if results for the phenotype exist, add the relevant column to the csv
	 if grep -qF "_${phe}_" <(grep "_F${fno}_" .gwas_results_all_buckets.txt)
	 then
	    # if more than one match per cohort/freeze/pheno, take one with latest timestamp
	    awk -v fno=$fno -v phe=$phe 'BEGIN{FS=","}NR==FNR{a[$1]=$2;next}FNR==1{print $0",F"fno"_"phe;next}{cc=""}($1 in a){cc=a[$1]}{print $0","cc}' <(grep "_F${fno}_" .gwas_results_all_buckets.txt | grep -F "_${phe}_" | awk 'function datenum(x){split(x,xx,"T");split(xx[1],y,"-");dn=((y[1]-2000)*366)+((y[2]-1)*31)+y[3];return dn}function timenum(x){split(x,xx,"T");split(xx[2],y,":");gsub("Z","",y[3]);tn=(3600*y[1])+(60*y[2])+y[3];return tn}{split($1,a,"/");split(a[4],b,"_")}!(b[1] in cc){cc[b[1]]=b[8]"/"b[9];d[b[1]]=datenum($3);t[b[1]]=timenum($3);next}(datenum($3)<d[b[1]]||(datenum($3)==d[b[1]]&&timenum($3)<t[b[1]])){cc[b[1]]=b[8]"/"b[9];d[b[1]]=datenum($3);t[b[1]]=timenum($3)}END{for(study in cc){print study","cc[study]}}') .study_contribution_table.csv > .tmp.csv && mv .tmp.csv .study_contribution_table.csv
	    # then add these files to the GWAS input list
	    grep "_F${fno}_" .gwas_results_all_buckets.txt | grep -F "_${phe}_" | awk 'function datenum(x){split(x,xx,"T");split(xx[1],y,"-");dn=((y[1]-2000)*366)+((y[2]-1)*31)+y[3];return dn}function timenum(x){split(x,xx,"T");split(xx[2],y,":");gsub("Z","",y[3]);tn=(3600*y[1])+(60*y[2])+y[3];return tn}{split($1,a,"/");split(a[4],b,"_")}!(b[1] in fn){fn[b[1]]=$1;d[b[1]]=datenum($3);t[b[1]]=timenum($3)}(datenum($3)<d[b[1]]||(datenum($3)==d[b[1]]&&timenum($3)<t[b[1]])){fn[b[1]]=$1;d[b[1]]=datenum($3);t[b[1]]=timenum($3)}END{for(study in fn){print fn[study]}}' >> .gwas_results_all_buckets_filtered.txt
         fi
      done
   fi
done


#############################
### GENERATE REPORT FILES ###
#############################

# get date as string
curdate=`date +'%Y.%m.%d'`

# create subdirectory for reports
mkdir -p Bucket_scan/$curdate/

# move pre-generated table
mv .study_contribution_table.csv Bucket_scan/$curdate/study_contribution_table_$curdate.csv

# move list of non-recognised files
mv .other_files_all_buckets.txt Bucket_scan/$curdate/unrecognised_files_all_buckets_$curdate.txt 

# create list of GWAS results to process through meta-analysis pipeline
mv .gwas_results_all_buckets_filtered.txt Bucket_scan/$curdate/gwas_results_for_munging.txt

# print details of available results to screen
print_tty " " "Summary of available GWAS results:"
phenolist=`head -n 1 Bucket_scan/$curdate/study_contribution_table_$curdate.csv | awk 'BEGIN{FS=","}{for(i=3;i<=NF;i++){if($i!~/sampleinfo/){split($i,a,"_");print a[2]}}}' | sort | uniq`
for pheno in $phenolist
do
   studystring=`awk -v pheno=$pheno 'BEGIN{FS=","}NR==1{for(i=3;i<=NF;i++){split($i,a,"_");if(a[2]==pheno){c[i]=a[1]}};next}{for(i=3;i<=NF;i++){if((i in c)&&length($i)>0){if(!($1 in s)){s[$1]=c[i]}else{s[$1]=s[$1]","c[i]}}}}END{str="";for(study in s){if(str==""){str=study" ("s[study]")"}else{str=str", "study" ("s[study]")"}};print str}' Bucket_scan/$curdate/study_contribution_table_$curdate.csv  | sort -k1`
   print_tty "${pheno}:    $studystring"
done

# print details of available sample summary data
print_tty " " "Summary of available sample summary data by data freeze:"
freezelist=`head -n 1 Bucket_scan/$curdate/study_contribution_table_$curdate.csv | awk 'BEGIN{FS=","}{for(i=3;i<=NF;i++){if($i~/sampleinfo/){split($i,a,"_");print a[1]}}}' | sed 's/F//g' | sort -n | awk '{print "F"$1}' | uniq`
for freeze in $freezelist
do
   studystring=`awk -v freeze=$freeze 'BEGIN{FS=","}NR==1{for(i=3;i<=NF;i++){if($i==freeze"_sampleinfo"){c[i]=freeze}};next}{for(i=3;i<=NF;i++){if((i in c)&&length($i)>0){if(!($1 in s)){s[$1]=c[i]}else{s[$1]=s[$1]","c[i]}}}}END{str="";for(study in s){if(str==""){str=study" ("s[study]")"}else{str=str", "study" ("s[study]")"}};print str}' Bucket_scan/$curdate/study_contribution_table_$curdate.csv  | sort -k1`
   print_tty "${freeze}:    $studystring"
done
echo

# print info messages
print_info "Output files created in folder \"Bucket_scan/$curdate/\"" "Please see \"study_contribution_table_$curdate.csv\" for case/control counts for identified results." " " "GWAS results to be used in the meta-analysis pipeline are listed in \"gwas_results_for_munging.txt\"." " " "A list of unrecognised files is stored in \"unrecognised_files_all_buckets_$curdate.txt\"." "Please check this list, as misnamed files will not be included in the study report and thus" "in the meta-analysis." " " "If you recognise GWAS results or sample summary files that should be included, please manually rename" "these files to match the specification in the Long COVID analysis plan and then rerun this script." " " "A full inventory can be found in \"full_inventory_all_buckets.txt\"." ""

################
### CLEAN UP ###
################

# move the full inventory file to the reports folder
mv .full_inventory_all_buckets.txt Bucket_scan/$curdate/full_inventory_all_buckets.txt

# remove hidden files 
rm -f .meta_files_all_buckets.txt .gwas_results_all_buckets.txt
