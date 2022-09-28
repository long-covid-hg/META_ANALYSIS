# COVID-19 HGI META ANALYSIS workflow for Long COVID

[google cloud platform for long-covid](https://console.cloud.google.com/home/dashboard?project=long-covid-hg)

For a new meta-analysis run on long-covid-hg-cromwell virtual machine, go to step 1. Meta-analysis

## 0. Connect to Google Virtual Machine (VM) and set up Cromwell

### 0.1 ssh connect 

`gcloud compute ssh --ssh-flag="-X" long-covid-hg-cromwell --zone us-central1-b`

### 0.2 Set up Google VM Cromwell (see [cromwell_google_setup](https://github.com/long-covid-hg/cromwell_google_setup) repository)

### 0.3 Pull Docker image from the main covid-19-hg project

example
```
gcloud docker -- pull gcr.io/covid-19-hg/plots:0.2
docker tag gcr.io/covid-19-hg/plots:0.2 gcr.io/long-covid-hg/plots:0.2
gcloud docker -- push gcr.io/long-covid-hg/plots:0.2
```

## 1. META-ANALYSIS
Mainly adopted from the [main HGI meta-analyses](https://github.com/covid19-hg/META_ANALYSIS)
Start a new meta-analysis run from this step (if the Cromwell has already been set up)

### 1.0 X11 forwarding setup

Run the whole pipeline as superuser (root)
```
sudo su
```

If you have issues with X11 forwarding, run first
```
/home/Analysis/META_ANALYSIS/scripts/x11_fix.sh
```

### 1.1 Format checking (format.wdl)

Run this pipeline from the META_ANALYSIS directory
```
cd /home/Analysis/META_ANALYSIS/
```

Check the summary statistic files submitted to the buckets by the contributing studies (run -> read the INFO printed and check the files created)
```
scripts/1_scan_buckets.sh -b scripts/bucketlist.txt
```
This will create in the /home/Analysis/Bucket_scan/ a directory called YYYY.MM.DD/ and a number of files. Check the "gwas_results_for_munging.txt" file to ensure that you are not missing any GWAS summary files - "unrecognised_files_all_buckets_YYYY.MM.DD.txt" will list files that don't match the required name format. If you see any GWAS summary files listed here, rename the files to match the LongCOVID filename specification and rerun this script, repeating the process until no more GWAS summary files are left unrecognised.

Check the summary statistics formats (columns) and create a tab-separated list of summary statistic files and their formats to include in the meta-analyses (step1_format.txt) 
In case some sum stat files have formats deviating from SAIGE / REGENIE formats, you can add their transformation in the scripts/format.wdl

Copy to the bucket, under a directory for this meta-analysis date [the date is hard-coded in the .json files so far]. First set the analysis date variable in YYYYMMDD format (e.g. AnalysisDate=20220430].
```
AnalysisDate=[YYYYMMDD]
DataFreeze=[DF3]
mkdir data/$DataFreeze
gsutil cp data/$DataFreeze/step1_format.txt gs://long-covid-hg-cromwell/$AnalysisDate/conf/
```

Connect to Cromwell 
(Did you mean zone [us-central1-b] for instance: [long-covid-hg-cromwell] (Y/n)? [Y])
```
lsof -ti:4999 | xargs kill -9
python3 CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell
```

Change the date in the format.json (line 2)
```
vi wdl/format.json
```

Run the format step
```
python3 CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/format.wdl --inputs wdl/format.json
```

Save the HEX job id you got after submitting format.wdl (e.g. jobid=02390119-ac77-467e-8bc6-824cb70666b0; save these also to yourself in case the connection is lost)
```
jobid={jobid}
```

For checking the status of the jobs, you can use e.g.
```
CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary | tail -n 3
CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary 
CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --failed
watch --interval=10 "CromwellInteract-master/cromwell_interact.py --port 4999 metadata $jobid --failed_jobs | tail -n 25"
```

### 1.2 Munging (munge.wdl)

First, copy the required scripts to the relevant bucket folder:
```
for sfile in harmonize.py meta_analysis.py qqplot.R
do
   gsutil cp scripts/$sfile gs://long-covid-hg-cromwell/$AnalysisDate/scripts/
done
```
Edit the dates in the munge.json file's options (on lines 2, 24, 28) to correspond to your AnalysisDate [YYYMMDD]
```
vi wdl/munge.json
```

Create a list of formatted files and ancestries (step2_munge.txt) 
(JOBID= HEX ID(s) of formatting job(s))
First update the current Data Freeze (DF) in the generate_munge_input.sh script output path (e.g. data/DF4/step2_munge.txt)
```
vi scripts/generate_munge_input.sh
scripts/generate_munge_input.sh $jobid
gsutil cp data/$DataFreeze/step2_munge.txt gs://long-covid-hg-cromwell/$AnalysisDate/conf/
```

(Connect to Cromwell again if the connection is not running anymore)
```
lsof -ti:4999 | xargs kill -9
CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell 
```

#Run the munging step
```
CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/munge.wdl --inputs wdl/munge.json
```

#Save the HEX job id you got after submitting munge.wdl 
```
jobid={jobid}
```

### 1.3 Meta-analysis (meta.wdl)

Use the same analysis date as with the fromat and munging steps
```
AnalysisDate={[YYYYMMDD]}
```

First, copy the required scripts to the bucket location for this run (if you haven't already in the previous munging step):
```
for sfile in harmonize.py meta_analysis.py qqplot.R
do
   gsutil cp scripts/$sfile gs://long-covid-hg-cromwell/$AnalysisDate/scripts/
done
```
Edit all the dates in the wdl/meta.json file's options (lines 2, 3, 5, 12) to correspond to your AnalysisDate [YYYMMDD]
```
vi wdl/meta.json
```

Make configuration files for meta.wdl

##OBS! This script is not prioritising the different phenotypes correctly for the meta-analyses - has to be fixed before using##
##Run a script generating a list of the munged summary stat files to meta-analyse (config_meta_F4.tsv) [$jobid is the HEX ID from the munge job, or if you have munged in several jobs, add each of those separated by spaces] 
[DataFreeze is still hard-coded in the scripts generate_makejson_input.sh and makejson.py - change if you're running other than DF4]##
```
##scripts/generate_makejson_input.sh ${jobid}##
```

Create a .json file for each meta-analysis phenotype
```
for pheno in `cut -f1 data/$DataFreeze/config_meta_F4.tsv | tail -n +2 | sort | uniq`
do
   python3 scripts/makejson.py --input data/$DataFreeze/config_meta_F4.tsv --output data/$DataFreeze/$pheno.json --pheno $pheno
done
```

Create a list of munged summary stat locations (step3_sumstats_loc.txt)
(First back up your old version of step3_sumstats_loc.txt if you want to keep it, as this command will replace it)
```
phenolist=`python3 scripts/makesumstats.py --input data/$DataFreeze/config_meta_F4.tsv --output data/$DataFreeze/step3_sumstats_loc.txt`
```

Create a list of meta-analysis phenotypes to analyse (step3_pheno_conf.txt) 

(Again, back up the old version of step3_pheno_conf.txt if you want to keep it)
```
{
   for pheno in $phenolist
   do
      echo "gs://long-covid-hg-cromwell/$AnalysisDate/conf/$pheno.json"
   done
} > data/$DataFreeze/step3_pheno_conf.txt
```

Copy the files from $DataFreeze to the cromwell bucket 
```
gsutil cp data/$DataFreeze/* gs://long-covid-hg-cromwell/$AnalysisDate/conf/
```


Run the meta-analysis step (with the correct analysis dates updated in the meta.json)
```
lsof -ti:4999 | xargs kill -9
CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell 

python3 CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/meta.wdl --inputs wdl/meta.json
```

Save the HEX job id you got after submitting meta.wdl 
```
jobid={jobid}
```

#Check status
```
CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary | tail -n 7
```

