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
/home/Analysis/x11_fix.sh
```

### 1.1 Format checking (format.wdl)

Check the summary statistic files submitted to the buckets by the contributing studies (run -> read the INFO printed and check the files created)
```
cd /home/Analysis/
./1_scan_buckets.sh -b bucketlist.txt
```

Check the summary statistics formats (columns) and create a tab-separated list of summary statistic files and their formats to include in the meta-analyses (META_ANALYSIS/data/DF2/step1_format.txt)

#...

Copy to the bucket, under a directory for this meta-analysis date [the date is hard-coded for now, but could be automated]
gsutil cp data/DF2/step1_format.txt gs://long-covid-hg-cromwell/20220331/conf/

Run this pipeline from META_ANALYSIS directory
```
cd /home/Analysis/META_ANALYSIS/
```

Connect to Cromwell 
(Did you mean zone [us-central1-b] for instance: [long-covid-hg-cromwell] (Y/n)? [Y])
```
lsof -ti:4999 | xargs kill -9
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell
```

Change the date in the format.json (line 2)
```
vi wdl/format.json
```

Run the format step
```
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/format.wdl --inputs wdl/format.json
```

Save the HEX job id you got after submitting format.wdl (e.g. jobid=02390119-ac77-467e-8bc6-824cb70666b0; save these also to yourself in case the connection is lost)
```
jobid={jobid}
```

For checking the status of the jobs, you can use e.g.
```
covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary | tail -n 3
covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary 
covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --failed
watch --interval=10 'covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 metadata $jobid --failed_jobs | tail -n 25'
```

### 1.2 Munging (munge.wdl)

Create a list of formatted files and ancestries (step2_munge.txt) [change date] [This step will be further automated]
```
gsutil ls gs://long-covid-hg-cromwell/format/$jobid/call-formatting/**/*.gz
vi data/DF2/step2_munge.txt
gsutil cp data/DF2/step2_munge.txt gs://long-covid-hg-cromwell/20220331/conf/
```

(Connect to Cromwell again if the connection is not running anymore)
```
lsof -ti:4999 | xargs kill -9
covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell 
```

Change the date in the munge.json (line 2)
```
vi wdl/munge.json
```

#Run the munging step
```
covid19-hgi/CromwellInteract-master/cromwell_interact.py submit --wdl wdl/munge.wdl --inputs wdl/munge.json
```

#Save the HEX job id you got after submitting munge.wdl 
```
jobid={jobid}
```

### 1.3 Meta-analysis (meta.wdl)

Make configuration files for meta.wdl
This manual step has been replaced by the script generating this list, see the next step
(Manually make a list of the summary stat files to meta-analyse) 
```
gsutil ls gs://long-covid-hg-cromwell/munge/$jobid/call-harmonize/shard-*/**/*.gz >> data/DF2/config_meta_F2.tsv
#vi data/DF2/config_meta_F2.tsv
```

Run a script generating a list of the munged summary stat files to meta-analyse (config_meta_F2.tsv) [$jobid is the HEX ID from the munge job, or if you have munged in several jobs, add each of those separated by spaces]

```
cd /home/Analysis/
/generate_makejson_input.sh {$jobid}

mv config_meta_F2.tsv META_ANALYSIS/data/DF2/

cd /home/Analysis/META_ANALYSIS/
```

Create a .json file for each meta-analysis phenotype
(Note that if you change this list of meta phenos, it has to be changed (and kept in the same order) also in scripts/makesumstats.py and the for-loop generating step3_pheno_conf.txt)
```
for pheno in W1.3 W2.3 WQ1.3 WQ2.3 N1.3 N2.3 NE1.3 NE2.3 NQ1.3 NQ2.3
do
python3 scripts/makejson.py --input data/DF2/config_meta_F2.tsv \
--output data/DF2/$pheno.json --pheno $pheno
done
```

Create a list of munged summary stat locations (step3_sumstats_loc.txt)
(First remove a possible old version of step3_sumstats_loc.txt so that the new does not get appended to it)
```
mv data/DF2/step3_sumstats_loc.txt data/DF2/step3_sumstats_loc.txt.bkp

python3 scripts/makesumstats.py --input data/DF2/config_meta_F2.tsv --output data/DF2/step3_sumstats_loc.txt
```

Create a list of meta-analysis phenotypes to analyse ($pheno.json) [change date]
(Note that this list should be the same (and in same order) as in the first for-loop creating the pheno.jsons and in the scripts/makesumstats.py)
(First remove a possible old version of step3_pheno_conf.txt so that the new does not get appended to it)
```
mv data/DF2/step3_pheno_conf.txt data/DF2/step3_pheno_conf.txt.bkp

for pheno in W1.3 W2.3 WQ1.3 WQ2.3 N1.3 N2.3 NE1.3 NE2.3 NQ1.3 NQ2.3
do
echo "gs://long-covid-hg-cromwell/20220331/conf/$pheno.json" >> data/DF2/step3_pheno_conf.txt
done
```

Copy the files from DF2 to the cromwell bucket [change date]
```
gsutil cp data/DF2/* gs://long-covid-hg-cromwell/20220331/conf/
```

Change the dates in the meta.json (lines 2 and 3)
```
vi wdl/meta.json
```

Run the meta-analysis step
```
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/meta.wdl --inputs wdl/meta.json
```

Save the HEX job id you got after submitting meta.wdl 
```
jobid={jobid}
```

#Check status
```
covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999  metadata $jobid --summary | tail -n 7
```

