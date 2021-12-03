# COVID-19 HGI META ANALYSIS workflow for long COVID freeze 1

[google cloud platform for long-covid](https://console.cloud.google.com/home/dashboard?project=long-covid-hg)

## 1. ssh connect 

`gcloud compute ssh --ssh-flag="-X" long-covid-hg-cromwell --zone us-central1-b`

## 2. META-ANALYSIS
Mainly adopted from the [main HGI meta-analyses](https://github.com/covid19-hg/META_ANALYSIS)

### 2.1. set up google VM cromwell (see cromwell_google_setup repository)

### 2.2. pull docker image from main covid-19-hg project

pulled images
`saige:0.36.3.2-2` `meta:1d50c` `plots:0.2`


example
```
gcloud docker -- pull gcr.io/covid-19-hg/plots:0.2
docker tag gcr.io/covid-19-hg/plots:0.2 gcr.io/long-covid-hg/plots:0.2
gcloud docker -- push gcr.io/long-covid-hg/plots:0.2
```

### 3.1. running munge_sumstats.wdl
```
lsof -ti:5000 | xargs kill -9
covid19-hgi/CromwellInteract-master/cromwell_interact.py connect long-covid-hg-cromwell 

covid19-hgi/CromwellInteract-master/cromwell_interact.py \
submit --wdl wdl/munge_sumstats.wdl \
--inputs wdl/munge_sumstats.json

covid19-hgi/CromwellInteract-master/cromwell_interact.py metadata \
${JOBID} --failed_jobs
```
