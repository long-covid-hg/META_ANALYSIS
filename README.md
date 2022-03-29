# COVID-19 HGI META ANALYSIS workflow for long COVID

[google cloud platform for long-covid](https://console.cloud.google.com/home/dashboard?project=long-covid-hg)

## 1. ssh connect 

`gcloud compute ssh --ssh-flag="-X" long-covid-hg-cromwell --zone us-central1-b`

## 2. META-ANALYSIS
Mainly adopted from the [main HGI meta-analyses](https://github.com/covid19-hg/META_ANALYSIS)

### 2.1. set up google VM cromwell (see [cromwell_google_setup](https://github.com/long-covid-hg/cromwell_google_setup) repository)

### 2.2. pull docker image from main covid-19-hg project

example
```
gcloud docker -- pull gcr.io/covid-19-hg/plots:0.2
docker tag gcr.io/covid-19-hg/plots:0.2 gcr.io/long-covid-hg/plots:0.2
gcloud docker -- push gcr.io/long-covid-hg/plots:0.2
```

### 2.3 format checking format.wdl

```{bash}
lsof -ti:4999 | xargs kill -9
covid19-hgi/CromwellInteract-master/cromwell_interact.py connect long-covid-hg-cromwell
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 connect long-covid-hg-cromwell
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/format.wdl --inputs wdl/format.json
```

### 2.4. running munge.wdl
```
lsof -ti:4999 | xargs kill -9
covid19-hgi/CromwellInteract-master/cromwell_interact.py connect long-covid-hg-cromwell 

covid19-hgi/CromwellInteract-master/cromwell_interact.py \
submit --wdl wdl/munge.wdl \
--inputs wdl/munge.json

covid19-hgi/CromwellInteract-master/cromwell_interact.py metadata \
${JOBID} --failed_jobs
```

### 2.5 make configuration files for meta.wdl

```
gsutil ls gs://long-covid-hg-cromwell/munge/757864ee-79dc-4be8-9f1c-84a280f19d44/call-harmonize/shard-*/*.gz >> data/DF2/config_meta_F2.tsv

vi data/DF2/config_meta_F2.tsv

for pheno in NQ1.3 NQ2.3 WQ1.3 WQ2.3 N1.3 N2.3 W1.3 W2.3
do
python3 scripts/makejson.py --input data/DF2/config_meta_F2.tsv \
--output data/DF2/$pheno.json --pheno $pheno
done

python3 scripts/makesumstats.py --input data/DF2/config_meta_F2.tsv \
--output data/DF2/step3_sumstats_loc.txt

for pheno in NQ1.3 NQ2.3 WQ1.3 WQ2.3 N1.3 N2.3 W1.3 W2.3
do
echo "gs://long-covid-hg-cromwell/20220301/conf/$pheno.json" >> data/DF2/step3_pheno_conf.txt
done

gsutil cp data/DF2/* gs://long-covid-hg-cromwell/20220301/conf/
```

### 2.6 step3 run meta.wdl

```
python3 covid19-hgi/CromwellInteract-master/cromwell_interact.py --port 4999 submit --wdl wdl/meta.wdl --inputs wdl/meta.json

```
