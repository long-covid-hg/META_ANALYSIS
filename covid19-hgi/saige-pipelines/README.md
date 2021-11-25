# saige-pipelines

## Running GWAS

How to run SAIGE GWAS with Cromwell  
This in an example scenario creating new phenotypes in R5 and running those

1. Create a covariate/phenotype file that contains the covariates, the PCs and the phenotypes for Analysis_A1-C2 (cases 1, controls 0, everyone else NA), and upload the file. (see script 'merge_pheno_PCs.R in my home on qc-main')
2. Create a text file with the phenotypes one per line, e.g.  
    my_phenos.txt
    ```
    Analysis_B2
    Analysis_C2
    ```
    and upload the file.
3. Change the paths to the different file in `saige.json`  
4. Connect to Cromwell server using cromwell interact:
`./CromwellInteract-master/cromwell_interact.py connect cromwell-covid`
5. Submit workflow (use a sensible label so you can retrieve the job id in ./CromwellInteract-master/workflows.log):  
`python3 ./CromwellInteract-master/cromwell_interact.py submit --wdl [..]/covid19-hgi/saige-pipelines/wdl/saige.wdl --inputs [..]/covid19-hgi/saige-pipelines/wdl/[..].json --deps [..]/covid19-hgi/saige-pipelines/wdl/saige_sub.zip --label `    
6. Check the workflow status/failed_jobs:  
`./CromwellInteract-master/cromwell_interact.py metadata <JOB_ID> --failed_jobs`  
`./CromwellInteract-master/cromwell_interact.py metadata <JOB_ID> --s`
7. Logs and results go under:  
`gs://cromwell-covid/saige/JOB_ID`, plots `gs://cromwell-covid/saige/JOB_ID/call-test_combine/shard-*/**/*.png`, summary stats and tabix indexes `gs://cromwell-covid/saige/JOB_ID/call-test_combine/shard-*/**/*.gz*`
