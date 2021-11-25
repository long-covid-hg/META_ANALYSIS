# covid19-hgi

## QC
1\. Run **QC** on the whole dataset using the `1preimp_qc.py` script. Here's an example below on how to run it:
```
hailctl dataproc submit hail 1preimp_qc.py \
--dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/ \
--basename COV_ILLUMINA_15102020.chr0.pos0.removed \
--phenofile gs://dsge-covid19-data/COV_ILLUMINA_15102020/cases_controls_phenotypes_ita_be_brazil_swe_ger_29_10_2020 \
--inputType plink [or hail]
```
 * **Input(s):** PLINK .bed, .bim, and .fam files OR a Hail MatrixTable 
 * **Output(s):** (1) One or Two (if you use plink files as input) Hail MatrixTables, use the one with suffix .qc.mt for steps 2 and after; (2) log file; (3) a tsv file with all the samples removed due to relatedness, if there are any; (4) a tsv file containing of samples that failed sex check, if there are any  
 * Some of the QC checks included in step 1 above are: (1) sample/ variant missingness rate; (2) MAF checks; (3) sample relatedness (IBD); (4) sample sex check  

2\. Run **ancestry PCA** using the `2ancestry_pca.py` script
```
hailctl dataproc submit hail 2ancestry_pca.py \
--intersect_ref --pca_project --overwrite \
--data_dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/qc_step1/ \
--data_basename COV_ILLUMINA_15102020.chr0.pos0.removed.qc \
--out_prefix gs://dsge-covid19-data/COV_ILLUMINA_15102020/
```
 * **Input(s):** qced Hail MatrixTable (from step1 above) 
 * **Output(s):** (1) Two PCA scores files, one reference (.txt.bgz) and one for the data (.tsv)
 
3\. Run **ancestry assignment** using the `3assign_pops.py` script. The script assigns population labels based on the results of PCA. Here is an example below:
 ```
hailctl dataproc submit hail 3assign_pops.py \
--ref_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/qc_step2/COV_ILLUMINA_15102020.chr0.pos0.removed.qc_data_scores.txt.bgz \
--data_scores gs://dsge-covid19-data/COV_ILLUMINA_15102020/qc_step2/data_COV_ILLUMINA_15102020.chr0.pos0.removed.qc_cases_controls_scores.tsv \
--dirname gs://dsge-covid19-data/COV_ILLUMINA_15102020/
--phenotype_file gs://dsge-covid19-data/phenotypes_main_2021_04_11_wo_mismatch.tsv
--cohort italy --cohort belgium --cohort egypt --cohort iran --cohort sweden --cohort germany
```
* **Input(s):** (1 --ref_scores)Reference and (2 --data_scores)data PCA scores files from step 2 above; (3 --cohort) name of cohorts you want to run HWE on (one cohort name each).
* **Output(s):** (1) Two scores files, one at p > 0.5 and another at p > 0.8, which can then be used for PCA plots.

4\. Run **Allele frequency check against gnomAD and HWE filtering** on each cohort separately using the `4allele_hwe_checks.py` script. Here is an example below:
 ```
hailctl dataproc submit hail 4allele_hwe_checks.py
--dirname gs://dsge-covid19-data/15042021/
--mt gs://dsge-covid19-data/15042021/cov_15042021.merged.qc.mt
--pcasamples gs://dsge-covid19-data/15042021/pca_samples_2021_04_12
--phenofile gs://dsge-covid19-data/phenotypes_main_2021_04_11_wo_mismatch.tsv
--cohort italy --cohort belgium --cohort egypt --cohort iran --cohort sweden --cohort germany
```
* **Input(s):** (1 --mt) QCed Hail MatrixTable (from step 1); (2 --pcasamples) file, without a header, containing list of samples from PCA to be used in HWE filtering (one sampleID per line); (3) --cohort name of cohorts you want to run HWE on (one cohort name each).
* **Output(s):** (1) HWE QCed Hail MatrixTable for each cohort specified in the cohorts file input; (2) a directory named after the cohort with log file and a png file showing the filtered SNPs during allele frequency check against gnomAD.


## Pre-imputation
https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/

The tools necessary for these steps are installed on a VM called **pre-imp**, in my home dir (/home/cordioli/). Add the following lines to your `~/.bashrc` file so that you can run the commands (e.g. plinks) without having to specify the whole path:

`alias plink='/home/cordioli/plink'`   
`alias bcftools='/home/cordioli/bcftools/bcftools'`   
`alias qctool='/home/cordioli/qctool/build/release/qctool_v2.0.7'`   

**Input:** .fam,.bed,.bim files for the Qc'ed data  
**Output:** VCFs (one per each chromosome) to submit to the imputation server

0. Rename chr names (from chrN to N, as in the ref panel):  
(if duplicates:`cut -f 2 <> | sort | uniq -d > dups`)
`plink --bfile <> --update-chr /home/cordioli/ucsc2ensembl.txt --exclude dups --make-bed --out <>`
1. Create frequency file:  
`plink --freq --bfile <input> --out <freq-file>` 
2. Execute check script:  
`perl /home/cordioli/HRC-1000G-check-bim.pl -b <bim file> -f <freq-file> -r /home/cordioli/PASS.Variantsbravo-dbsnp-all.tab.gz -n -h`   
3. The perl script ran in the previous step generates a bash file with a set of plink commands to update or remove SNPs, split the bfile and create a VCF per each chromosome.
4. Sort VCFs  
`for i in {1..23}; do bcftools sort <PREFIX>-chr${i}.vcf -Oz -o <PREFIX>-chr${i}.vcf.gz & done`
5.  _"If your input data is GRCh38/hg38 please ensure chromosomes are encoded with prefix 'chr' (e.g. chr20)"_ :  
`for i in {1..23}; do bcftools annotate -Oz --rename-chrs /home/cordioli/ensembl2ucsc.txt <PREFIX>-chr${i}.vcf.gz > <PREFIX>_chr${i}.vcf.gz & done`
6. Submit to TopMed server, selecting "rsq Filter 0.3" to have the output already filtered for imputation score.


## Convert imputed VCF to bgen
0. Instead of downloading the imputed data to the VM and then transfer to the bucket, it's easier to just mount the bucket as an external disk:    
`gcsfuse --only-dir path/to/dir bucket-name mount-point`
e.g. `gcsfuse --only-dir 20211030/germany/imputed dsge-covid19-data /home/cordioli/20211030/germany/imputed/`
1. Download (in parallel) results to the 'pre-imp' VM using the wget links provided: ~10min    
`readarray -t urls < wgets.txt`  
`for i in "${urls[@]}"; do ${i} & done`  
Unzip (in parallel): `for i in *.zip; do unzip -P 'pwd' ${i} & done`
2. Recode chrX using Hail:  
2.a. run the recode_chrX.py script  
`gcloud dataproc jobs submit pyspark recode_chrX.py --cluster hail --region europe-west1 -- --dirname gs://dsge-covid19-data/ita_14082020/data_imputed/ --vcfname chrX.dose.vcf.gz`  
2.b. the output will have extension .bgz, so we need to rename both the original one (to keep it) and the recoded one:
`mv chrX.dose.vcf.gz chrX.dose.ori.vcf.gz`
`mv chrX.dose.vcf.bgz chrX.dose.vcf.gz`
3. Create tabix for each file (in parallel): ~20min  
`for i in *.dose.vcf.gz; do tabix ${i} & done`
4. Run ConvertVCF pipeline (see readme in the ConvertVCF subfolder): ~1h
5. Copy bgen chunks to bucket:   
`gsutil -m cp gs://cromwell-covid/convert_bgen/<JOB_ID>/call-chrom_convert/shard-*/glob-*/*.bgen* gs://dsge-covid19-data/<DEST_FOLDER>/bgen/`
6. Create file listing all bgen chunks:  
`gsutil ls gs://dsge-covid19-data/<DEST_FOLDER>/_*.*.bgen > bgen_pieces.txt`
and use script create_bgen_pieces_file.txt to create a file with 16 chunks per line.
7. Create file with samples IDs in the order they appear in the bgen:  
`gsutil cat <ANY OF THE ".bgen.sample" FILE> | tail -n +3 | awk '{print $1}' > bgen.samples.txt && gsutil mv bgen.samples.txt gs://dsge-covid19-data/<DEST-FOLDER>/conf/`


##  Create plink file for GRM
(should take ~2h)
1. Use only autosomes SNPs with info>0.95:  
`for i in {1..22}; do zcat chr${i}.info.gz | awk '$7>=0.95 {print $1}' >> variants_info_0.95; done`
2. Create a plink file for each chr (run in parallel to save time):  
`for i in {1..22};do plink --memory 3000 --vcf <VCF> --extract variants_info_0.95 --make-bed --const-fid --allow-extra-chr --out <OUT> & 
done`
3. Merge plink files and prune  
```
plink \
--bfile chr1 \
--merge-list merge_list.txt \
--make-bed \
--out chr1_22_info_0.95 && \
plink \
--memory 40000 \
--allow-extra-chr \
--snps-only \
--bfile chr1_22_info_0.95 \
--geno 0.03 \
--maf 0.01 \
--indep-pairwise 1000000 1000 0.2 \
--make-bed \
--out chr1_22_info_0.95.pruned && \
plink \
--memory 40000 \
--allow-extra-chr \
--bfile chr1_22_info_0.95.pruned \
--extract chr1_22_info_0.95.pruned.prune.in \
--make-bed \
--out <OUT>
```
