task clean_filter {

    String docker
    File sumstat_file
    String af_col
    Float min_af
    String info_col
    Float min_info
    String outfile = sub(basename(sumstat_file, ".gz"), "\\.bgz$", "") + ".munged.AF." + min_af + ".INFO." + min_info + ".gz"
    String dollar = "$"

    command <<<

        #TODO use pandas or something

        echo "COVID-19 HGI meta-analysis - clean and filter sumstats"
        echo "${sumstat_file}"
        echo ""

        catcmd="cat"
        if [[ ${sumstat_file} == *.gz ]] || [[ ${sumstat_file} == *.bgz ]]
        then
            catcmd="zcat"
        fi

        echo "`date` original number of variants"
        $catcmd ${sumstat_file} | tail -n+2 | wc -l

        chr_col=$($catcmd ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "CHR" | head -1 | cut -d ':' -f1)
        pos_col=$($catcmd ${sumstat_file} | head -1 | tr '\t ' '\n' | grep -nx "POS" | head -1 | cut -d ':' -f1)
        printf "`date` col CHR "${dollar}{chr_col}" col POS "${dollar}{pos_col}"\n"

        $catcmd ${sumstat_file} | awk ' \
          BEGIN{FS="\t| "; OFS="\t"}
          NR==1 {
              for (i=1;i<=NF;i++) { sub("^INFO$", "${info_col}", $i); sub("^Rsq$", "${info_col}", $i); sub("^CHR", "#CHR", $i); a[$i]=i; if ($i=="POS") pos=i }
              gsub("Pvalue", "p.value", $0);
              print $0
          } NR>1 {
              sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]);
              if ($a["#CHR"] ~ /^[0-9]+$/ && $a["p.value"] != 0 && $a["BETA"] < 1e6 && $a["BETA"] > -1e6 && $a["${af_col}"]>${min_af} && (1-$a["${af_col}"])>${min_af} && $a["${info_col}"]>${min_info}) {
                  printf $1
                  for (i=2; i<=NF; i++) {
                      if (i==pos) {
                          printf "\t%d", $i
                      } else {
                          printf "\t"$i
                      }
                  }
                  printf "\n"
              }
          }' | \
        sort -k$chr_col,${dollar}{chr_col}g -k$pos_col,${dollar}{pos_col}g -u | \
        bgzip > ${outfile}
        tabix -S 1 -s $chr_col -b $pos_col -e $pos_col ${outfile}

        echo "`date` new number of variants"
        gunzip -c ${outfile} | tail -n+2 | wc -l
        echo "`date` headers"
        gunzip -c ${outfile} | head -1 | tr '\t' '\n'

        gunzip -c ${outfile} | tail -n+2 | cut -f $chr_col | uniq > chr.tmp
        echo "`date` $(wc -l chr.tmp | cut -d' ' -f1) chromosomes"
        cat chr.tmp

        echo "`date` unique number of fields"
        gunzip -c ${outfile} | awk 'BEGIN{FS="\t"} {print NF}' | sort -u > n.tmp
        cat n.tmp
        if [ $(wc -l n.tmp | cut -d' ' -f1) != 1 ]; then echo "file not square"; exit 1; fi
        if [ $(wc -l chr.tmp | cut -d' ' -f1) -lt 22 ]; then echo "less than 22 chromosomes"; exit 1; fi

        echo "`date` done"

    >>>

    output {
        File out = outfile
        File tbi = outfile + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

task lift {

    String docker
    File sumstat_file
    File tbi_file = sumstat_file + ".tbi"
    String base = basename(sumstat_file)
    File b37_ref
    File b38_ref
    String dollar = "$"

    command <<<

        echo "COVID-19 HGI meta-analysis - lift over sumstats if needed"
        echo "${sumstat_file}"
        echo "${b37_ref}"
        echo "${b38_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${tbi_file} ${base}.tbi

        tabix -R ${b37_ref} ${base} | wc -l > b37.txt
        tabix -R ${b38_ref} ${base} | wc -l > b38.txt

        echo "`date` `cat b37.txt` chr 21 variants build 37"
        echo "`date` `cat b38.txt` chr 21 variants build 38"

        if ((`cat b37.txt` == 0 && `cat b38.txt` == 0)); then
            echo "`date` no chr 21 variants found in either build, quitting"
            exit 1
        fi

        if ((`cat b37.txt` > `cat b38.txt`)); then
            echo "`date` lifting to build 38"
            time /META_ANALYSIS/scripts/lift.py -chr "#CHR" -pos POS -ref Allele1 -alt Allele2 \
            -chain_file /liftover/hg19ToHg38.over.chain.gz -tmp_path /cromwell_root/ \
            ${base} > ${base}.lift.out 2> ${base}.lift.err
            gunzip -c ${base}.lifted.gz | \
            cut -f2- | awk '
            BEGIN { FS=OFS="\t" }
            NR==1 { for (i=1;i<=NF;i++) a[$i]=i; print $0 }
            NR>1 {
                temp=$a["#CHR"]; $a["#CHR"]=$a["anew_chr"]; $a["anew_chr"]=temp; temp=$a["POS"]; $a["POS"]=$a["anew_pos"]; $a["anew_pos"]=temp;
                sub("^0", "", $a["#CHR"]); sub("^chr", "", $a["#CHR"]); sub("^X", "23", $a["#CHR"]);
                if ($a["#CHR"] ~ /^[0-9]+$/) {
                    print $0
                }
            }' | bgzip > ${base}
        else
            echo "`date` presumably already in build 38"
        fi

    >>>

    output {
        File out = base
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

task harmonize {

    String docker
    File sumstat_file
    String base = basename(sumstat_file)
    File gnomad_ref
    String gnomad_ref_base = basename(gnomad_ref)
    Int n
    String options

    command <<<

        echo "COVID-19 HGI meta-analysis - harmonize sumstats to reference"
        echo "${sumstat_file}"
        echo "${gnomad_ref}"
        echo ""

        mv ${sumstat_file} ${base}
        mv ${gnomad_ref} ${gnomad_ref_base}

        echo "`date` harmonizing stats with gnomAD"
        python3 /META_ANALYSIS/scripts/harmonize.py ${base} ${gnomad_ref_base} ${n} ${options}\
        | bgzip > ${base}.${gnomad_ref_base} && \
        tabix -S 1 -s 1 -b 2 -e 2 ${base}.${gnomad_ref_base} && \
        echo "`date` done"

    >>>

    output {
        File out = base + "." + gnomad_ref_base
        File out_tbi = base + "." + gnomad_ref_base + ".tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

task flip_mahalanobis {

    String docker
    File sumstat_file
    Int sd
    String outfile = basename(sumstat_file, ".gz") + ".mahalanobis_filt_" + sd + "sd"

    command <<<

        python3 - <<EOF

        from scipy.spatial.distance import mahalanobis
        import scipy as sp
        import pandas as pd
        import numpy as np
        import re

        data = pd.read_csv('${sumstat_file}', sep='\t')
        gnomad_col = list(filter(re.compile('AF_gnomad_v3_b38_ref_').match, data.columns))[0]

        # flip frequencies and betas
        data['flipped'] = np.where(((data['AF_Allele2'] > 0.6) & (data[gnomad_col] < 0.4)) | ((data['AF_Allele2'] < 0.4) & (data[gnomad_col] > 0.6)), 1, 0)
        data['AF_Allele2'] = np.where(data['flipped'] == 1, 1-data['AF_Allele2'], data['AF_Allele2'])
        data['AF_fc'] = np.where(data['flipped'] == 1, data['AF_Allele2']/data[gnomad_col], data['AF_fc'])
        data['BETA'] = np.where(data['flipped'] == 1, -data['BETA'], data['BETA'])

        # calculate squared mahalanobis distance based on non-X variants
        try:
            afs = data[data['#CHR'] != 23][['AF_Allele2', gnomad_col]]
            inv_cov = sp.linalg.inv(afs.cov().values)
            mean = afs.mean().values
            data['mahalanobis2'] = afs.apply(lambda row: mahalanobis(row, mean, inv_cov) ** 2, axis=1)
            # filter data based on mahalanobis except X chr
            mean = np.mean(data['mahalanobis2'])
            sd = np.std(data['mahalanobis2'])
            data = data[(data['#CHR'] == 23) | ((data['mahalanobis2'] > mean - ${sd} * sd) & (data['mahalanobis2'] < mean + ${sd} * sd))]
        except np.linalg.LinAlgError: # singular matrix if AF is fixed
            data['mahalanobis2'] = np.nan

        data.to_csv('${outfile}', index=False, na_rep="NA", sep='\t')
        EOF

        bgzip -@4 ${outfile}
        tabix -S 1 -s 1 -b 2 -e 2 ${outfile}.gz

    >>>

    output {
        File out = outfile + ".gz"
        File out_tbi = outfile + ".gz.tbi"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        # TODO memory by file size
        memory: "20 GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

task plot {

    File sumstat_file
    String base = basename(sumstat_file)
    String prefix = sub(base, "\\.txt.*", "")
    Int loglog_ylim
    String pop
    String docker

    command <<<

        gunzip -c ${sumstat_file} | awk '
        BEGIN {FS=OFS="\t"}
        NR==1 {for(i=1;i<=NF;i++) { a[$i]=i; if ($i=="#CHR" || $i=="POS" || $i=="p.value" || $i~"AF_") b[i]=1}}
        {sep=""; for(i=1;i<=NF;i++) if (b[i]==1) { printf sep""$i; sep="\t"} printf "\n"}
        ' | bgzip > ${base} && \

        Rscript - <<EOF
        require(ggplot2)
        require(data.table)
        options(bitmapType='cairo')
        data <- fread("${base}")
        png("${base}_AF.png", width=1000, height=1000, units="px")
        p <- ggplot(data, aes_string(x="AF_Allele2", y="AF_gnomad_v3_b38_ref_${pop}")) +
          geom_point(alpha=0.1) +
          xlab("AF_${base}") +
          theme_minimal(base_size=18)
        print(p)
        dev.off()
        EOF

        qqplot.R --file ${base} --bp_col "POS" --chrcol "#CHR" --pval_col "p.value" --loglog_ylim ${loglog_ylim}
        [[ ! "${base}" =~ "HOSTAGE" && ! "${base}" =~ "Ancestry" ]] && qq_plot.R --input=${base} --prefix=${prefix} --af=AF_Allele2 --pvalue=p.value || echo "qq_plot.R not run"
    >>>

    output {
        Array[File] pngs = glob("*.png")
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: 10*ceil(size(sumstat_file, "G")) + " GB"
        disks: "local-disk 200 HDD"
        zones: "us-east1-d"
        preemptible: 2
        noAddress: true
    }
}

workflow munge_sumstats {

    File sumstats_loc
    Array[Array[String]] sumstat_files = read_tsv(sumstats_loc)
    String gnomad_ref_template

    scatter (sumstat_file in sumstat_files) {
        call clean_filter {
            input: sumstat_file=sumstat_file[0]
        }
        call lift {
            input: sumstat_file=clean_filter.out
        }
        call harmonize {
            input: sumstat_file=lift.out, gnomad_ref=sub(gnomad_ref_template, "POP", sumstat_file[1]), n=sumstat_file[2]
        }
        #call flip_mahalanobis {
        #    input: sumstat_file=harmonize.out
        #}
        #call plot {
        #    input: sumstat_file=flip_mahalanobis.out, pop=sumstat_file[1]
        #}
        call plot {
            input: sumstat_file=harmonize.out, pop=sumstat_file[1]
        }
    }
}
