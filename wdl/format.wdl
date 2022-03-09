version 1.0

workflow format {

    input {
        File sumstats_loc
        Array[Array[String]] summary_stats = read_tsv(sumstats_loc)
    }

    scatter (i in range(length(summary_stats))) {
        call formatting {
            input: summary_stat = summary_stats[i][0], format_info = summary_stats[i][1]
        }
    }

    output {
        Array[File] formated_sumstats = formatting.out
    }
}

# Filter bad quality variants
task formatting {

    input {
        File summary_stat

        String docker
        String format_info
	String outfile = sub(basename(summary_stat, ".gz"), "\\.bgz$", "") + ".formatted.gz"
  }

    command <<<

        set -eux

        echo "GWAS formatting"
        echo "~{summary_stat}"
        echo ""

        catcmd="cat"
        if [[ ~{summary_stat} == *.gz ]] || [[ ~{summary_stat} == *.bgz ]]
        then
            catcmd="zcat"
        fi

        if [[ ~{format_info} == "SAIGE" ]]
        then
          $catcmd ~{summary_stat} | bgzip > ~{outfile}
        fi

        if [[ ~{format_info} == "EXCEED" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["CHR"],$a["BP"],$a["Allele1"],$a["Allele2"],$a["AF_Allele2"],$a["imputationInfo"],$a["BETA"],$a["SE"],$a["P"],$a["N"]}' | \
          bgzip -@4 > ~{outfile}
        fi

        if [[ ~{format_info} == "REGENIE" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["CHROM"],$a["GENPOS"],$a["ALLELE0"],$a["ALLELE1"],$a["A1FREQ"],$a["INFO"],$a["BETA"],$a["SE"],10^-$a["LOG10P"],$a["N"]}' | \
          bgzip -@4 > ~{outfile}
        fi


        if [[ ~{format_info} == "JAPAN" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["CHR"],$a["POS"],$a["Allele1"],$a["Allele2"],$a["AF_Allele2"],$a["imputationInfo"],$a["BETA"],$a["SE"],$a["p.value"],$a["N.Cases"]+$a["N.Controls"]}' | \
          bgzip -@4 > ~{outfile}
        fi

        if [[ ~{format_info} == "PMBB" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["#CHROM"],$a["POS"],$a["ALLELE0"],$a["ALLELE1"],$a["A1FREQ"],$a["INFO"],$a["BETA"],$a["SE"],10^-$a["LOG10P"],$a["N"]}' | \
          bgzip -@4 > ~{outfile}
        fi

        if [[ ~{format_info} == "US" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["CHR"],$a["POS"],$a["Allele1"],$a["Allele2"],$a["AF_Allele2"],$a["imputationInfo"],$a["Beta"],$a["SE"],$a["p.value"],$a["N"]}' | \
          bgzip -@4 > ~{outfile}
        fi

        if [[ ~{format_info} == "HELIX" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {FS=OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {split($22,b,";|="); print $a["Chr"],$a["Pos"],$a["Ref"],$a["Alt"],$a["AAF"],b[6],b[2],b[4],$a["Pval"],($a["Num_Cases"]+$a["Num_Controls"])}' | \
          awk '$8 != "NA" && $9 != "NA"' | \
          bgzip -@4 > ~{outfile}
        fi

        if [[ ~{format_info} == "GENOTEK" ]]
        then
          $catcmd ~{summary_stat} | awk 'BEGIN {OFS="\t"}
          NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
          NR >1 {print $a["CHR"],$a["POS"],$a["Allele1"],$a["Allele2"],$a["AF_Allele2"],$a["imputationInfo"],$a["Beta"],$a["SE"],$a["p.value"],$a["N.Cases"]+$a["N.Controls"]}' | \
          bgzip -@4 > ~{outfile}
        fi

        echo "`date` done"

    >>>

    output {
        File out = sub(basename(summary_stat, ".gz"), "\\.bgz$", "") + ".formatted.gz"
    }

    runtime {
        docker: "~{docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk " + 5*ceil(size(summary_stat, "G")) + " HDD"
        zones: "us-central1-b"
        preemptible: 2
        noAddress: true
    }
}

