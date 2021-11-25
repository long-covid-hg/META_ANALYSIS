workflow extracts_snps_vcf {

    Array[String] chrom_list 
    String docker

    scatter (chrom in chrom_list){
        call chrom_subset {
            input :
            chrom = chrom,
            docker = docker
        }
    }

    call gather {
        input :
        chromosomes = chrom_subset.out,
        docker = docker
    }
}


task chrom_subset {

    String chrom
    String docker

    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)
    File cTbi = sub(chromPath,"CHROM",chrom) + ".tbi"

    File variants_extract

    command <<<

        echo "`date` extract variants ..."
        time bcftools view -Oz \
        -i "ID = @${variants_extract}" \
        ${cFile} \
        > chr${chrom}.extract.vcf.gz
        echo "`date` tabix ..."
        time tabix chr${chrom}.extract.vcf.gz
        echo "`date` DONE"

    >>>

    runtime {
        docker: "${docker}"
        cpu: "8"
        disks: "local-disk 200 HDD"
        memory: "8 GB"
	    preemptible: 0
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }

    output {    
        File out = "chr${chrom}.extract.vcf.gz"
        File out_tbi = "chr${chrom}.extract.vcf.gz.tbi"
    }
}


task gather {

    String docker

    Array[File] chromosomes

    command <<<

        echo "${chromosomes}"

        echo "`date` merging ..."
        time bcftools concat -Oz ${chromosomes} > all_chroms.extract.vcf.gz
        echo "`date` tabix ..."
        time tabix all_chroms.extract.vcf.gz
        echo "`date` DONE"

    >>>

    runtime {
        docker: "${docker}"
        cpu: "8"
        disks: "local-disk 200 HDD"
        memory: "8 GB"
	    preemptible: 0
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }

    output {    
        File out = "all_chroms.extract.vcf.gz"
        File out_tbi = "all_chroms.extract.vcf.gz.tbi"
    }
}