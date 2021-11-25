workflow merge_vcf {

    Array[String] chrom_list 
    String docker

    scatter (chrom in chrom_list){
        call chrom_subset {
            input :
            chrom = chrom,
            docker = docker
        }
    }

    # output of the module
    output {
        Array[File] out = chrom_subset.out
        Array[File] out_tbi = chrom_subset.out_tbi
    }   
}


task chrom_subset {

    String chrom
    String docker

    String chromPathA
    File cFileA = sub(chromPathA,"CHROM",chrom)

    String chromPathB
    File cFileB = sub(chromPathB,"CHROM",chrom)


    command <<<

        echo "`date` indexing ..."
        time bcftools index ${cFileA}
        time bcftools index ${cFileB}
        echo "`date` merging ..."
        time bcftools merge --threads 64 -Oz --merge id --info-rules AF:avg,MAF:avg,R2:min ${cFileA} ${cFileB} > chr${chrom}.dose.vcf.gz
        echo "`date` tabix ..."
        time tabix chr${chrom}.dose.vcf.gz
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
        File out = "chr${chrom}.dose.vcf.gz"
        File out_tbi = "chr${chrom}.dose.vcf.gz.tbi"
    }
}
