workflow subset_vcf {

    Array[String] chrom_list 
    String docker
    File samples

    scatter (chrom in chrom_list){
        call chrom_subset {
            input :
            chrom = chrom,
            docker = docker,
            samples = samples
        }
    }

    # output of the module
    output {
        Array[File] out = chrom_subset.out
    }   
}


task chrom_subset {

    String chrom
    String docker

    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)

    File samples

    command <<<

    echo "tabix"
    echo "${cFile}"
    tabix ${cFile}
    echo "subset"
    bcftools view --threads 16 -Oz --force-samples -S ${samples} ${cFile} > chr${chrom}.sub.vcf.gz

    >>>

    runtime {
        docker: "${docker}"
        cpu: "8"
        disks: "local-disk 50 HDD"
        memory: "16 GB"
	    preemptible: 0
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }

    output {    
        File out = "chr${chrom}.sub.vcf.gz"
    }
}
