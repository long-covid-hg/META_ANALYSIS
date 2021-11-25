workflow convert_plink {


    Array[String] chrom_list 
    String name
    String docker
    Int disk_factor


    scatter (chrom in chrom_list){
       call chrom_convert {
	   input :
           docker = docker,
	   chrom = chrom,
           disk_factor = disk_factor
       }     
    }
     
    # gathers the results from the file_size scatter
    call sum {input: values = chrom_convert.file_size,docker = docker}

    call merge_plink {
        input:
	bed_files = chrom_convert.bed,
	bim_files = chrom_convert.bim,
	fam_files = chrom_convert.fam,
        docker = docker,
        name = name,
        disk_factor = disk_factor,
        total_size  = ceil(sum.out)
        }
        
}



task merge_plink {

    Array[File] bed_files 
    Array[File] bim_files
    Array[File] fam_files 
    
    String name
    String pargs

    Int total_size
    Int disk_factor
    Int mem
    String docker
    
    Int disk_size = total_size*disk_factor + 20
    Int plink_mem = mem*1000 - 2000

    command <<<
    cat ${write_lines(bed_files)} | sed -e 's/.bed//g' > merge_list.txt
    plink --merge-list merge_list.txt ${pargs} --keep-allele-order --memory ${plink_mem} --make-bed --out ${name}
    plink2 --bfile ${name} --keep-allele-order --freq --out ${name}
    >>>

    runtime {
        docker: "${docker}"
	cpu: 16
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory:  "${mem}" + " GB"
        preemptible: 0

    }

    output {    
       File out_plink_merged_bed  = "${name}.bed"
       File out_plink_merged_bim  = "${name}.bim"
       File out_plink_merged_fam  = "${name}.fam"
       File out_plink_merged_freq = "${name}.afreq"
    }
}

task chrom_convert {

    String chrom
    String docker 
    String pargs
    Int disk_factor
    Int mem
    Int cpu
    
    # get path to vcf
    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)
    Int chrom_size = ceil(size(cFile,"GB"))

    
    Int disk_size = disk_factor * chrom_size + 20
    Int plink_mem = mem*1000 - 2000

    
    command <<<
    plink2 --vcf ${cFile} \
    ${pargs} --vcf-half-call h \
    --memory ${plink_mem} \
    --make-bed --out ${chrom} 
    >>>

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
	memory:  "${mem}" + " GB"
        preemptible: 0

    }

    output {	
	File bed         = "${chrom}.bed"
	File bim         = "${chrom}.bim"
	File fam         = "${chrom}.fam"
        Int file_size = chrom_size
    }
}
task sum {

  Array[Float] values
  String docker

  command <<<
    python -c "print(${sep="+" values})"
  >>>

  output {
    Float out = read_float(stdout())
  }

  runtime {
    docker: "${docker}"
    memory: "1 GB"
    cpu: 1
    maxRetries: 1
  }
}
