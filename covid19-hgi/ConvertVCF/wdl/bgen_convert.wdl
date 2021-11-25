workflow convert_bgen {

    Array[String] chrom_list 
    String docker
    String name

    scatter (chrom in chrom_list){
        call chrom_convert {
            input :
            chrom = chrom,
            docker = docker,
            name = name
        }
    }

    call sum {
        input:
        values = chrom_convert.file_size,
        docker = docker}

    #call merge_bgen{
    #    input:
    #    name=name,
	#    docker = docker,
    #    total_size  = ceil(sum.out),
    #    bgen_files = chrom_convert.out_bgen,
    #    bgen_sample = chrom_convert.out_bgen_sample
    #    }

    # output of the module
    output {
        #File merged_bgen = merge_bgen.merged_bgen
        #File merged_bgen_index = merge_bgen.merged_bgen_index
        #File bgen_samples = merge_bgen.bgen_samples
        Array[File] out_bgen                 = chrom_convert.out_bgen
        Array[File] out_bgen_index           = chrom_convert.out_bgen_index
        Array[File] out_chrom_convert_log    = chrom_convert.out_chrom_convert_log
        Array[Array[File]] bgen_chunks       = chrom_convert.bgen_chunks
    }   
}

task chrom_convert {

    String chrom
    String name
    String docker

    # get path to vcf and tabix file
    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)
    String name_chrom =  name + "_" + chrom
    File tbiFile = cFile + '.tbi'

    Int disk_factor
    # Specify variant or positions (optional)
    # String variants_root
    # File vcf_variants = sub(variants_root,'CHROM',chrom)

    # get chrom size and convert to disk size & runtime options
    Int chrom_size = ceil(size(cFile,"GB"))
    String disk_size =disk_factor * chrom_size
    
    #bgen conversion options 
    String bargs 

    # vcf spliting options
    String sep

    String annotation
    String annotate  = if annotation == "" then annotation else " --annotate "
     
    String missingness
    String set_missingness  = if missingness == "" then missingness else " --set-missingness " + missingness  

    String chunk 
    String chunk_size = if chunk == "" then chunk else "--chunk-size " + chunk  

    # CPU AND MEM
    Int tmp_cpu =  if chrom_size < 48 then chrom_size*2 else 95
    Int split_cpu = if tmp_cpu < 8 then 8 else tmp_cpu
    
    Int tmp_mem = if chrom_size < 64 then chrom_size else 64
    Int mem = if tmp_mem < 16 then 16 else tmp_mem
   
    Int chrom_int = if chrom =="X" then 23 else chrom
    Int cpu = if chunk =="" then split_cpu else 30 - chrom_int # use less cpus in case of chunks
    
    command <<<
    echo ${disk_size}  ${mem} ${cpu} && df -h 
    python3 /Scripts/annotate.py \
    --sep ${sep} \
    --cFile ${cFile} \
    --tbiFile ${tbiFile} \
    --oPath "/cromwell_root/" \
    --name ${name_chrom} \
    --split \
    ${set_missingness} \
    ${annotate} \
    ${chunk_size} \
    -b \
    --bargs ${bargs} \
    | tee  chrom_convert_${name_chrom}.log        
    df -h >> chrom_convert_${name_chrom}.log
    >>>

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem}" + " GB"
	preemptible: 0
    }

    output {    
        File out_bgen = "/cromwell_root/${name_chrom}/${name_chrom}.bgen"
        File out_bgen_sample = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.sample"
        File out_bgen_index = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.bgi"
        Array[File] bgen_chunks = glob("/cromwell_root/${name_chrom}/bgen/*")
        File out_chrom_convert_log  = "chrom_convert_${name_chrom}.log"        
        Int file_size = chrom_size
    }
}



task merge_bgen {

    Array[File] bgen_files 
    Array[File] bgen_sample
    String name

    Int total_size
    # Int total_size = read_int(sum_chrom_sizes)
    Int disk_factor
    Int disk_size = total_size*disk_factor

    String mem = if total_size < 64 then total_size else 64
     
    String docker
    command <<<
     	
        python3 /Scripts/merge_bgen.py \
     	--fileList ${write_lines(bgen_files)} \
        --oPath "/cromwell_root/" \
        --name ${name} | tee merge_bgen_${name}.log
        
        df -h >> merge_bgen_${name}.log

    >>>

    runtime {
        docker: "${docker}"
	cpu: 64
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory:  "${mem}" + " GB"
	preemptible: 0

    }

    output {    
        File merged_bgen  = "/cromwell_root/${name}.bgen"
        File merged_bgen_index  = "/cromwell_root/${name}.bgen.bgi"
        File bgen_samples  = bgen_sample[0]
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
    disks: "local-disk 10 HDD"
    cpu: 1
    maxRetries: 1
  }

}
