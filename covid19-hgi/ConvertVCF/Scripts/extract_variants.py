from utils import make_sure_path_exists,mapcount,tmp_bash,basic_iterator,pretty_print
import os,subprocess,shlex

def do_variants(args):
   
    extract_variants(args)
    position_split(args)


def extract_variants(args):
    """
    Creates final list of positions for splitting.
    First variants are extracted (if not provided).
    Then variants are subset (if required).
    Then unique positions are extracted from variant ID
    """
    # GET LIST OF VARIANTS IN VCF IF NOT PROVIDED
    if not args.vcf_variants or mapcount(args.vcf_variants) ==0:
        args.force = True
        print('extracting variants from vcf.')
        args.vcf_variants = os.path.join(args.variant_path, args.name +"_variants_vcf.txt")
        if not os.path.isfile(args.vcf_variants) or mapcount(args.vcf_variants) == 0:           
            cmd = f"zcat {args.cFile}  | grep -v '^#' |  gawk -v OFS='\t' '{{print $3 }}'  > {args.vcf_variants}"
            tmp_bash(cmd)

    print(f'{mapcount(args.vcf_variants)} variants in the vcf file')
    
    # GET SHARED VARIANTS IF SUBSETTING IF REQUESTED. In this way it doesn't matter if splitting or not, the final list of positions will be equally valid (i.e. shorter in case of subset provided).
    
    if args.variant_file:
        print(f'{mapcount(args.variant_file)} variants in the variant file')
        final_ids = os.path.join(args.variant_path, args.name +"_variants_ids.txt")
        if not os.path.isfile(final_ids) or args.force :
            args.force = True
            cmd = f"""comm -12 <(sort {args.vcf_variants}) <(sort {args.variant_file})  > {final_ids}"""
            tmp_bash(cmd)
        print(f'{mapcount(final_ids)} variants are shared among the two lists')
        args.vcf_variants = final_ids # update final list of ids to be the shared one!

    # args.vcf_variants contains the final list of variants that need to be returned

    # we extract the positions based on the final list of variants
    args.final_positions = os.path.join(args.variant_path, args.name +"_positions_vcf.txt")       
    if not os.path.isfile(args.final_positions) or args.force:
        args.force = True
        cmd = f"""cat {args.vcf_variants} |  awk -F '{args.sep}' '{{print $1"\\t"$2  }}' | uniq > {args.final_positions}"""
        tmp_bash(cmd)
    print(f'{mapcount(args.final_positions)} unique positions to be used.')    



def position_split(args):
    '''
    Given the final position file, it splits it into chunks using regions based on whether it's split or whole file.
    '''
    
    args.variant_split_path = os.path.join(args.variant_path, 'chrom_split/')
    make_sure_path_exists(args.variant_split_path)
    
    out_root = os.path.join(args.variant_split_path, args.name)
    
    # here i keep track of the alternative chrom variants
    chunk_size = get_chunk_size(args)    

    if args.test:
        chunk_size = min(50,chunk_size)
        
    fid = 0   #file count
    count = 1  #line count
    start_pos = 1
    f = open(f"{out_root}_{fid}_positions.txt", 'w')

    for chrom,position in basic_iterator(args.final_positions) :
        # in the variant scenario, we want chrom_pos_pos for subsetting
        if args.variant_file:
            f.write(f"{chrom}\t{position}\t{position}\n")

        if count%chunk_size ==0:
            if args.split:
                final_pos = int(position) -1
                f.write(f"{chrom}\t{start_pos}\t{final_pos}\n")
            f.close()
            fid += 1 #increase file count
            f = open(f"{out_root}_{fid}_positions.txt", 'w')
            start_pos = int(position)
        count +=1

    if args.split:
        final_pos = int(position)
        f.write(f"{chrom}\t{start_pos}\t{final_pos}\n")    
    f.close()
    
    args.chunks = fid +1
    if args.test:
        args.chunks = args.cpus

    print(f"chunks : {args.chunks}")



def get_chunk_size(args):
    file_lines = mapcount(args.final_positions)

    if not args.chunk_size:
        chunks = args.chunks       

    else:
        # add one chunk to spread evenly the "extra" variants
        chunks = 1 + (file_lines//(args.chunk_size))
        
    chunk_size = int(file_lines/int(chunks)) + 1

    print(f"chunk size : {chunk_size}")
    return chunk_size
