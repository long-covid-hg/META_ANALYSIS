from utils import tmp_bash,make_sure_path_exists,file_exists,valid_string,pretty_print,pad,natural_sort,get_filepaths,mapcount
import os,multiprocessing,argparse,subprocess,shlex
from collections import defaultdict as dd
from itertools import product
from tempfile import NamedTemporaryFile

mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
mem_mib = mem_bytes/(1024.**2) 
plink_mem = max(mem_mib - 4000,int(0.8*mem_mib)) #leave some memory available

from extract_variants import do_variants


################################
#--------VCF SPLITTING --------#
################################

def split_vcf(args,proc):
    '''
    Splits the original chrom file in chunks based on the corresponding chunk region
    '''
    pos_file =  os.path.join(args.variant_split_path, f'{args.name}_{proc}_positions.txt')
    vcf_file = args.vcf_file_path.replace('CHUNK',str(proc))

    if not os.path.isfile(vcf_file) or args.force:
        args.force = True
        
        sample =f' -S {args.samples}' if args.samples else ''
        print(f'generating {vcf_file}')
        basic_cmd = f"bcftools view {args.cFile} -R {pos_file} -T {pos_file}  {sample} {args.vargs}"
        
        if args.variant_file and mapcount(args.variant_file) > 0:
            basic_cmd += f" --i 'ID=@{args.final_ids}' "

        #write add annotation to basic command
        if args.annotate:
            basic_cmd +=  " -Ou | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' "
            
        if args.set_missingness:
            gp = str(args.set_missingness)
            cmd = f""" {basic_cmd} -Ou | bcftools +setGT  -Oz -o  {vcf_file} -- -i 'FORMAT/GP[:0]< {gp} & FORMAT/GP[:1]< {gp} & FORMAT/GP[:2]< {gp}' -n '.' -t q """
            if proc ==0:
                print(cmd)
            tmp_bash(cmd)
            return

        # in the output from minimac4, chrX is represented as haploid genotypes for males, however, SAIGE can only
        # handle diploid genotypes. The command here produces unphased diploid genotypes.But since the haploid genotypes
        # are in diploid format as REF / REF or ALT / ALT, we can simply set the phase for those alleles with a
        # simple sed replacement.
        if args.fix_ploidy_x and (args.chr_x_id in vcf_file):
            basic_cmd += " -Ou | bcftools +fixploidy  -Ov -- -f 2 | sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | bcftools view "

        cmd = f"{basic_cmd}  -Oz -o { vcf_file}"
        if proc ==0:
            print(cmd)
        tmp_bash(cmd)
    else:
        print(f'{vcf_file}  already generated')
        
def vcf_wrapper(args):
    split_vcf(*args)

def vcf_multiproc(args):
    '''
    Function that calls all the operations to be done in parallel.
    '''

    args.vcf_file_path = os.path.join(args.vcf_path,args.name + '.CHUNK.vcf.gz')
    # the number of parallel processes is the smallest of the two
    pool = multiprocessing.Pool(args.cpus)
    paramList = list(product([args],range(args.chunks)))
    pool.map(vcf_wrapper,paramList)
    pool.close()

    print('vcf splitting done')
#############################################
#---------BGEN CONVERSION AND MERGE---------#
#############################################

'''
Pipeline for converting each annotate vcf.gz chunk in bgen format
'''

def bgen_convert(args,proc):
    '''
    Function that converts the chunk to bgen
    '''
    vcf_file = args.vcf_file_path.replace('CHUNK',str(proc))
    bgen_file =  args.bgen_file_path.replace('CHUNK',str(proc))    
    log_file = os.path.join(args.log_path,str(proc) + '_bgen.log')

    if not os.path.isfile(bgen_file) or args.force:
        args.force = True
        print(f'generating {bgen_file}')
        cmd = f'qctool -g {vcf_file} {args.bargs} -og {bgen_file} -os {bgen_file}.sample'
        if proc ==0: print(cmd)            
        subprocess.call(shlex.split(cmd),stderr = open(log_file,'w'))
        cmd = f' bgenix -g  {bgen_file} -clobber -index'
        subprocess.call(shlex.split(cmd))
      
    else:
        print(bgen_file + ' already generated')
    
    return None

def bgen_wrapper(args):
    bgen_convert(*args)
    
def bgen_multiproc(args):

    pool = multiprocessing.Pool(args.cpus)
    paramList = list(product([args],range(args.chunks)))
    pool.map(bgen_wrapper,paramList)
    pool.close()

def bgen_merge(args):
    '''
    bgen_merger
    '''
    args.bgen_file_path = os.path.join(args.bgen_path,args.name + '.CHUNK.bgen')
    merged_bgen = os.path.join(args.oPath, args.name + '.bgen')
    if not os.path.isfile(merged_bgen) or args.force:
        args.force = True
        # convert chunks
        bgen_multiproc(args)
        # remove vcf chunks if no other operation will take place
        if args.cleanup and not (args.plink or args.vcf):
            cleanup(args)
               
        #collect list of chunks and concatenate them
        bgen_files = [f for f in natural_sort(get_filepaths(args.bgen_path)) if f.endswith('.bgen')]
        cmd = 'cat-bgen -g ' + ' '.join(bgen_files)  +  ' -og ' + merged_bgen
        subprocess.call(shlex.split(cmd))
        #index final bgen
        cmd = f' bgenix -g  {merged_bgen} -clobber -index'
        subprocess.call(shlex.split(cmd))
        # create sample list
        sample_file = [f for f in natural_sort(get_filepaths(args.bgen_path)) if f.endswith('.sample')][0]
        cmd = f" cp {sample_file} {merged_bgen}.sample"
        subprocess.call(shlex.split(cmd))    
        #remove old bgen chunks
        if args.cleanup:
            for f in bgen_files:os.remove(f)
    else:
        print(f"bgen merged file already generated")
        subprocess.call(shlex.split(f' qctool -g  {merged_bgen}'))

        
##############################################
#---------PLINK CONVERSION AND MERGE---------#
##############################################

def plink_multiproc(args):

    
    pool = multiprocessing.Pool(args.cpus)
    paramList = list(product([args],range(args.chunks)))
    pool.map(plink_wrapper,paramList)
    pool.close()

def plink_wrapper(args):
    plink_convert(*args)
    
def plink_convert(args,proc):
    '''
    Convert vcf.gz chunk to plink
    '''
    vcf_file = args.vcf_file_path.replace('CHUNK',str(proc))    
    plink_root =  os.path.join(args.plink_path,args.name + "." + str(proc))
    memory = f" -memory {int( mem_mib / (args.cpus))}"
    threads = ' --threads 1'
    log_file = os.path.join(args.log_path,str(proc) + '_plink.log')
    
    if not os.path.isfile(plink_root + '.bed') or args.force:
        args.force = True
        print('generating ' + plink_root)
        cmd = f"plink2 --vcf {vcf_file}  {args.pconvargs} {memory} {threads}   --make-bed --out {plink_root}"
        if proc ==0 :
            print(cmd)
        subprocess.call(shlex.split(cmd),stdout = open(log_file,'w'))
    else:
        print(plink_root + ' already generated')

    return None

def plink_merge(args):
    '''
    Merges chunks together
    '''    
    merged_plink= os.path.join(args.oPath, args.name)
    if not os.path.isfile(merged_plink + '.bed') or args.force:
        args.force = True 
        #runk plink chunk conversion
        plink_multiproc(args)
        # cleanup chunks if needed
        if args.cleanup and not args.vcf:
            cleanup(args)
            
        #merge chunks
        merge_list = NamedTemporaryFile(delete=True)
        with open(merge_list.name,'w') as o:
            plink_files = [elem.split('.bed')[0] for elem in natural_sort(get_filepaths(args.plink_path)) if elem.endswith('.bed')]
            for f in plink_files:o.write(f + '\n')
        cmd =  f'plink --merge-list {merge_list.name} {args.pargs}  --memory {int(plink_mem)} --make-bed --out {merged_plink}'        
        subprocess.call(shlex.split(cmd))

        freq_file = merged_plink+'.afreq'
        cmd = f"plink2 --bfile {merged_plink} --freq --out {merged_plink}"
        subprocess.call(shlex.split(cmd))
        # remove old plink chunks
        if args.cleanup:
            for f in get_filepaths(args.plink_path):os.remove(f)
    
#############################
#---------VCF MERGE---------#
#############################
def bcf_merge(args):
    '''
    Merges vcf chunks
    '''

    # get list of chunk files in natural order (by chunk order) else the tabix fails
    vcf_files = natural_sort(get_filepaths(args.vcf_path))
    print(vcf_files)
    #write mergelist
    print('merging vcf.gz files...')
    vcf_file =  os.path.join(args.oPath,args.name + '.vcf.gz')
    
    merge_list  =  NamedTemporaryFile(delete=True)
    with open(merge_list.name,'w') as o:
        for f in vcf_files:o.write(f +'\n')
    cmd =  'bcftools concat -n -f ' + merge_list.name + '  -Oz -o ' + vcf_file
    subprocess.call(shlex.split(cmd))

    #index final vcf.gz
    print('creating index...')
    cmd = f'tabix -f {vcf_file}'
    subprocess.call(shlex.split(cmd))  
    print('done')
    
    if args.check_vcf is True:
        print('checking output..')
        #count number of variants of original chrom file
        cmd = 'bcftools index -s '
        files = [args.cFile,vcf_file]
        for f in files:
            subprocess.call(shlex.split(cmd +f))
    return None

def cleanup(args):
    files = get_filepaths(args.vcf_path)
    cmd = ' rm ' + ' '.join(files)
    subprocess.call(shlex.split(cmd))
    
def main(args):
  # REQUIRED PART: IT SPLITS THE  VCFS AND ANNOTATES THEM
     # convert cols to stringf
    args.log_path =  os.path.join(args.oPath , 'logs')
    make_sure_path_exists(args.log_path)
    
    args.variant_path = os.path.join(args.oPath ,'variants')
    make_sure_path_exists(args.variant_path)
    pretty_print('VARIANTS EXTRACTION')
    do_variants(args)

    args.vcf_path = os.path.join(args.oPath,'vcf')
    make_sure_path_exists(args.vcf_path)
    pretty_print('VCF SPLITTING')
    vcf_multiproc(args)

    #bgen conversion
    if args.bgen:
        pretty_print("BGEN CONVERSION")
        args.bgen_path = os.path.join(args.oPath,'bgen')
        make_sure_path_exists(args.bgen_path)
        bgen_merge(args)
        
    if args.plink:
        pretty_print("PLINK CONVERSION")
        args.plink_path = os.path.join(args.oPath,'plink')
        make_sure_path_exists(args.plink_path)
        plink_merge(args)

    #vcf merging
    if args.vcf:
        bcf_merge(args)
        if args.cleanup:
            cleanup(args)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description ="Annotation and conversion of VCF files")

    #Basic inputs
    parser.add_argument("--cFile", type= file_exists,
                        help="gs chrom file path to import",required = True)
    parser.add_argument("--tbiFile", type= file_exists,
                        help="gs tbi file path to import",required  = False )

    parser.add_argument("-o","--oPath",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument("--name", type= str,
                        help="chromosome number or output name in general ",required = True)
    parser.add_argument("--sep", type= str,
                        help="variant separator",default = "_")

    parser.add_argument("--fix_ploidy_x", action="store_true")
    parser.add_argument("--chr_x_id", type=str,
                        help="chrX code (X/23)", default="X")

    parser.add_argument("--chunks",type = int, help = "number of chunks", default =  None)
    parser.add_argument("--cpus",type = int, help = "number of cpus to use", default =  multiprocessing.cpu_count())

    parser.add_argument('--vcf-variants',type = file_exists,help = 'List of vcf variants',default = None)
    parser.add_argument('--samples',type= file_exists,help = 'Flag to choose what samples to keep')

    split_type = parser.add_mutually_exclusive_group(required = True)
    split_type.add_argument('--variant-file',type = file_exists,help = 'List of variants. If variatn file is empty, then all variants are processed.',default = None)
    split_type.add_argument('--split',action = 'store_true',help = 'Flag to split the vcf file in chunks')
    
    # output choices
    parser.add_argument('-v','--vcf',help = "output vcf",action = "store_true")
    parser.add_argument('-b','--bgen',help = "output bgen",action = "store_true")
    parser.add_argument('-p','--plink',help = "output plink",action = "store_true")

    # matching options for vcf/bgen/plink output
    parser.add_argument('--vargs',type = str,default = '',help = 'String with kwargs to pass to bcftools')

    parser.add_argument('--bargs',type = str,default = '',help = 'String with kwargs to pass to the conversion script, based on the desired output')
    parser.add_argument('--pargs',type = str,default = '',help = 'String with kwargs to pass to all plink calls, based on the desired output')
    parser.add_argument('--pconvargs',type = str,default = '',help = 'String with kwargs to pass to the conversion plink script, based on the desired output')

    
    parser.add_argument('--set-missingness',type = float,default = False,help = 'Annotates missingess based on GP') #only vcf flag at the moment
    parser.add_argument('--annotate',action = 'store_true',help = 'Flag to annotate variants')
    
    #optional tags
    parser.add_argument('--chunk-size',type = int,help = 'Force chunk size')
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force')
    parser.add_argument('--check-vcf',action = 'store_true',help = 'Check that output vcf files are identical')
    parser.add_argument('--test',action = 'store_true',help = 'Flag to run small chunks')
    parser.add_argument('--cleanup',action = 'store_true',help = 'Flag to delete vcf chunks at the end of the process')
    

    args = parser.parse_args()
    if args.oPath.endswith('/'):
        args.oPath = args.oPath[:-1]

    # creating outputpath
    args.oPath +=  '/' + args.name +'/'
    make_sure_path_exists(args.oPath)
    
    #adding padding to kwargs
    args.bargs =  pad(args.bargs)
    args.pargs =  pad(args.pargs)
    args.pconvargs = pad(args.pconvargs) +args.pargs

    if not args.chunks:
        args.chunks = args.cpus
    print(args)
    main(args)


#test purposes
#python3 annotate.py -vbp --cFile /mnt/disks/tera/vcf_files/R2_chr22.vcf.gz --oPath /mnt/disks/tera/annotate/ --chrom 22 --bargs ' -filetype vcf -bgen-bits 8 -bgen-compression zlib -vcf-genotype-field "GP" -bgen-permitted-input-rounding-error 0.005 -ofiletype "bgen_v1.2" ' --pconvargs ' --vcf-half-call h ' --pargs ' --biallelic-only strict --allow-extra-chr' --set-missingness 0.95 --cleanup --check-vcf --test --force
