#!/usr/bin/python

import json
import re
import argparse

def tsv2json(input_file,output_file, pheno):
    arr = []
    file = open(input_file, 'r')
    a = file.readline()
      
    # The first line consist of headings of the record 
    # so we will store it in an array and move to 
    # next line in input_file.
    titles = [t.strip() for t in a.split('\t')]
    nfields = len(titles)
    for line in file:
        if re.search(pheno, line.split('\t')[0]):
            d = {}
            for t, f in zip(titles[1:nfields], line.split('\t')[1:nfields]):
                d[t] = f.strip().replace('gs://','/cromwell_root/')
            arr.append(d)
          
        # we will append all the individual dictionaires into list 
        # and dump into file.
    arr_up = {}
    arr_up['meta'] = arr
    
    with open(output_file, 'w', encoding='utf-8') as output_file:
        output_file.write(json.dumps(arr_up, indent=4))

def main(args):
    tsv2json(args.input,args.output, args.pheno)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input_filename",type=str, required=True)
    parser.add_argument("--output", help="output_json",type=str, required=True)
    parser.add_argument("--pheno", help="phenotype",type=str, required=True)
    args = parser.parse_args()
    main(args)
