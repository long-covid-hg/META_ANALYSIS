#!/usr/bin/python

import json
import re
import argparse

def makesumstats(input_file,output_file):
    for pheno in ['NQ1.3', 'NQ2.3', 'WQ1.3', 'WQ2.3', 'N1.3', 'N2.3', 'W1.3', 'W2.3']:
        arr = []
        file = open(input_file, 'r')
        for line in file:
            if re.search(pheno, line.split('\t')[0]):
                arr.append(line.split('\t')[2])

        with open(output_file, 'a', encoding='utf-8') as output:
            output.write('\t'.join(arr) + '\n')

def main(args):
    makesumstats(args.input,args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input_filename",type=str, required=True)
    parser.add_argument("--output", help="output_sumstats",type=str, required=True)
    args = parser.parse_args()
    main(args)
