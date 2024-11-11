#!/usr/bin/python

import json
import re
import argparse
import os

def getphenolist(input_file):
    file = open(input_file, 'r')
    phenolist = []
    next(file)
    for line in file:
        pheno=line.split('\t')[0]
        if pheno not in phenolist:
            phenolist.append(pheno)
    return phenolist


def makesumstats(input_file,phenolist,output_file):
    if os.path.exists(output_file):
        os.remove(output_file)
    for pheno in phenolist:
        arr = []
        print(pheno)
        file = open(input_file, 'r')
        for line in file:
            if re.search(pheno, line.split('\t')[0]):
                arr.append(line.split('\t')[2])

        with open(output_file, 'a', encoding='utf-8') as output:
            output.write('\t'.join(arr) + '\n')

def main(args):
    makesumstats(args.input,getphenolist(args.input),args.output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input_filename",type=str, required=True)
    parser.add_argument("--output", help="output_sumstats",type=str, required=True)
    args = parser.parse_args()
    main(args)
