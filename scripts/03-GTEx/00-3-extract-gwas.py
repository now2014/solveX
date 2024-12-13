#!/usr/bin/env python3

import gzip
import os

def extract_gwas(gz_file, snp_file, out_file):
    with open(snp_file, 'r') as f:
        snps = set([line.strip() for line in f])
    
    with gzip.open(gz_file, 'rt') as fhI, \
        open(out_file, 'w') as fhO:
        fhO.write(fhI.readline())

        for i, line in enumerate(fhI, 1):
            snp = line.split('\t', 2)[1]
            if snp in snps:
                fhO.write(line)
            if i % 1000000 == 0:
                print(f'{i/1000000} M lines processed.')


def main():
    # change working directory
    os.chdir('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')
    gz_file = '../GTEx_Analysis_v7_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'
    snp_file = 'data/selected-snps.tsv'
    out_file = 'data/gwas.tsv'
    extract_gwas(gz_file, snp_file, out_file)

if __name__ == '__main__':
    main()