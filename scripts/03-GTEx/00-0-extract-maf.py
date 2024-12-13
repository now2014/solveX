#!/usr/bin/env python3

import gzip
import os

def read_maf(gz_file='GTEx_Analysis_v7_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'):
    dsnp_maf = dict()
    with gzip.open(gz_file, 'rt') as f:
        f.readline()
        for i, line in enumerate(f, 1):
            # snp = line.split('\t', 2)[1]
            snp, *_, maf = line.split('\t', 6)[1:6]
            dsnp_maf[snp] = maf
            if i % 1000000 == 0:
                print(f'{i/1000000} M lines processed.')
    return(dsnp_maf)

def main():
    out_dir = './data'
    gz_file = 'GTEx_Analysis_v7_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, 'snp-maf.tsv')
    dsnp_maf = read_maf(gz_file)
    with open(out_file, 'wt') as f:
        f.write('SNP\tMAF\n')
        for snp, maf in dsnp_maf.items():
            f.write(f'{snp}\t{maf}\n')

if __name__ == '__main__':
    main()