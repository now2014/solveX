#!/usr/bin/env python3

import gzip
import os
from collections import Counter

def count_snp(gz_file='GTEx_Analysis_v7_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'):
    dsnp_cnt = Counter()
    with gzip.open(gz_file, 'rt') as f:
        f.readline()
        for i, line in enumerate(f, 1):
            snp = line.split('\t', 2)[1]
            dsnp_cnt[snp] += 1
            if i % 1000000 == 0:
                print(f'{i/1000000} M lines processed.')
    return(dsnp_cnt)

def main():
    out_dir = './data'
    gz_file = 'GTEx_Analysis_v7_eQTL_all_associations/Whole_Blood.allpairs.txt.gz'
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, 'snp-cnt.tsv')
    dsnp_cnt = count_snp(gz_file)
    with open(out_file, 'wt') as f:
        f.write('SNP\tn.genes\n')
        for snp, cnt in dsnp_cnt.items():
            f.write(f'{snp}\t{cnt}\n')

if __name__ == '__main__':
    main()