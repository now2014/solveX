#!/usr/bin/env python3

import os
import gzip
import argparse
import multiprocessing as mp

def read_idx(idx_file='data/04-sample-identification/idx.random-SNPs-MAF-0.05.tsv'):
    dchrom_kept_idx = {}
    with open(idx_file, 'rt') as fhI:
        fhI.readline()
        for line in fhI:
            chrom, i = line.strip().split('\t', 2)[:2]
            dchrom_kept_idx.setdefault(chrom, []).append(int(i))
    return(dchrom_kept_idx)

def subset_vcf(chrom, data_dir, out_dir, kept_idx=[]):
    prefix = f'ALL.chr{chrom}.'
    vcf_file = [fn for fn in os.listdir(data_dir) if fn.endswith('.vcf.gz') \
                and fn.startswith(prefix)]
    if len(vcf_file) != 1:
        print(f'Error: {vcf_file}')
        return
    vcf_file = os.path.join(data_dir, vcf_file[0])

    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f'{chrom}.tsv')

    geno2cnt = {'0|0': '0', '0|1': '1', '1|0': '1', '1|1': '2'}

    kept_idx = sorted(set(kept_idx))
    with gzip.open(vcf_file, 'rt') as fhI, \
        open(out_file, 'wt') as fhO:
        ## skip header
        for line in fhI:
            if line.startswith('#CHROM'):
                break
        line = line[1:]
        fhO.write(f'i\t{line}') # write header

        if len(kept_idx) == 0:
            return
        ## extract SNPs
        for i, line in enumerate(fhI, 1):
            if i==kept_idx[0]:
                # convert geno to cnt
                line = line.strip().split('\t')
                genos = line[9:]
                # genos = [g.split(':')[0] for g in genos]
                genos = [geno2cnt.get(g, 'NA') for g in genos]
                line[9:] = genos
                line = '\t'.join(line)
                fhO.write(f'{i}\t{line}\n')

                # del i from kept_idx
                kept_idx = kept_idx[1:]
                if len(kept_idx) == 0:
                    break

def main():
    args = argparse.ArgumentParser(description='Subset VCF files')
    args.add_argument('-d', '--data', help='Directory with 1KG data.', default='data/1kg/phase3-GRCh37')
    args.add_argument('-o', '--out', help='Output directory.', default='data/04-sample-identification/vcf-random-SNPs-MAF-0.05')
    args.add_argument('--idx_file', help='Path to the idx file', default='data/04-sample-identification/idx.random-SNPs-MAF-0.05.tsv')
    args = args.parse_args()

    dchrom_kept_idx = read_idx(args.idx_file)
    pool = mp.Pool()
    for chrom, kept_idx in dchrom_kept_idx.items():
        pool.apply_async(subset_vcf, args=(chrom, args.data, args.out, kept_idx))
    pool.close()
    pool.join()

if __name__ == '__main__':
    main()