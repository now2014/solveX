#!/usr/bin/env python3
#SBATCH -J extract1kg
#SBATCH -n 11
#SBATCH -w c002
#SBATCH --nodes=1

import os
import gzip
import argparse
import multiprocessing as mp

def append_to_file(out_file, lines, mode='a'):
    out_dir = os.path.dirname(out_file)
    if out_dir != '':
        os.makedirs(out_dir, exist_ok=True)
    with open(out_file, mode) as f:
        for line in lines:
            f.write(line)

def read_phase3_populations(data_dir):
    pop_file = os.path.join(data_dir, 'integrated_call_samples_v3.20130502.ALL.panel')
    related_file = os.path.join(data_dir, '20140625_related_individuals.txt')
    
    skipped_samples = set()
    with open(related_file, 'rt') as fhI:
        fhI.readline()
        for line in fhI:
            skipped_samples.add(line.split('\t', 1)[0])

    sample2pop = {}
    with open(pop_file, 'rt') as fhI:
        fhI.readline()
        for line in fhI:
            # sample, pop = line.split('\t', 3)[:2]
            sample, _, pop = line.split('\t', 4)[:3]
            if sample in skipped_samples:
                continue
            sample2pop[sample] = pop
    return(sample2pop)

def calc_af(genos, pop2i, pops):
    cnts0 = [g.count('0') for g in genos]
    cnts1 = [g.count('1') for g in genos]

    out_str = []
    for pop in pops:
        ixs = pop2i.get(pop, [])
        c0 = sum(cnts0[ix] for ix in ixs)
        c1 = sum(cnts1[ix] for ix in ixs)
        # out_str.append(f'{c0}:{c1}')
        try:
            af = c1 / (c0 + c1)
        except ZeroDivisionError:
            af = 0
        out_str.append(f'{af:.4f}')
    out_str = '\t'.join(out_str)
    return(out_str)

def extract_snps(chrom, data_dir, out_dir, sample2pop):
    prefix = f'ALL.chr{chrom}.'
    vcf_file = [fn for fn in os.listdir(data_dir) if fn.endswith('.vcf.gz') \
                and fn.startswith(prefix)]
    if len(vcf_file) != 1:
        print(f'Error: {vcf_file}')
        return
    vcf_file = os.path.join(data_dir, vcf_file[0])

    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f'{chrom}.tsv')

    pops = sorted(list(set(sample2pop.values())))
    pop2i = {pop:list() for pop in pops}

    with gzip.open(vcf_file, 'rt') as fhI, \
        open(out_file, 'wt') as fhO:
        ## skip header
        for line in fhI:
            if line.startswith('#CHROM'):
                break

        ## map samples to populations
        samples = line.strip().split('\t')[9:]
        for i, sample in enumerate(samples):
            if sample in sample2pop:
                pop = sample2pop[sample]
                pop2i[pop].append(i)
        fhO.write('\t'.join(['POS', 'REF', 'ALT'] + pops) + '\n')

        ## extract SNPs
        for line in fhI:
            pos = line.split('\t', 2)[1]
            contents = line.strip().split('\t')
            ref, alt = contents[3:5]
            af = calc_af(contents[9:], pop2i, pops)
            fhO.write('\t'.join([pos, ref, alt, af]) + '\n')

def main():
    args = argparse.ArgumentParser('Extract SNPs from 1KG.')
    args.add_argument('-d', '--data', help='Directory with 1KG data.', default='data/1kg/phase3-GRCh37')
    args.add_argument('-o', '--out', help='Output directory.', default='data/1kg-extracted')
    args = args.parse_args()

    chroms = [str(chrom) for chrom in range(1, 23)]
    sample2pop = read_phase3_populations(args.data)

    with mp.Pool(11) as pool:
        for chrom in chroms:
            pool.apply_async(extract_snps, (chrom, args.data, args.out, sample2pop))
        pool.close()
        pool.join()

if __name__ == '__main__':
    main()