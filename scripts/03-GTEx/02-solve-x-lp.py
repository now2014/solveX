#!/usr/bin/env python3
#SBATCH -J sHDL
#SBATCH -n 48
#SBATCH -w c002
#SBATCH --nodes=1

from scipy import stats
import gzip
import os
import sys
import numpy as np
from scipy.optimize import linprog
import time
import multiprocessing as mp

ncpus = 1
nthreads = 1
os.environ["OMP_NUM_THREADS"] = '%d'%ncpus
os.environ["OPENBLAS_NUM_THREADS"] = '%d'%ncpus
os.environ["MKL_NUM_THREADS"] = '%d'%ncpus
os.environ["VECLIB_MAXIMUM_THREADS"] = '%d'%ncpus
os.environ["NUMEXPR_NUM_THREADS"] = '%d'%nthreads

def solve_x(beta, Y, err=0.001, disp=False, maxiter=10000):

    options = {"disp": disp, "maxiter":maxiter}
    N = Y.shape[1]
    idx = np.isfinite(beta)

    M = sum(idx)
    if M == 0:
        return(np.ones(N) * np.nan)
    c = np.ones(N)

    # integrality = [1] * M

    A = Y[idx, :]
    beta = beta[idx]

    A_ub, b_ub = get_param(A, beta, err)
    res_linprog = linprog(c, A_ub = A_ub, b_ub = b_ub, bounds = (0, 2), options = options)
    status = res_linprog.status
    x = res_linprog.x

    # status:
    # 0 : Optimization proceeding nominally.
    # 1 : Iteration limit reached.
    # 2 : Problem appears to be infeasible.
    # 3 : Problem appears to be unbounded.
    # 4 : Numerical difficulties encountered.

    while x is None and err < 1:
        err = err + 0.001
        # print(err, status, x[:5])
        A_ub, b_ub = get_param(A, beta, err)
        res_linprog = linprog(c, A_ub = A_ub, b_ub = b_ub, bounds = (0, 2), options = options)
        status = res_linprog.status
        x = res_linprog.x

    if x is None:
        x = np.ones(N) * np.nan
    else:
        x = np.abs(np.rint(x))

    return(err, x)

def get_param(A, beta, err):
    M = len(beta)
    bu, bl = (1 + err) * beta, (1 - err) * beta
    b_ub = np.concatenate([(1 + err) * beta, (1 - err) * beta])
    
    pos_idx = b_ub > np.concatenate([beta, beta])
    pos_idx = pos_idx * 2 - 1
    b_ub = pos_idx * b_ub
    A_ub = np.concatenate([A, A]).T.dot(np.diag(pos_idx)).T
    return(A_ub, b_ub)

def read_b_Y(b_file='data/b.tsv.gz', Y_file='data/Y.tsv.gz'):
    dgene_Y = {}
    with gzip.open(Y_file, 'rt') as fhI:
        UIDs = fhI.readline().strip().split('\t')[1:]
        for line in fhI:
            gene, *Y = line.strip().split('\t')
            Y = np.array(Y, dtype=float)
            dgene_Y[gene] = Y
    
    with gzip.open(b_file, 'rt') as fhI:
        genes = fhI.readline().strip().split('\t')[1:]
        for line in fhI:
            line = line.replace('NA', 'nan')
            SNP, *b = line.strip().split('\t')
            b = np.array(b, dtype=float)
            idx = np.isfinite(b)
            if sum(idx) == 0:
                continue
            b = b[idx]
            kept_genes = [genes[i] for i in range(len(genes)) if idx[i]]
            kept_Y = np.array([dgene_Y[gene] for gene in kept_genes])
            yield SNP, b, kept_Y

def solve_genome(b_file='data/b.tsv.gz', Y_file='data/Y.tsv.gz', out_dir='data/solved-X'):
    out_tsv = 'X' + os.path.basename(b_file)[1:]
    out_tsv = os.path.join(out_dir, out_tsv)
    with open(out_tsv, 'wt') as fhO:
        for SNP, b, Y in read_b_Y(b_file, Y_file):
            err, x = solve_x(b, Y)
            # convert x to int
            x = np.rint(x).astype(int)
            fhO.write('%s\t%s\t%s\n'%(SNP, err, '\t'.join(map(str, x))))

def main():
    os.chdir('/share/home/lanao/projects/06-solveX/GTEx/v7/04-solve')
    out_dir = 'data/solved-X'
    os.makedirs(out_dir, exist_ok=True)
    
    Y_file = 'data/Y.tsv.gz'
    b_files = [os.path.join('./data', f) for f in os.listdir('./data') if f.endswith('.tsv.gz') and f.startswith('b')]

    with mp.Pool(48) as pool:
        pool.starmap(solve_genome, [(b_file, Y_file, out_dir) for b_file in b_files])


if __name__ == '__main__':
    main()
