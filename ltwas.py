#!/usr/bin/env python
'''
Local TWAS

lTWAS

Created on 2021-1-16

@author: Yiliang
'''

import argparse, os.path, sys
import pandas as pd
import numpy as np
from prep import prep
from calculate import calculate


try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('SUPERGNOVA requires pandas version > 0.15.2')


pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


# returns whether the parent directory of path exists
def parent_dir_exists(path):
    return os.path.exists(os.path.abspath(os.path.join(path, os.pardir)))

def pipeline(args):
    pd.options.mode.chained_assignment = None

    # Sanity check args
    if not parent_dir_exists(args.out):
        raise ValueError('--out flag points to an invalid path.')

    print('Preparing files for analysis...')
    gwas_snps, N1, N2 = prep(args.bfile, args.chr, args.start, args.end, args.sumstats1, args.sumstats2, args.N1, args.N2)
    print('Calculating local TWAS...')
    out = calculate(args.bfile, gwas_snps, N1, N2, args.h1, args.h2)
    out.to_csv(args.out, sep=' ', na_rep='NA', index=False)


parser = argparse.ArgumentParser()

parser.add_argument('sumstats1',
    help='The first sumstats file.')
parser.add_argument('sumstats2',
    help='The second sumstats file.')

parser.add_argument('--bfile', required=True, type=str,
    help='Prefix for Plink .bed/.bim/.fam file.')
parser.add_argument('--chr', required=True, type=str,
    help='Chromosome of the region')
parser.add_argument('--start', required=True, type=str,
    help='Start position of the genomic region')
parser.add_argument('--end', required=True, type=str,
    help='End position of the genomic region')
parser.add_argument('--N1', type=int,
    help='N of the sumstats1 file. If not provided, this value will be inferred '
    'from the sumstats1 arg.')
parser.add_argument('--N2', type=int,
    help='N of the sumstats2 file. If not provided, this value will be inferred '
    'from the sumstats2 arg.')
parser.add_argument('--h1', type=int,
    help='Local heritability of the first trait.')
parser.add_argument('--h2', type=int,
    help='Local heritability of the second trait.')

parser.add_argument('--out', required=True, type=str,
    help='Location to output results.')

if __name__ == '__main__':
    pipeline(parser.parse_args())

