#!/usr/bin/python
from __future__ import division, print_function
import numpy as np
import pandas as pd
import numpy.linalg as linalg
from math import sqrt
import ld.ldscore as ld
import ld.parse as ps
from ldsc_thin import __filter_bim__
from scipy.stats import norm
from collections import OrderedDict


def nearest_Corr(input_mat):
    d, v = linalg.eigh(input_mat)
    A = (v * np.maximum(d, 0)).dot(v.T)
    A = (A + A.T) / 2
    multiplier = 1 / np.sqrt(np.diag(A))
    A = A * multiplier
    A = (A.T * multiplier).T
    return A


def calculate(bfile, gwas_snps, N1, N2, h1, h2):
    m = len(gwas_snps)

    snp_file, snp_obj = bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = bfile+'.bed', ld.PlinkBEDFile

    array_snps = snp_obj(snp_file)
    keep_snps = __filter_bim__(gwas_snps, array_snps)

    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    keep_indivs = None

    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    max_dist = 1
    block_left = ld.getBlockLefts(coords, max_dist)

    blockLD = geno_array.ldCorrVarBlocks(block_left)
    local_LD = nearest_Corr(blockLD)

    d, v = linalg.eigh(local_LD)
    order = d.argsort()[::-1]
    d = d[order]
    v = v[:,order]
    
    sub_d = d[d>0]
    sub_v = v[:,d>0]

    tz1 = np.dot(sub_v.T, gwas_snps['Z_x'])
    tz2 = np.dot(sub_v.T, gwas_snps['Z_y'])
    y = tz1 * tz2

    Z_x = gwas_snps['Z_x']
    Z_y = gwas_snps['Z_y']

    
    threshold = 1
    cur_d = sub_d[sub_d>threshold]
    cur_y = y[sub_d>threshold]
    cur_dsq = cur_d ** 2
    denominator = (wh1 * cur_d / m0 + 1 / n1) * (wh2 * cur_d / m0 + 1 / n2)
    cur_v1 = np.sum(cur_dsq / denominator)
    cur_v2 = np.sum(cur_y / sqrt(n1 * n2) / denominator)
    cur_v3 = np.sum(cur_y ** 2 / (n1 * n2) / (denominator * cur_dsq))

    emp_var = [(cur_v3 - (cur_v2 ** 2) / cur_v1) / (cur_v1 * (len(cur_d) - 1))]
    theo_var = [1 / cur_v1]

    for K in range(len(cur_d), len(sub_d)):
        eig = sub_d[K]
        tmp_y = y[K]
        cur_v1 += eig ** 2 / ((wh1 * eig / m0 + 1 / n1) * (wh2 * eig / m0 + 1 / n2))
        cur_v2 += tmp_y / sqrt(n1 * n2) / ((wh1 * eig / m0 + 1 / n1) * (wh2 * eig / m0 + 1 / n2))
        cur_v3 += tmp_y ** 2 / (n1 * n2) / ((wh1 * eig ** 2 / m0 + eig / n1) * (wh2 * eig ** 2 / m0 + eig / n2))
        emp_var.append((cur_v3 - (cur_v2 ** 2) / cur_v1) / (cur_v1 * K))
        theo_var.append(1 / cur_v1)
    
    max_emp_theo = np.maximum(emp_var, theo_var)
    min_idx = np.argmin(max_emp_theo)

    y = y[:(len(cur_d)+min_idx-1)]
    sub_d = sub_d[:(len(cur_d)+min_idx-1)]
    sub_dsq = sub_d ** 2

    var_rho = m0 ** 2 * min(max_emp_theo)
    q = (wh1 * sub_d / m0 + 1 / n1) * (wh2 * sub_d / m0 + 1 / n2)
    v4 = np.sum(sub_d/q)/np.sum(sub_dsq/q)
    var_phencorr = pheno_corr_var / (n1 * n2) * m0 ** 2 * v4 ** 2
    var_rho += var_phencorr

    se_rho = sqrt(var_rho)
    p_value = norm.sf(abs(Localrho / se_rho)) * 2

    if Localh1 < 0 or Localh2 < 0:
        corr = np.nan
    else:
        corr = Localrho / sqrt(Localh1 * Localh2)

    df = pd.DataFrame(OrderedDict({"chr":[CHR], "start":[START], "end":[END], "rho":[Localrho], "corr":[corr], "h2_1":[Localh1], "h2_2":[Localh2], "var":[var_rho], "p":[p_value], "m":[m0]}))

    return df

def _supergnova(bfile, partition, thread, gwas_snps, ld_scores, n1, n2, pheno_corr, pheno_corr_var):

    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]
    bps = np.array(array_snps.df['BP'])[geno_array.kept_snps]

    ## Calculating local genetic covariance
    
    results = []
    def collect_results(result):
        results.append(result)
    pool = multiprocessing.Pool(processes = thread)
    for i in range(blockN):
        pool.apply_async(calLocalCov, args=(i, tmp_partition, geno_array, coords, 
            bps, tmp_gwas_snps, tmp_ld_scores, n1, n2, pheno_corr, pheno_corr_var),
            callback=collect_results)
    pool.close()
    pool.join()
    df = pd.concat(results, ignore_index=True)
    #df = pd.DataFrame(results)
    #df.columns = ["chr", "start", "end", "rho", "corr", "h1", "h2", "var", "p", "m"]
    convert_dict = {"chr": int, "start": int, "end":int, "m":int}
    df = df.astype(convert_dict)
    return df
