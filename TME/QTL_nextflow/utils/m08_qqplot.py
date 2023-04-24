'''
Created by Shangzhong.Li@pfizer.com on 2022/02/28
This file calculates QQ plot for the q values and check
beta distribution
'''

import argparse
import glob
import os
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='QQplot')


parser.add_argument('-p','--pre',action='store',dest='prefix',help='prefix of study')

args = parser.parse_args()
pre = args.prefix
qtl_fns = sorted(glob.glob(f'{pre}*.txt.gz'))

def GClambda(pval):
    pval = pval[~np.isnan(pval)]
    z = scipy.stats.norm.ppf(pval / 2)
    y = np.round(np.median(np.square(z)/0.456), 3)
    return y


def ppoints(n):
    # Generates the sequence of probability points 
    # (1:m - a)/(m + (1-a)-a) where m is either n, if length(n)==1, or length(n).
    a = 3./8. if n <= 10 else 1./2
    try:
        n = np.float(len(n))
    except TypeError:
        n = np.float(n)
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)

def qqplot(pval, figure, CI=0.95, topn = 10000):
    pval = pval[~np.isnan(pval)]
    n = len(pval)
    p = ppoints(n)
    pval_sort = np.sort(pval)

    # get lambda value
    lmda = GClambda(pval)
    ## qqline
    Qx = np.quantile(pval_sort, [0.25,0.75])
    Qz = np.array([0.25,0.75])
    b = (Qx[1] - Qx[0]) / (Qz[1] - Qz[0])
    a = Qx[0] - b * Qz[0]
    y = a + b * (-np.log10(p))
    
    # CI for topn values
    z = (1 - CI) / 2
    l = np.empty(topn)
    u = np.empty(topn)
    for i in range(1,topn+1):
        l[i-1] = -np.log10(scipy.stats.beta.ppf(1-z, i, n+1-i))
        u[i-1] = -np.log10(scipy.stats.beta.ppf(z, i, n+1-i))
    # plot
    plt.scatter(-np.log10(p), -np.log10(pval_sort))
    plt.xlabel('theoretical values')
    plt.ylabel('-log10(pval)')
    plt.plot(-np.log10(p), y, color='r')
    plt.plot(-np.log10(p)[:topn], l, color='b')
    plt.plot(-np.log10(p)[:topn], u, color='b')
    plt.text(0.1, 0.1, f'lambda={lmda}')
    plt.savefig(figure,dpi=300)

for qtl_fn in qtl_fns:
    figure = '.'.join(qtl_fn.split('/')[-1].split('.')[:-3]) + '.png'
    pval = pd.read_csv(qtl_fn, sep='\t', header=0, compression='gzip',usecols=[6])['p-value'].to_numpy()
    qqplot(pval, figure, CI=0.95, topn = 10000)