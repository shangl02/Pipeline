'''
Created by Shangzhong.Li@pfizer.com in 2022/02
'''

import pandas as pd
import os,glob
import gzip
import argparse
import numpy as np
import scipy.stats as st

parser = argparse.ArgumentParser(description='prepare input for matrixQTL')

parser.add_argument('-p','--pre',action='store',dest='prefix',help='prefix of files')

args = parser.parse_args()
pre = args.prefix



fam_fn = f'{pre}.germ_dedup.fam'
cova_fn = f'{pre}.covariates.tsv'
cyto_fn = f'{pre}.cyto.tsv'
final_cyto_fn = f'{pre}.final_cyto.txt'
final_cov_fn = f'{pre}.final_covariates.txt'

# 1. plink samples
fam_df = pd.read_csv(fam_fn, header=None, sep='\t')
gt_samples = fam_df[1].tolist()

# 2. cytoreason file
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
# get cyto samples
cyto_samples = cyto_df.index.tolist()

# 3. get covariate file
cov_df = pd.read_csv(cova_fn, sep='\t',header=0,index_col=0)
cov_samples = cov_df.index.tolist()

common_samples = list(set(gt_samples).intersection(cyto_samples).intersection(cov_samples))

# 4. update cyto
t_cyto_df = cyto_df.loc[common_samples,:].transpose()
t_cyto_df.index.name = 'cell'
t_cyto_df = t_cyto_df[t_cyto_df.isna().sum(axis=1) < t_cyto_df.shape[1] * 0.1]
mode = 'norm'
if mode == 'log':
    t_cyto_df = t_cyto_df.apply(np.log)
elif mode == 'norm':
    rank_cyto_df = (t_cyto_df.rank(axis=1)-0.5).divide(t_cyto_df.rank(axis=1).max(axis=1),axis='index')
    rank_cyto_df = rank_cyto_df.apply(st.norm.ppf)
    t_cyto_df = rank_cyto_df

t_cyto_df[common_samples].to_csv(final_cyto_fn, sep='\t')

# 5. update covariates
# cov_df.replace({'.':'NA'},inplace=True)
cov_df.fillna('NA',inplace=True)
cov_df = cov_df.loc[common_samples,:].T
cov_df.index.name = 'id'
cov_df.to_csv(final_cov_fn, sep='\t')


def gt_file4matrixQTL(out_gt_fn,analyze_samples):
    '''
    gt_path: folder that have all traw files
    out_gt_fn: out genotype file for matrixQTL
    '''
    in_gt_fns = sorted(glob.glob('*.traw.gz'))
    header = True
    if os.path.exists(out_gt_fn):
        os.remove(out_gt_fn)
    for in_gt in in_gt_fns:
        # prepare SNPfiles
        with gzip.open(in_gt,'rt') as in_f, open(out_gt_fn,'a') as out_f:
            head = in_f.readline().strip().split('\t')
            if header:
                out_f.write('\t'.join(['SNP'] + analyze_samples)+'\n')
                header = False
            indexes = [head.index(s) for s in analyze_samples]
            for line in in_f:
                item = line.strip().split('\t')
#                 gts = [float(item[i]) for i in indexes]
                out_f.write('\t'.join([item[1]]+[item[i] for i in indexes])+'\n')

out_gt_fn = f'{pre}.SNP.txt'
gt_file4matrixQTL(out_gt_fn, common_samples)