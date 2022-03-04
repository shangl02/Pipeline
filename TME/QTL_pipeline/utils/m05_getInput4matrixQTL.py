'''
Created by Shangzhong.Li@pfizer.com in 2022/02
'''

import pandas as pd
import os,glob
import gzip
import argparse

parser = argparse.ArgumentParser(description='prepare input for matrixQTL')

# parser.add_argument('-c','--cov',action='store',dest='cov',help='covariate file')
# parser.add_argument('-d','--decon',action='store',dest='decon',help='cell deconvolution file')
parser.add_argument('-p','--path',action='store',dest='path',help='path')

args = parser.parse_args()

# cova_fn = args.cov
# cyto_fn = args.decon
work_path = args.path


plink_path = f'{work_path}/plink'
matrix_path = f'{work_path}/matrixQTL'
os.makedirs(matrix_path, exist_ok=True)

fam_fn = f'{plink_path}/germ_newSp38.fam'
cova_fn = f'{work_path}/covariates.tsv'
cyto_fn = f'{work_path}/cyto.tsv'
final_cyto_fn = f'{matrix_path}/cyto.txt'
final_cov_fn = f'{matrix_path}/covariates.txt'

# 1. plink samples
fam_df = pd.read_csv(fam_fn, header=None, sep=' ')
gt_samples = fam_df[1].tolist()

# 2. cytoreason file
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
# get cyto samples
cyto_samples = cyto_df.index.tolist()

# 3. get covariate file
cov_df = pd.read_csv(cova_fn, sep='\t',header=0,index_col=0)
cov_df['Sex'] = cov_df['Sex'].map(lambda x: 0 if x=='M' else 1)
cov_df = cov_df
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
    import scipy.stats as st
    rank_cyto_df = (t_cyto_df.rank(axis=1)-0.5).divide(t_cyto_df.rank(axis=1).max(axis=1),axis='index')
    rank_cyto_df = rank_cyto_df.apply(st.norm.ppf)
    t_cyto_df = rank_cyto_df

t_cyto_df[common_samples].to_csv(final_cyto_fn, sep='\t')

# 5. update covariates
cov_df.replace({'.':'NA'},inplace=True)
cov_df = cov_df.loc[common_samples,:].T
cov_df.index.name = 'id'
cov_df.to_csv(final_cov_fn, sep='\t')


def gt_file4matrixQTL(gt_path,out_gt_fn,analyze_samples):
    '''
    gt_path: folder that have all traw files
    out_gt_fn: out genotype file for matrixQTL
    '''
    in_gt_fns = sorted(glob.glob(gt_path + '/*38.traw.gz'))
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

out_gt_fn = f'{matrix_path}/SNP.txt'
gt_file4matrixQTL(plink_path, out_gt_fn, common_samples)