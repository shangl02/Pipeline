'''
Created by Shangzhong.Li@pfizer.com in 2022.03.25.
This file plot the boxplot of snp affect on the cell type abundance.
Input parameters:
snp: needs to be in the format chr:pos
'''
import argparse
import os
import pandas as pd
import tabix
import gzip
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='plot genotye count for a snp and gene expression')


parser.add_argument('-s','--snp',action='store',dest='snp',help='snp id')
parser.add_argument('-g','--gid',action='store',dest='gene_id',help='ensembl gene id')
parser.add_argument('-o','--out',action='store',dest='out',help='output file')
parser.add_argument('-r','--tpm',action='store',dest='tpm',help='tpm file')
parser.add_argument('--na',action='store',dest='na',help='include cyto scores for people with no na')
parser.add_argument('--gt',action='store',dest='gt',help='genotype file')

args = parser.parse_args()
snp = args.snp
gid = args.gene_id
out_fn = args.out
tpm_fn = args.tpm
na = args.na
gt_fn = args.gt

def get_genotype(snp,gt_fn,samples):
    # 1. find the chromosome file with genotype
    chrom, pos= snp.split(':')[:2]
    # 2. get the header
    with gzip.open(gt_fn,'rt') as f:
        head = f.readline().strip().split('\t')
        samples = list(set(samples).intersection(head))
        iids = head[6:]
        indexes = [iids.index(s) for s in samples]
        tb = tabix.open(gt_fn)
        records = tb.query(chrom,int(pos), int(pos)+1)
        for record in records:
            record = record[6:]
            gt = [record[i] for i in indexes]
            gt_df = pd.DataFrame({'gt':gt},index=samples)
    return gt_df

# 1. get samples for the run
tpm_df = pd.read_csv(tpm_fn, sep='\t',header=0, index_col=0)
samples = tpm_df.columns.to_list()
tpm_df = tpm_df.T

# 2. get genotype
gt_df = get_genotype(snp, gt_fn, samples)
samples = gt_df.index.tolist()

# 3. get gene expression values
tpm_gt_df = pd.merge(tpm_df[gid],gt_df,left_index=True,right_index=True)

if na == 'no':
    tpm_gt_df = tpm_gt_df[tpm_gt_df['gt'].isin(['0','1','2'])]

# 4. plot for a single plot
tpm_gt_df = tpm_gt_df.sort_values('gt',ascending=False)
ax = sns.boxplot(y=gid,x="gt",data=tpm_gt_df,boxprops=dict(alpha=.5))
ax = sns.swarmplot(y=gid,x="gt",data=tpm_gt_df)
ax.set_title(snp)
ax.set_ylabel(f'{gid} (tpm)')
plt.savefig(out_fn)