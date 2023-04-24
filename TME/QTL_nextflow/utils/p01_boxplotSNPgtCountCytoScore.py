'''
Created by Shangzhong.Li@pfizer.com in 2022.03.25.
This file plot the boxplot of snp affect on the cell type abundance.
Input parameters:
snp: needs to be in the format chr:pos
'''

import argparse
import os
import subprocess
import pandas as pd
import tabix
import gzip
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='plot genotye count vs gene expression for a snp')


parser.add_argument('-s','--snp',action='store',dest='snp',help='snp id in format chr:pos:ref:alt or rsid, depends on the id in the SNP.txt file')
parser.add_argument('-c','--cell',action='store',dest='cell',help='cell type')
parser.add_argument('-o','--out',action='store',dest='out',help='output file')
parser.add_argument('--na',action='store',dest='na',help='include cyto scores for people with no na')
parser.add_argument('--cyto',action='store',dest='cyto',help='deconvolusion file')
parser.add_argument('--gt',action='store',dest='gt',help='genotype file')

args = parser.parse_args()
snp = args.snp
cell = args.cell
cyto_fn = args.cyto
gt_fn = args.gt

out_fn = args.out
na = args.na

# gt_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WES/plink/germ_newSp38.traw.gz'
# cyto_fn = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WES/matrixQTL/cyto.txt'

def get_genotype(snp,gt_fn,samples):
    # 1. find the chromosome file with genotype
    chrom, pos= snp.split(':')[:2]
    # 2. get the header
    with gzip.open(gt_fn,'rt') as f:
        head = f.readline().strip().split('\t')
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
cyto_df = pd.read_csv(cyto_fn,sep='\t',header=0,index_col=0)
samples = cyto_df.columns.to_list()
cyto_df = cyto_df.T
# 2. get genotype
gt_df = get_genotype(snp,gt_fn,samples)
# 3. get cell type cyto scores
cell_gt_df = pd.merge(cyto_df[cell],gt_df,left_index=True,right_index=True)
cell_gt_df = cell_gt_df.sort_values('gt')
if na == 'no':
    cell_gt_df = cell_gt_df[cell_gt_df['gt'].isin(['0','1','2'])]
# 4. plot for a single plot
cell_gt_df = cell_gt_df.sort_values('gt',ascending=False)
ax = sns.boxplot(y=cell,x="gt",data=cell_gt_df,boxprops=dict(alpha=.5))
ax = sns.swarmplot(y=cell,x="gt",data=cell_gt_df)
ax.set_title(snp)
plt.savefig(out_fn)
