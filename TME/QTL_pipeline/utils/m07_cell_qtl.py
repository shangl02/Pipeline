'''
Created by Shangzhong.Li@pfizer.com on 2021/01/29.
This file split the matrixQTL results into split results
'''

import argparse
import glob
import os
import pandas as pd
import glob
import gzip

parser = argparse.ArgumentParser(description='split qtl results into cell types')

parser.add_argument('-p','--path',action='store',dest='path',help='work path')

args = parser.parse_args()
work_path = args.path
bim_fn = f'{work_path}/plink/germ_newSp38.bim'


qtl_fn = f'{work_path}/matrixQTL/qtl.txt.gz'
out_path = f'{work_path}/matrixQTL/cell_qtl'
os.makedirs(out_path,exist_ok=True)

def add_snp_info(qtl_fn, bim_fn):
    bim_df = pd.read_csv(bim_fn, sep='\t',header=None,names=['chr','SNP','n','pos','ref','alt'])
    for qtl_df in pd.read_csv(qtl_fn,sep='\t',header=0,compression='gzip',chunksize=1e6):
        merge_df = pd.merge(qtl_df,bim_df[['SNP','chr','pos','ref','alt']],
                                on=['SNP'],how='left')
        merge_df['std'] = merge_df['beta'].div(merge_df['t-stat'])
        pre = ['chr','pos']
        columns = ['chr','pos'] + [c for c in merge_df.columns if c not in pre]
        for cell, df in merge_df.groupby('gene'):
            out_fn = f'{out_path}/{cell}.txt.gz'
            if not os.path.exists(out_fn):
                df[columns].to_csv(out_fn,sep='\t',index=False,compression='gzip')
            else:
                df[columns].to_csv(out_fn,sep='\t',index=False,
                    compression='gzip',header=None,mode='a')


def sort_bgzip(qtl_fn, out_path):
    temp_path = f'{out_path}/temp'
    os.makedirs(temp_path,exist_ok=True)
    # get header
    head_fn = f'{out_path}/head.txt'
    with gzip.open(qtl_fn,'rt') as fin, open(head_fn,'w') as fout:
        fout.write(fin.readline())
    # sort and bgzip
    sort_fn = f'{qtl_fn[:-7]}.sort.txt.gz'
    cmd = f'gunzip -c {qtl_fn} | tail -n +2 | \
           sort -k1,1 -k2,2n -T {temp_path} - | cat {head_fn} - | \
           bgzip > {sort_fn} && \
           tabix -p bed -f -S 1 -b 2 -e 2 {sort_fn}'
    os.system(cmd)
    os.remove(head_fn)
    os.remove(qtl_fn)
#----- 1. split files for each cell type -------------- 
# add_snp_info(qtl_fn, bim_fn)

# #------2. generate index for each cell type -----------
qtl_fns = [f for f in glob.glob(f'{out_path}/*.txt.gz') 
           if not f.endswith('.sort.txt.gz')]
for qtl_fn in qtl_fns:
    sort_bgzip(qtl_fn, out_path)