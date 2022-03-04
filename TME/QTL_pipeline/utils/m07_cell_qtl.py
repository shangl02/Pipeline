'''
Created by Shangzhong.Li@pfizer.com on 2021/01/29.
This file split the matrixQTL results into split results
'''

import argparse
import glob
import os
import pandas as pd


parser = argparse.ArgumentParser(description='split qtl results into cell types')

parser.add_argument('-p','--path',action='store',dest='path',help='work path')
# parser.add_argument('-i','--input',action='store',dest='input',help='input fn')
# parser.add_argument('-o','--outpath',action='store',dest='out_path',help='out_path')

args = parser.parse_args()
work_path = args.path
# qtl_fn = args.input
# out_path = args.out_path

qtl_fn = f'{work_path}/matrixQTL/qtl.txt.gz'
out_path = f'{work_path}/matrixQTL/cell_qtl'
os.makedirs(out_path,exist_ok=True)


for qtl_df in pd.read_csv(qtl_fn,sep='\t',header=0,compression='gzip',chunksize=1e6):
    cell_types = qtl_df['gene'].unique().tolist()
    for cell in cell_types:
        out_fn = f'{out_path}/{cell}.txt.gz'
        if not os.path.exists(out_fn):
            qtl_df.query('gene == @cell').to_csv(out_fn,sep='\t',index=False,compression='gzip')
        else:
            qtl_df.query('gene == @cell').to_csv(out_fn,sep='\t',index=False,
                compression='gzip',header=None,mode='a')