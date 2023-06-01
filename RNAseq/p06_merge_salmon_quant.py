'''
Created on 2022/01/21 by Shangzhong.Li@pfizer.com
This merge salmon quant information
'''
import glob, os
from natsort import natsorted
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='merge quant into a table')

parser.add_argument('-i','--input',action='store',dest='input',help='input path')
parser.add_argument('-c','--count',action='store',dest='count',help='count type')
# parser.add_argument('-o','--out',action='store',dest='output',help='out_fn')


args       = parser.parse_args()
path    = args.input
count_type = args.count
# out_fn   = args.output
if count_type == 'count':
    cname = 'NumReads'
    out = 'counts'
elif count_type == 'TPM':
    cname = 'TPM'
    out = 'tpm'

quant_fns = [p for p in natsorted(glob.glob(f'{path}/*/quant.genes.sf'))]
dfs = []
for fn in quant_fns:
    df = pd.read_csv(fn,sep='\t',header=0,index_col=0)
    sp = '_'.join(fn.split('/')[-2].split('-'))
    df.rename(columns={cname:sp}, inplace=True)
    if count_type == 'count':
        df[sp] = df[sp].round(0).astype(int)
    dfs.append(df[[sp]])
    
dfs = pd.concat(dfs,axis=1)
dfs.to_csv(f'{path}/{out}.tsv',sep='\t')