'''
This file aggregate transcript expression values to gene level
'''

import glob, os
import pandas as pd
import re
import argparse


parser = argparse.ArgumentParser(description='aggregate transcript expression level to gene level')

parser.add_argument('-i','--input',action='store',dest='input',help='input rna expression file')
parser.add_argument('-o','--out',action='store',dest='output',help='out gene level expression file')
parser.add_argument('-g','--gff',action='store',dest='gff',help='gff file')


args      = parser.parse_args()
rna_fn    = args.input
gene_fn   = args.output
gff_fn    = args.gff


def get_id(anno, dtype='rna'):
    try:
        if dtype == 'rna':
            return re.search('(?<=transcript_id=).+?(?=[;,])',anno).group(0)
        elif dtype == 'gene':
            return re.search('(?<=Parent=gene:).+?(?=[;,])',anno).group(0)
    except:
        return ''

# build {gene_id:rna_id} from gff file
gff_df = pd.read_csv(gff_fn, sep='\t', header=None, comment='#', low_memory=False)
gff_df['gene_id'] = gff_df[8].map(lambda x: re.search('',x))
gff_df['gene_id'] = gff_df[8].map(lambda x: get_id(x,'gene'))
gff_df['rna_id'] = gff_df[8].map(lambda x: get_id(x,'rna'))
gff_df = gff_df[['gene_id','rna_id']].drop_duplicates()
gff_df.index = gff_df['rna_id']


df = pd.read_csv(rna_fn, sep='\t', header=0, index_col=0)
df.index = df.index.map(lambda x: x.split('.')[0])

df = pd.merge(df, gff_df['gene_id'].to_frame(), how='left', left_index=True, right_index=True)
gene_expr_df = df.groupby('gene_id').sum()
gene_expr_df.to_csv(gene_fn, sep='\t')