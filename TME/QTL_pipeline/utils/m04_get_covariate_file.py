import pandas as pd
import os,glob
import gzip
import argparse

parser = argparse.ArgumentParser(description='prepare covaraite files')

parser.add_argument('-p','--path',action='store',dest='path',help='path')
parser.add_argument('--id_map',action='store',dest='id_map',help='id mapping files')
parser.add_argument('-d','--deconv',action='store',dest='deconv',help='deconvoluted cell type abundance file')

args = parser.parse_args()
work_path = args.path
meta_fn = args.id_map
raw_cyto_fn = args.deconv

pca_fn = f'{work_path}/plink/f02_pca/bpca.txt'
cov_fn = f'{work_path}/covariates.tsv'
final_cyto_fn = f'{work_path}/cyto.tsv'
#----------- prepare covariate files -----------------------
meta_df = pd.read_csv(meta_fn, sep='\t',header=0)
meta_df = meta_df[['sid','Sex','Age','Race','cytoid']].drop_duplicates()
meta_df['sid'] = meta_df['sid'].astype(str)
pca_df = pd.read_csv(pca_fn,header=0,sep='\t')
pca_df['sid'] = pca_df['IID'].map(lambda x: str(x.split('_')[1]))
pca_df = pd.merge(pca_df, meta_df, how='left',on='sid')
pca_df.drop(['sid','FID','Race','cytoid'],axis=1,inplace=True)
pca_df.to_csv(cov_fn,sep='\t',index=False)
#----------- prepare cytoreason files ----------------------
cyto_df = pd.read_csv(raw_cyto_fn, sep='\t', header=0,index_col=0)
cyto_df['cytoid'] = cyto_df.index.map(lambda x:x.split()[1])
cyto_df = pd.merge(cyto_df,meta_df[['cytoid','sid']],how='left',on='cytoid')
cyto_df.index = 'DNA_' + cyto_df['sid'].astype(str)
cyto_df.drop(['cytoid','sid'],axis=1,inplace=True)
cyto_df.to_csv(final_cyto_fn,sep='\t')


