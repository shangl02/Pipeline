import os, gzip, subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='get samples that are duplicated')

parser.add_argument('-p','--pre',action='store',dest='prefix',help='path')
parser.add_argument('--id_map',action='store',dest='id_map',help='id mapping file between vcf, covariates and cyto reason')

args = parser.parse_args()
pre  = args.prefix
id_map_fn = args.id_map

bfile = f'{pre}.germ_rmSNP'
def run(cmd):
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError:
        raise
#---------- remove the duplicated samples ---------------------
smiss_fn = f'{pre}.smiss'
smiss_df = pd.read_csv(smiss_fn, header=0, delim_whitespace=True)
smiss_df['id'] = smiss_df['IID'].map(lambda x: '_'.join(x.split('_')[:2]))
smiss_df = smiss_df[smiss_df.groupby('id')['F_MISS'].transform(min) == smiss_df['F_MISS']]
# keep samples
keep_sp = f'{pre}.keep_sps.txt'
smiss_df[['#FID','IID']].to_csv(keep_sp, sep='\t',index=False)

rmSp_bfile = f'{pre}.germ_rmSNP_rmSp'
cmd = f'plink2 --bfile {bfile} --keep {keep_sp} --make-bed --out {rmSp_bfile}'
run(cmd)
#---------------- update sample names ---------------------
update_sp_fn = f'{pre}.update_sample.txt'
id_map_df = pd.read_csv(id_map_fn,sep='\t',header=0)
id_map_df['sid'] = 'DNA_' + id_map_df['sid'].astype('str')
vcfid_sid_dic = id_map_df.set_index('vcfid')['sid'].to_dict()
smiss_df['id'] = smiss_df['IID'].map(lambda x: vcfid_sid_dic[x])
smiss_df[['IID','IID','id','id']].to_csv(update_sp_fn, sep='\t',index=False,header=None)
newSp_bfile = f'{pre}.germ_newSp'
cmd = f'plink2 --bfile {rmSp_bfile} --update-ids {update_sp_fn} \
         --make-bed --out {newSp_bfile}'
run(cmd)