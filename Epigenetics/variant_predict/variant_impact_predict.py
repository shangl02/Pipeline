'''
Created by Shangzhong.Li@pfizer.com on 2022/06/27
'''

import argparse
import os
import json
import sys
import pandas as pd

code_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(description='predict snp effect on openess of chromatine')


parser.add_argument('-m','--model',action='store',dest='model',help='path of the pytorch model file')
parser.add_argument('-s','--snp',action='store',dest='snp_bed',help='bed file with snp information')
parser.add_argument('-r','--ref',action='store',dest='ref_fa',help='reference fa file')
parser.add_argument('-o','--out',action='store',dest='out',help='output of prediction results')

args    = parser.parse_args()
model   = args.model
snp_fn  = args.snp_bed
ref_fa  = args.ref_fa
out_fn  = args.out


model_name = os.path.splitext(os.path.basename(model))[0]
#============== predict the impact of variants ======================
# prepare configure
if not model:
    raise 'model file  does not exist'
if not snp_fn:
    raise 'snp file does not exist'
if not ref_fa:
    raise 'reference fa file does not exist'

model_path = os.path.dirname(model)
# create configure file for variant impact prediction
config_fn = f'{model_path}/variant_config.json'
config = {
  "peak_snp_bed": snp_fn,
  "ref_fasta": ref_fa,
  "model_name": "scEpiLock",
  "model_weight_path": model,
  "n_class": 1,
  "result_path": f"{model_path}/snp_pred_result"
}

with open(config_fn,'w') as f:
    json.dump(config, f)
epilock_variant_code = f'{code_path}/scEpiLock/variant_impact_module/main.py'

# variant impact
cmd = f'python {epilock_variant_code} {config_fn}'
os.system(cmd)


#================== reformat the results =======================
def get_fwd_rev_pred(df):
    # split foward and reverse prediction into two columns
    n = int(df.shape[0] / 2)
    fwd_df = df.loc[0:n-1,:]
    rev_df = df.loc[n:,:].reset_index(drop=True)
    df = pd.concat([fwd_df, rev_df], axis=1, ignore_index=True)
    return df

def merge_wt_mut_pred(path):
    # diff
    pred_snp_fn = f'{path}/diff_yhat.csv'
    pred_snp_df = pd.read_csv(pred_snp_fn,sep=',',header=0,index_col=0)
    pred_snp_df = get_fwd_rev_pred(pred_snp_df)
    pred_snp_df.columns = ['diff_fwd_','diff_fwd']
    pred_snp_df['diff_mean'] = pred_snp_df.mean(axis=1)
    # wt
    wt_snp_fn = f'{path}/wt_yhat.csv'
    wt_snp_df = pd.read_csv(wt_snp_fn,sep=',',header=0,index_col=0)
    wt_snp_df = get_fwd_rev_pred(wt_snp_df)
    wt_snp_df.columns = ['wt_fwd','wt_rev']
    wt_snp_df['wt_mean'] = wt_snp_df.mean(axis=1)
    # mut
    mut_snp_fn = f'{path}/mutant_yhat.csv'
    mut_snp_df = pd.read_csv(mut_snp_fn,sep=',',header=0,index_col=0)
    mut_snp_df = get_fwd_rev_pred(mut_snp_df)
    mut_snp_df.columns = ['mut_fwd','mut_fwd']
    mut_snp_df['mut_mean'] = mut_snp_df.mean(axis=1)
    # merge
    df = pd.concat([wt_snp_df,mut_snp_df,pred_snp_df], axis=1)
    return df


snp_df = pd.read_csv(snp_fn,sep='\t',header=None,names=['peak_chr','peak_s','peak_e',
                'snp_chr','snp_s','snp_e','rsid','ref','alt'])
df = merge_wt_mut_pred(f'{model_path}/snp_pred_result')
merge_df = pd.concat([snp_df,df], axis=1)
merge_df.to_csv(out_fn,sep='\t',index=False)



