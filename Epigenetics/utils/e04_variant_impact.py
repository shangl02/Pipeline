'''
Created by Shangzhong.Li@pfizer.com on 2022/06/27
'''

import argparse
import os
import json
import sys
import pandas as pd
from pybedtools import BedTool

code_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(description='predict snp effect on openess of chromatine')


parser.add_argument('-c','--cell',action='store',dest='cell',help='cell name')
parser.add_argument('-p','--path',action='store',dest='path',help='root path')
parser.add_argument('-s','--snp',action='store',dest='snp_fn',help='snp file with each line as a snp id in format chr:pos:ref:alt')
parser.add_argument('-o','--out',action='store',dest='out_fn',help='out path')

args        = parser.parse_args()

cell        = args.cell
root_path   = args.path
snp_fn      = args.snp_fn
out_fn    = args.out_fn

if os.path.isabs(out_fn):
    out_path = os.path.dirname(out_fn)
else:
    out_path = os.getcwd()

model_path = os.path.join(root_path, 'p03_cell_model') 
peak_path = os.path.join(root_path, 'f04_cell_peak','adult') 
ref_fa = os.path.join(root_path, 'genome', 'hg38.fasta') 
peak_file = os.path.join(peak_path, cell + '.bed.gz') 
model_file = os.path.join(model_path, cell + '_model.pt') 


# prepare configure
if not model_file:
    raise 'model file  does not exist'
if not snp_fn:
    raise 'snp file does not exist'
if not ref_fa:
    raise 'reference fa file does not exist'
if not peak_file:
    raise 'peak file does not exist'


#============== prepare input for the model given snpid ======================
def prepare_pred_input(snp_fn, peak_file, out_path):
    # prepare snp bed file
    snp_df = pd.read_csv(snp_fn, names=['snp'], header=None)
    snp_df[['chr','pos','ref','alt']] = snp_df.apply(lambda x: pd.Series(x['snp'].split(':')), axis=1)
    snp_df['pos'] = snp_df['pos'].astype(int)
    snp_df['start'] = snp_df['pos'] - 1
    snp_bed = BedTool.from_dataframe(snp_df[['chr','start','pos','ref','alt','snp']])
    # prepare peak bed file
    peak_df = pd.read_csv(peak_file, sep='\t', header=None, compression='gzip')
    peak_bed = BedTool.from_dataframe(peak_df[[0,1,2]])
    # overlap between peak and snp
    inter_df = peak_bed.intersect(snp_bed, wa=True, wb=True).to_dataframe(names=['peak_chr','peak_s','peak_e','snp_chr','snp_s','snp_e','ref','alt','snp'])
    # output file
    if inter_df.empty:
        print('no overlap between query snps and human atlas peaks')
        sys.exit(1)
    else:
        inter_df['peak_s'] -= 300
        inter_df['peak_e'] += 300
        cols = ['peak_chr','peak_s','peak_e','snp_chr','snp_s','snp_e','snp','alt','ref']
        pred_in_fn = os.path.join(out_path, cell + '_snp4prediction.tsv')
        inter_df[cols].to_csv(pred_in_fn,sep='\t',header=None,index=False)
    return pred_in_fn


pred_in_fn = prepare_pred_input(snp_fn, peak_file, out_path)
model_name = os.path.splitext(os.path.basename(model_file))[0]
#============== predict the impact of variants ======================
# create configure file for variant impact prediction
config_fn = os.path.join(out_path, cell + '_variant_config.json')
config = {
  "peak_snp_bed": pred_in_fn,
  "ref_fasta": ref_fa,
  "model_name": "scEpiLock",
  "model_weight_path": model_file,
  "n_class": 1,
  "result_path": os.path.join(out_path,cell + '_snp_pred_result')
}

with open(config_fn,'w') as f:
    json.dump(config, f)
epilock_variant_code = os.path.join(code_path, 'scEpiLock', 'variant_impact_module', 'main.py')

# variant impact
cmd = f'python {epilock_variant_code} {config_fn}'
print(cmd)
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
    pred_snp_fn = os.path.join(path, 'diff_yhat.csv')
    pred_snp_df = pd.read_csv(pred_snp_fn,sep=',',header=0,index_col=0)
    pred_snp_df = get_fwd_rev_pred(pred_snp_df)
    pred_snp_df.columns = ['diff_fwd','diff_rev']
    pred_snp_df['diff_max'] = pred_snp_df.apply(lambda x: x[x.abs().idxmax()],axis=1)
    pred_snp_df *= -1
    # wt
    wt_snp_fn = os.path.join(path, 'wt_yhat.csv')
    wt_snp_df = pd.read_csv(wt_snp_fn,sep=',',header=0,index_col=0)
    wt_snp_df = get_fwd_rev_pred(wt_snp_df)
    wt_snp_df.columns = ['wt_fwd','wt_rev']
    wt_snp_df['wt_max'] = wt_snp_df.apply(lambda x: x[x.abs().idxmax()],axis=1)
    # mut
    mut_snp_fn = os.path.join(path, 'mutant_yhat.csv')
    mut_snp_df = pd.read_csv(mut_snp_fn,sep=',',header=0,index_col=0)
    mut_snp_df = get_fwd_rev_pred(mut_snp_df)
    mut_snp_df.columns = ['mut_fwd','mut_rev']
    mut_snp_df['mut_max'] = mut_snp_df.apply(lambda x: x[x.abs().idxmax()],axis=1)
    # merge
    df = pd.concat([wt_snp_df,mut_snp_df,pred_snp_df], axis=1)
    return df


snp_df = pd.read_csv(pred_in_fn,sep='\t',header=None,names=['peak_chr','peak_s','peak_e',
                'snp_chr','snp_s','snp_e','snpid','alt','ref'])
res_path = os.path.join(out_path, cell + '_snp_pred_result')
df = merge_wt_mut_pred(res_path)
merge_df = pd.concat([snp_df,df], axis=1)
merge_df = merge_df.drop_duplicates()
merge_df.to_csv(out_fn,sep='\t',index=False)