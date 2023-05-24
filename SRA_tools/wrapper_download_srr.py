import pandas as pd
from pysradb import SRAweb
import sys
import os
import argparse
import glob

code_path = os.path.dirname(os.path.realpath(__file__))
download_srr_code = f'{code_path}/download_srr.py'


parser = argparse.ArgumentParser(description='download srr')

parser.add_argument('-i','--gse_design',action='store',dest='input',help='gse design file')
parser.add_argument('-o','--out_dir',action='store',dest='out_dir',help='out directory')

args = parser.parse_args()
gse_design_fn = args.input
out_dir = args.out_dir

os.chdir(out_dir)

db = SRAweb()

def get_srx_srr_id_map(gse_design_fn):
    '''get srr id from the gse design fn'''
    df = pd.read_excel(gse_design_fn, sheet_name='Design')
    srx_gsm_dict = df.set_index('SrxID')['GSM'].to_dict()
    dfs = []
    for srx_id in df['SrxID']:
        dfs.append(db.sra_metadata(srx_id))
    meta_df = pd.concat(dfs)

    return meta_df, srx_gsm_dict
    

def download_wrapper(out_dir, srx_srr_dict, srx_gsm_dict):
    for srx, srrs in srx_srr_dict.items():
        gsm = srx_gsm_dict[srx]
        srr_ids = ','.join(srrs)
        cmd = f'python {download_srr_code} -s {srr_ids} -g {gsm} -o {out_dir}'
        bsub = f'bsub -q medium -o {out_dir}/{gsm}.log \"{cmd}\"'
        os.system(bsub)

            

# get srx and srr id mappings
print('get srr id from srx id')
srx_srr_df, srx_gsm_dict = get_srx_srr_id_map(gse_design_fn)
print('finished: get srr id from srx id')
# output meta_df
srx_srr_id_map_fn = os.path.splitext(gse_design_fn)[0] + '.srx2srr.tsv'
srx_srr_df.to_csv(srx_srr_id_map_fn, sep='\t', index=False)
# build {srx: [srr ids]}
srx_srr_dict = {k:list(v) for k,v in srx_srr_df.groupby('experiment_accession')['run_accession']}

# download
print('sub mit download')
download_wrapper(out_dir, srx_srr_dict, srx_gsm_dict)








