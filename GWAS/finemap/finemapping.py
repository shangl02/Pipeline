'''
Created by Shangzhong.Li@pfizer.com on 2022/06/01
This file does finemapping using SusieR

-------- parameter explanations ------------------
# -i: input gwas file, can be gzip or not, needs to be full path
# -r: reference genotype file for calculating LD in plink format, needs to be the binary file which ends with bed/bgen.
# -s: column name indicating SNP id in the input gwas file.
# -o: output path that will save all the results
# --stats: column name in gwas file indicating the stats of the results
'''
import os
import gzip
import pandas as pd
import glob
import argparse


parser = argparse.ArgumentParser(description='finemapping')


parser.add_argument('-i','--input',action='store',dest='input',help='gwas result file name')
parser.add_argument('-r','--ref',action='store',dest='ref',help='reference genotype file for calculating LD in plink format, needs to be the binary file which ends with bed/bgen')
parser.add_argument('-s','--snp',action='store',dest='snp',help='column in gwas file indicating the SNP id',default='SNP')
parser.add_argument('-o','--out',action='store',dest='out',help='output path')
parser.add_argument('--stats',action='store',dest='stats',help='column name in gwas file indicting the stats of the results')

args        = parser.parse_args()
gwas_fn     = f'/media/{args.input}'
gt_ref      = f'/media/{args.ref}'
out_path    = f'/media/{args.out}'
os.makedirs(f'{out_path}',exist_ok=True)
snp         = args.snp
stats       = args.stats

plink1      = '/opt/plink1'
plink2      = '/opt/plink2'


def get_ld_matrix(gt_ref, snp_fn, out_ld):
    '''this file generate ld matrix of the file'''
    if gt_ref.endswith('bed'):
        cmd = f'{plink1} --bfile {gt_ref[:-4]} --r square gz \
              --extract {snp_fn} --out {out_ld}'
        print(cmd)
        os.system(cmd)
    elif gt_ref.endswith('bgen'):
        sample_fn = gt_ref[:-4] + 'sample'
        # 1. transfer to bed file format
        bed_pre = out_ld
        cmd = f'{plink2} --bgen {gt_ref} ref-first \
               --sample {sample_fn} --extract {snp_fn} \
               --make-bed --out {bed_pre}'
        os.system(cmd)
        # 2. calculate ld matrix
        cmd = f'{plink1} --bfile {bed_pre} --r square gz \
              --out {out_ld}'
        # print(cmd)
        os.system(cmd)


def run_susie(gwas_fn, gt_ref, out_path, snp, stats):
    if gt_ref.endswith('bed'):
        prefix = gt_ref[:-4]
    elif gt_ref.endswith('bgen'):
        prefix = gt_ref[:-5]
    ref_snp_fn = f'{prefix}.bim'
    sub_gwas_fn = f'{out_path}/sub_gwas.txt.gz'
    common_snp_fn = f'{out_path}/snp.txt'
    out_ld = f'{out_path}/ld_mtx'

    # 1. get ref snps
    print('get ref snps')
    ref_snps = pd.read_csv(ref_snp_fn,sep='\t',header=None,usecols=[1])[1].tolist()
    # 2. get gwas snps and overlap with ref snps
    print('get gwas and ref overlapped snps')
    if gwas_fn.endswith('.gz'):
        gwas_df = pd.read_csv(gwas_fn,sep='\t',header=0,compression='gzip')
    else:
        gwas_df = pd.read_csv(gwas_fn,sep='\t',header=0)
    snps = gwas_df[snp].tolist()
    common_snps = list(set(snps).intersection(ref_snps))
    print(len(common_snps))
    # 3. subset of the gwas results
    print('get subset ofthe gwas results')
    sub_df = gwas_df[gwas_df[snp].isin(common_snps)]
    sub_df = sub_df[~sub_df[snp].duplicated()]
    sub_df.to_csv(sub_gwas_fn,sep='\t',compression='gzip')
    with open(common_snp_fn,'w') as f:
        f.write('\n'.join(sub_df[snp]))
    # 4. get ld matrix
    print('get ld matrix')
    get_ld_matrix(gt_ref, common_snp_fn, out_ld)
    # 5. run susieR
    print('run susieR')
    cmd = f'Rscript /opt/susie.R {sub_gwas_fn} {out_ld}.ld.gz {out_path}/result.png {out_path}/susie_result.txt {stats}'
    os.system(cmd)


run_susie(gwas_fn, gt_ref, out_path, snp, stats)