import pandas as pd
from pysradb import SRAweb
import sys
import os
import argparse
import glob



parser = argparse.ArgumentParser(description='download srr')

parser.add_argument('-s','--srr',action='store',dest='srr',help='srr id, can be multiple ids separated by ,')
parser.add_argument('-g','--gsm',action='store',dest='gsm',help='gsm id')
parser.add_argument('-o','--out_dir',action='store',dest='out_dir',help='out directory')

args = parser.parse_args()
srr_ids = [s.strip() for s in (args.srr).split(',')]
gsm     = args.gsm
out_dir = args.out_dir

os.chdir(out_dir)

db = SRAweb()



def download_srr(srr_ids, gsm, out_dir):
    temp = f'{out_dir}/temp'
    os.makedirs(temp, exist_ok=True)
    all_fq_files = []
    # download
    for srr in srr_ids:
        # fasterq-dump
        cmd = f'fasterq-dump {srr} -O {out_dir} -t {temp} --progress --split-files'
        os.system(cmd)
        # gzip 
        fastq_files = sorted(glob.glob(f'{out_dir}/{srr}*.fastq'))
        for fq in fastq_files:
            cmd = f'gzip {fq}'
            os.system(cmd)
        gzip_fq_files = [f'{fq}.gz' for fq in fastq_files]
        all_fq_files.append(gzip_fq_files)
    # merge fastq files if one srx has more than 1 srr ids, this case usually applies to WGS level data.
    if len(all_fq_files[0]) == 1:
        fst_fq_files = ' '.join([f[0] for f in all_fq_files])
        merge_fst_fq_file = f'{out_dir}/{gsm}.fq.gz'
        os.system(f'cat {fst_fq_files} > {merge_fst_fq_file}')
    if len(all_fq_files[0]) == 2:
        # first read
        fst_fq_files = ' '.join([f[0] for f in all_fq_files])
        merge_fst_fq_file = f'{out_dir}/{gsm}_R1.fq.gz'
        os.system(f'cat {fst_fq_files} > {merge_fst_fq_file}')
        # second read
        snd_fq_files = ' '.join([f[1] for f in all_fq_files])
        merge_snd_fq_file = f'{out_dir}/{gsm}_R2.fq.gz'
        os.system(f'cat {snd_fq_files} > {merge_snd_fq_file}')
    for fq_lst in all_fq_files:
        for fq in fq_lst:
            os.remove(fq)

    
download_srr(srr_ids, gsm, out_dir)
