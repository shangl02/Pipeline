'''
Created by Shangzhong.Li@pfizer.com on 2021/11/18
This file does liftover for GWAS summary stats

-------- parameter explanations ------------------
# -i: input gwas file, can be gzip or not, needs to be full path
# -o: output gwas file, output file will be gzip to save space, needs to be full path.
# -t: temparay folder to store some intermediate files, needs to be full path
# -c: header in the input gwas file that indicates chromosome
# -p: header in the input gwas file that indicates position
# --conversion: 38to37 or 37to38
'''
import os,gzip
import pandas as pd
import glob
import argparse


parser = argparse.ArgumentParser(description='merge GWAS results')

parser.add_argument('-i','--input',action='store',dest='input',help=' file name')
parser.add_argument('-t','--temp',action='store',dest='temp',help='temprary path')
parser.add_argument('-o','--output',action='store',dest='out',help='out file')
parser.add_argument('-c','--chr',action='store',dest='chrom',help='header in file that indicates chromosome')
parser.add_argument('-s','--start',action='store',dest='start',help='header in file that indicates start')
parser.add_argument('-e','--end',action='store',dest='end',help='header in file that indicates end')
parser.add_argument('--conversion',action='store',dest='conver',help='conversion type')
parser.add_argument('--chain',action='store',dest='chain',help='chain type',default='ensembl')

args        = parser.parse_args()
fn          = args.input
temp_path   = args.temp
out_fn      = args.out

chrom_head  = args.chrom
start       = args.start
end         = args.end
conver_type = args.conver
chain_sour  = args.chain

if fn.startswith('/'):
    fn = f'/media/{fn}'
if out_fn.startswith('/'):
    out_fn = f'/media/{out_fn}'
if temp_path.startswith('/'):
    temp_path = f'/media/{temp_path}'
os.makedirs(temp_path, exist_ok = True)


if chain_sour == 'ensembl':
    if conver_type == '37to38':
        chain_fn = '/GRCh37_to_GRCh38.chain.gz'
    elif conver_type == '38to37':
        chain_fn = '/GRCh38_to_GRCh37.chain.gz'
elif chain_sour == 'UCSC':
    if conver_type == '19to38':
        chain_fn = '/hg19ToHg38.over.chain.gz'
    elif conver_type == '38to19':
        chain_fn = '/hg38ToHg19.over.chain.gz'


def get_new_col(columns):
    new_columns = [chrom_head, start, end] + \
      [h for h in columns if h not in [chrom_head, start,end]]
    return new_columns

def gwas2bed(gwas_fn, temp_path):
    head = True
    bed_fn = f'{temp_path}/temp_' + gwas_fn.split('/')[-1] + '.gz'
    if gwas_fn.endswith('.gz'):
        for df in pd.read_csv(gwas_fn,sep='\t',header=0,compression='gzip',chunksize=1e5):
            columns = df.columns.tolist()
            new_cols = get_new_col(columns)
            if head:
                df[new_cols].to_csv(bed_fn,sep='\t',index=False,compression='gzip')
                head = False
            else:
                df[new_cols].to_csv(bed_fn,sep='\t',index=False,compression='gzip',mode='a',header=None)
    else:
        for df in pd.read_csv(gwas_fn,sep='\t',header=0,chunksize=1e5):
            columns = df.columns.tolist()
            new_cols = get_new_col(columns)
            if head:
                df[new_cols].to_csv(bed_fn,sep='\t',index=False,compression='gzip')
                head = False
            else:
                df[new_cols].to_csv(bed_fn,sep='\t',index=False,compression='gzip',mode='a',header=None)
    return bed_fn


def run_crossmap(bed_fn, out_fn, chain):
    if '/' in bed_fn:
        head_fn = os.path.dirname(bed_fn) + '/' + bed_fn.split('/')[-1].split('.')[0] + '.head.txt'
    else:
        head_fn = bed_fn + '.head.txt'
    # 1. get head
    with gzip.open(bed_fn,'rt') as in_h, open(head_fn,'w') as out_h:
        fst_line = in_h.readline()
        out_h.write(fst_line)
    # 2. cross map
    if out_fn.endswith('.gz'):
        out_fn = out_fn[:-3]
    cmd = ('CrossMap.py bed {chain} {infn} {out}').format(chain=chain,infn=bed_fn,out=out_fn)
    os.system(cmd)
    # print(cmd)
    # 3. sort temp folder
    temp = head_fn[:-9]
    os.makedirs(temp, exist_ok=True)
    cmd = ('sort -k1,1 -k2,2n -T {t} {out} | \
           cat {head} - | bgzip > {out}.gz && \
           tabix -p bed -f -S 1 -b 2 -e 2 {out}.gz && \
           rm {infn} {out} {head} && rm -r {t}').format(t=temp,out=out_fn,head=head_fn,infn=bed_fn)
    os.system(cmd)
    # print(cmd)



print('step1: transfer to bed format')
temp_fn = gwas2bed(fn, temp_path)
print('step1 finished')
print('step2: liftover')
run_crossmap(temp_fn, out_fn, chain_fn)
print('step2 finished')