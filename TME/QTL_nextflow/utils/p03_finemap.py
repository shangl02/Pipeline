'''
Created by Shangzhong.Li@pfizer.com on 2022/03/25
Do finemappig for specific SNP in a specific QTL using Susie
'''

import argparse
import os,gzip
import tabix
import pandas as pd
from collections import defaultdict
from m10_annotate_clump_qtl import anno_clump


parser = argparse.ArgumentParser(description='finemapping using susieR')

parser.add_argument('-i','--snpid',action='store',dest='snp_id',help='significant snp id in format of chr:pos')
parser.add_argument('-c','--cell',action='store',dest='cell',help='cell type that you want to do finemapping')
parser.add_argument('-n','--name',action='store',dest='name',help='any name to indicate the trait of the loci', default='test')
parser.add_argument('-p','--path',action='store',dest='path',help='work path')
parser.add_argument('-f','--fmp',action='store',dest='finemap_path',help='finemap path')

parser.add_argument('--plink1',action='store',dest='plink1',
    help='plink path',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink')
parser.add_argument('--plink2',action='store',dest='plink2',
    help='plink path',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink2')

parser.add_argument('--susie',action='store',dest='susie',
    help='susie docker path',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/container/susie_0.11.42.sif')
parser.add_argument('--ref_anno',action='store',dest='ref_anno',
    help='ref annotation file')

args = parser.parse_args()

snp = args.snp_id
chrom,pos = snp.split(':')[:2]
pos = int(pos)
start = max(pos-5e5, 0)
end = pos + 5e5
cell = args.cell
name = args.name

work_path = args.path
finemap_path = args.finemap_path

raw_gwas_fn = f'{work_path}/matrixQTL/cell_qtl/{cell}.sort.txt.gz'
gt_ref_bed = f'{work_path}/plink/germ_newSp38.bed'

plink1 = args.plink1
plink2 = args.plink2
susie = args.susie
susie_code = f'{os.path.dirname(os.path.realpath(__file__))}/p03_susie.R'
ref_anno = args.ref_anno


def get_sub_gwas(chrom, lo, hi, gwas_fn, out_fn):
    '''extract the gwas hits in a region from the raw gwas
    summary stats file'''
    tb = tabix.open(gwas_fn)
    with gzip.open(gwas_fn,'rt') as in_f:
        head = in_f.readline().strip().split('\t')
    with gzip.open(out_fn,'wb') as out_h:
        line = '\t'.join(head) + '\n'
        out_h.write(line.encode())
        records = tb.query(str(chrom), int(lo), int(hi))
        for record in records:
            record_dic = {k:v for k,v in zip(head, record)}
            line = '\t'.join(record) + '\n'
            out_h.write(line.encode())

def get_ld_matrix(gt_ref, snp_fn, out_ld):
    '''this file generate ld matrix of the file'''
    if gt_ref.endswith('bed'):
        cmd = f'{plink1} --bfile {gt_ref[:-4]} --r square gz \
              --extract {snp_fn} --out {out_ld}'
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

def run_susie(susie, susie_code, gwas_fn, ld_fn, plot, out_fn):
    col = 't.stat'
    cmd = f'singularity run -B /:/media {susie} Rscript /media/{susie_code} \
            /media/{gwas_fn} /media/{ld_fn} /media/{plot} /media/{out_fn} {col}'
    os.system(cmd)
    
def run_finemap(snp, start, end, name, raw_gwas_fn, ref_anno):
    # get chrom
    chrom, pos = snp.split(':')[:2]
    nid = '_'.join(snp.split(':'))
    # get raw plink bed file
    lo, hi = start, end
    # 1. get sub gwas region
    snp_path = f'{finemap_path}/{name}_{nid}'
    os.makedirs(snp_path,exist_ok=True)
    sub_gwas_fn = snp_path + '/' + nid + '.tsv.gz'
    get_sub_gwas(chrom,lo,hi,raw_gwas_fn, sub_gwas_fn)

    # 2. get snps file
    snp_id_fn = snp_path + '/' + '_'.join(snp.split(':')) + '.snp.txt'
    df = pd.read_csv(sub_gwas_fn, sep='\t',header=0,compression='gzip')
    df['SNP'].to_csv(snp_id_fn,sep='\t',header=0,index=False)

    # 3. calculate LD matrix                        
    ld_pre = snp_path + '/' + '_'.join(snp.split(':'))
    get_ld_matrix(gt_ref_bed, snp_id_fn, ld_pre)

    # 4. run susie
    ld_fn = ld_pre + '.ld.gz'
    plot = snp_path + '/' + '_'.join(snp.split(':')) + '.jpg'
    cs_fn = snp_path + '/' + '_'.join(snp.split(':')) + '.cs.txt'
    run_susie(susie, susie_code, sub_gwas_fn, ld_fn, plot, cs_fn)
    # 5. annotate sub gwas results
    tb = tabix.open(ref_anno)
    with gzip.open(ref_anno, 'rt') as in_f:
        columns = in_f.readline().strip().split('\t')
        anno_chrom = in_f.readline()

    cs_df = pd.read_csv(cs_fn,sep='\t',header=0)
    cs_dict = {}
    for idx, row in cs_df.iterrows():
        for i in str(row['variable']).split(','):
            cs_dict[int(float(i))-1] = row['cs']
    idxes = list(cs_dict.keys())

    df = pd.read_csv(sub_gwas_fn,sep='\t',header=0,compression='gzip')
    df[['Gene','near_50kgene','Consequence','rsid','gnomad','CADD']] = \
             df.apply(lambda row:anno_clump(row, tb, columns, anno_chrom), axis=1)

    anno_fn = sub_gwas_fn.split('.')[0] + '.anno.tsv'
    df.to_csv(anno_fn,sep='\t')


run_finemap(snp, start, end, name, raw_gwas_fn, ref_anno)