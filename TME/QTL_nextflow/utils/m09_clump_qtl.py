'''
Created by Shangzhong.Li@pfizer.com on 2022/02/28
Clump the QTL results
'''
import os
import argparse
import glob

parser = argparse.ArgumentParser(description='Clump the QTL results')


parser.add_argument('--pval',type=float,action='store',dest='pval',help='pvalue threshold',default=5e-8)
parser.add_argument('--pre',action='store',dest='prefix',help='study prefix')

# parser.add_argument('--ref',action='store',dest='ref',help='ref plink file for clumping',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WGS/plink/germ_newSp38')
args = parser.parse_args()

pre = args.prefix
pval = float(args.pval)


qtl_fns = sorted(glob.glob(f'{pre}*.txt.gz'))
ref_bfile = f'{pre}.germ_dedup'


for qtl_fn in qtl_fns:
    out_fn = '.'.join(qtl_fn.split('/')[-1].split('.')[:-3])
    cmd = f'plink --bfile {ref_bfile} --clump {qtl_fn} \
          --clump-snp-field SNP \
          --clump-field p-value \
          --clump-p1 {pval} --clump-p2 1 \
          --clump-r2 0.2 --clump-kb 250 \
          --out {out_fn}'
    os.system(cmd)
    # print(cmd)

