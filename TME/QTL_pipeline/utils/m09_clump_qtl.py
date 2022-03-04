'''
Created by Shangzhong.Li@pfizer.com on 2022/02/28
Clump the QTL results
'''
import os
import argparse
import glob

parser = argparse.ArgumentParser(description='Clump the QTL results')

# parser.add_argument('-i','--qtl',action='store',dest='qtl',help='path that have qtl files')
# parser.add_argument('-o','--out',action='store',dest='out',help='out path')

parser.add_argument('-p','--path',action='store',dest='path',help='path of project')
parser.add_argument('--pval',type=float,action='store',dest='pval',help='pvalue threshold',default=5e-8)

# parser.add_argument('--ref',action='store',dest='ref',help='ref plink file for clumping',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WGS/plink/germ_newSp38')
parser.add_argument('--plink2',action='store',dest='plink2',help='plink2',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink2')
parser.add_argument('--plink1',action='store',dest='plink1',help='plink1',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink')

args = parser.parse_args()


work_path = args.path
qtl_path = f'{work_path}/matrixQTL/cell_qtl'
qtl_fns = sorted(glob.glob(f'{qtl_path}/*.txt.gz'))
out_path = f'{work_path}/matrixQTL/clump'
os.makedirs(out_path, exist_ok=True)
pval = args.pval
ref_bfile = f'{work_path}/plink/f02_pca/germ_rmSNP38'

plink1 = args.plink1
plink2 = args.plink2

for qtl_fn in qtl_fns:
    out_fn = f'{out_path}/' + '.'.join(qtl_fn.split('/')[-1].split('.')[:-2])
    cmd = f'{plink1} --bfile {ref_bfile} --clump {qtl_fn} \
          --clump-snp-field SNP \
          --clump-field p-value \
          --clump-p1 5e-8 --clump-p2 1 \
          --clump-r2 0.2 --clump-kb 250 \
          --out {out_fn}'
    os.system(cmd)
    # print(cmd)

