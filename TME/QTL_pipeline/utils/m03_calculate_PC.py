'''
Created by Shangzhong.Li@pfizer.com on 
'''
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='calculate PC for the cohort')

parser.add_argument('-r','--ref',action='store',dest='ref',help='ref 1k bfile',default='/hpc/grid/hgcb/workspace/projects/P002_reference_information/plink/1000G_6PoP_MAF0.05')
parser.add_argument('--plink2',action='store',dest='plink2',help='plink2',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink2')
parser.add_argument('--plink1',action='store',dest='plink1',help='plink1',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink')
parser.add_argument('-p','--path',action='store',dest='path',help='path')
parser.add_argument('-f','--flashPCA',action='store',dest='flashPCA',help='flashPCA',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/software/flashpca_x86-64')


args = parser.parse_args()
work_path = args.path
k1_bfile = args.ref

plink2 = args.plink2
plink1 = args.plink1
flashPCA = args.flashPCA

bfile = f'{work_path}/plink/germ_newSp38'
# bfile = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WGS/plink/germ_newSp38'
# path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/tumor/Clinical/BLCA/WGS'
# k1_bfile = '/hpc/grid/hgcb/workspace/projects/P002_reference_information/plink/1000G_6PoP_MAF0.05'

#---------- calculate PC for the cohort ---------------------
plink_path = f'{work_path}/plink'
pca_path = f'{plink_path}/f02_pca'
os.makedirs(pca_path, exist_ok=True)

b_snps = pd.read_csv(f'{bfile}.bim',sep='\t',header=None)[1].to_list()
k1_snps = pd.read_csv(f'{k1_bfile}.bim',sep='\t',header=None)[1].to_list()
# output common variants
common_snps = list(set(b_snps).intersection(k1_snps))
common_snp_fn = f'{pca_path}/common_snps.txt'
with open(common_snp_fn, 'w') as f:
    f.write('\n'.join(common_snps))
# #---------- extract bfile and k1 bfile -------------
inter_bfile = f'{pca_path}/germ38'
cmd = f'{plink2} --bfile {bfile} --extract {common_snp_fn} --make-bed --out {inter_bfile}'
os.system(cmd)

inter_1kbfile = f'{pca_path}/k1'
cmd = f'{plink2} --bfile {k1_bfile} --extract {common_snp_fn} --make-bed --out {inter_1kbfile}'
os.system(cmd)
#------------ Dry run to identify SNPs that will fail ----------
fail = f'{pca_path}/merge_fail'
cmd = f'{plink1} --bfile {inter_1kbfile} --bmerge {inter_bfile} --merge-mode 6 \
         --out {fail}'
os.system(cmd)
#-------------- Remove multi allels snps ------------------------
inter_rmSnp_bfile = f'{pca_path}/germ_rmSNP38'
cmd = f'{plink2} --bfile {inter_bfile} --exclude {fail}.missnp \
      --make-bed --out {inter_rmSnp_bfile}'
os.system(cmd)
#-------------- Merge with 1kg data ---------------------------
merge = f'{pca_path}/merge'
cmd = f'{plink1} --bfile {inter_rmSnp_bfile} --bmerge {inter_1kbfile} \
        --out {merge}'
os.system(cmd)
#-------------- Filter missing variants, rare variants and HWE , exclude high LD region --------------------
merge4prune = f'{pca_path}/merge4prune'
cmd = f'{plink2} --bfile {merge} --geno 0.01 --maf 0.01 --hwe 0.000001 --exclude range /lustre/workspace/projects/ECD/bzhang2/data/high-LD-regions.txt --make-bed --out {merge4prune}'
os.system(cmd)
#-------------- LD pruning -----------------
prune = f'{pca_path}/prune'
cmd = f'{plink2} --bfile {merge4prune} \
      --indep-pairwise 1500 150 0.2 --out {prune}'
os.system(cmd)
cmd = f'{plink2} --bfile {merge4prune} \
        --extract {prune}.prune.in --make-bed --out {prune}'
os.system(cmd)
#------------- extract 1kg for flashpca --------
k1g4pca = f'{pca_path}/k1g4pca'
cmd = f'{plink2} --bfile {prune} --chr 1-22 --keep {k1_bfile}.fam \
        --make-bed --out {k1g4pca}'
os.system(cmd)
#------------- extract cohort for flashpca --------
b4pca = f'{pca_path}/b4pca'
cmd = f'{plink2} --bfile {prune} --chr 1-22 --keep {bfile}.fam \
       --make-bed --out {b4pca}'
os.system(cmd)
#------------ run flashpca ------------------------
os.chdir(pca_path)
pca = f'{pca_path}/pca'
cmd = f'{flashPCA} --bfile {k1g4pca} --suffix _1kg --outload {pca}_loadings.txt --outmeansd {pca}_1kg_meansd.txt'
os.system(cmd)
#------------ project to cohort -------------------
print('project to cohort')
bpca = f'{pca_path}/bpca.txt'
cmd = f'{flashPCA} --bfile {b4pca} --project \
    --inmeansd {pca}_1kg_meansd.txt --outproj {bpca} \
    --inload {pca}_loadings.txt --suffix _study -v'
os.system(cmd)






