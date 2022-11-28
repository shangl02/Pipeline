'''
Created by Shangzhong.Li@pfizer.com on 
'''
import os
import argparse
import pandas as pd
import subprocess
import gzip

parser = argparse.ArgumentParser(description='calculate PC for the cohort')

parser.add_argument('-r','--ref',action='store',dest='ref',help='ref 1k bfile')
parser.add_argument('--hld',action='store',dest='high_ld',help='high ld file')
parser.add_argument('--pre',action='store',dest='prefix',help='prefix of bfiles')

args = parser.parse_args()
k1_bfile = args.ref
high_LD_fn = args.high_ld
pre = args.prefix

def run(cmd):
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError:
        raise

#---------- calculate PC for the cohort ---------------------
bfile = f'{pre}.germ_newSp'
b_snps = pd.read_csv(f'{pre}.bim',sep='\t',header=None)[1].to_list()
k1_snps = pd.read_csv(f'{k1_bfile}.bim',sep='\t',header=None)[1].to_list()
# output common variants
common_snps = list(set(b_snps).intersection(k1_snps))
common_snp_fn = f'{pre}.common_snps.txt'
with open(common_snp_fn, 'w') as f:
    f.write('\n'.join(common_snps))
# #---------- extract bfile and k1 bfile -------------
inter_bfile = f'{pre}.germ_inter'
cmd = f'plink2 --bfile {bfile} --extract {common_snp_fn} --make-bed --out {inter_bfile}'
os.system(cmd)

inter_1kbfile = f'{pre}.inter_k1'
cmd = f'plink2 --bfile {k1_bfile} --extract {common_snp_fn} --make-bed --out {inter_1kbfile}'
os.system(cmd)
#------------ Dry run to identify SNPs that will fail ----------
fail = f'{pre}.merge_fail'
cmd = f'plink --bfile {inter_1kbfile} --bmerge {inter_bfile} --merge-mode 6 \
         --out {fail}'
os.system(cmd)
#-------------- Remove multi allels snps ------------------------
inter_rmSnp_bfile = f'{pre}.germ_inter_rmSNP'
cmd = f'plink2 --bfile {inter_bfile} --exclude {fail}.missnp \
      --make-bed --out {inter_rmSnp_bfile}'
os.system(cmd)
#-------------- Merge with 1kg data ---------------------------
merge = f'{pre}.merge'
cmd = f'plink --bfile {inter_rmSnp_bfile} --bmerge {inter_1kbfile} \
        --out {merge}'
os.system(cmd)
#-------------- Filter missing variants, rare variants and HWE , exclude high LD region --------------------
merge4prune = f'{pre}.merge4prune'
cmd = f'plink2 --bfile {merge} --geno 0.01 --maf 0.01 --max-maf 0.99 \
             --hwe 0.000001 --exclude bed0 {high_LD_fn} \
             --make-bed --out {merge4prune}'
os.system(cmd)
#-------------- LD pruning -----------------
prune = f'{pre}.prune'
cmd = f'plink2 --bfile {merge4prune} \
      --indep-pairwise 50 5 0.5 --out {prune}'
os.system(cmd)
cmd = f'plink2 --bfile {merge4prune} \
        --extract {prune}.prune.in --make-bed --out {prune}'
os.system(cmd)
#------------- extract 1kg for flashpca --------
k1g4pca = f'{pre}.k1g4pca'
cmd = f'plink2 --bfile {prune} --chr 1-22 --keep {k1_bfile}.fam \
        --make-bed --out {k1g4pca}'
os.system(cmd)
#------------- extract cohort for flashpca --------
b4pca = f'{pre}.b4pca'
cmd = f'plink2 --bfile {prune} --chr 1-22 --keep {bfile}.fam \
       --make-bed --out {b4pca}'
os.system(cmd)
#------------ run flashpca ------------------------
pca = f'{pre}.pca'
cmd = f'flashPCA --bfile {k1g4pca} --suffix _1kg --outload {pca}.loadings.txt --outmeansd {pca}.1kg_meansd.txt'
os.system(cmd)
#------------ project to cohort -------------------
print('project to cohort')
bpca = f'{pre}.bpca.txt'
cmd = f'flashPCA --bfile {b4pca} --project \
    --inmeansd {pca}.1kg_meansd.txt --outproj {bpca} \
    --inload {pca}.loadings.txt --suffix _study -v'
os.system(cmd)


#===== change SNP id to chr:pos:ref:alt to make them unambiguous ====
#---------------- update snp id to chr:pos:ref:alt ---------------------
bim_fn = f'{bfile}.bim'
bim_df = pd.read_csv(bim_fn,sep='\t',header=None)
bim_df[1] = bim_df.apply(lambda x:':'.join([str(x[0]),str(x[3]),x[5],x[4]]),axis=1)

newSp_bfile = f'{pre}.germ_newSp_id'
bim_df.to_csv(f'{newSp_bfile}.bim',sep='\t',index=False,header=None)
cmd = f'cp {bfile}.bed {newSp_bfile}.bed'
run(cmd)
cmd = f'cp {bfile}.fam {newSp_bfile}.fam'
run(cmd)
# get allele frequency
cmd = f'plink2 --bfile {newSp_bfile} \
     --freq cols=chrom,pos,ref,alt,altfreq,nobs --out {newSp_bfile}'
run(cmd)
#------ remove duplicated SNPS -----------------------
dedup_bfile = f'{pre}.germ_dedup'
cmd = f'plink2 --bfile {newSp_bfile} --rm-dup list force-first --make-bed --out {dedup_bfile}'
run(cmd)
# os.rename(f'{dedup_bfile}.bed',f'{newSp_bfile}.bed')
# os.rename(f'{dedup_bfile}.bim',f'{newSp_bfile}.bim')
# os.rename(f'{dedup_bfile}.fam',f'{newSp_bfile}.fam')
#-------------- transfer the plink format into table format ----
traw_pre = f'{pre}.germ_dedup'
cmd = f'plink2 --bfile {newSp_bfile} --recode A-transpose --out {traw_pre}'
run(cmd)
with open(f'{traw_pre}.traw') as in_f, gzip.open(f'{traw_pre}.traw.gz','wt') as out_f:
    head = in_f.readline().strip().split('\t')
    head = ['_'.join(c.split('_')[:2]) if '_' in c else c for c in head]
    out_f.write('\t'.join(head) + '\n') 
    for line in in_f:
        out_f.write(line)
os.remove(f'{traw_pre}.traw')
#-------------- transfer the traw.gz to bgzip format -----------
cmd = f'gunzip {traw_pre}.traw.gz' 
run(cmd)
cmd = f'bgzip {traw_pre}.traw'
run(cmd)
cmd = f'tabix -p bed -f -S 1 -b 4 -e 4 {traw_pre}.traw.gz'
run(cmd)




