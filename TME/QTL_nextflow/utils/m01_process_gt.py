import os
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='process genotype for downstream QTL run')

parser.add_argument('-i','--vcf',action='store',dest='vcf',help='input vcf file')
parser.add_argument('-p','--pre',action='store',dest='out_pre',help='prefix of outfiles')
parser.add_argument('--af',type=float,action='store',dest='af',help='allele frequency threshold',default=0.01)


args = parser.parse_args()
vcf = args.vcf
pre = args.out_pre
af_thred = args.af

def run(cmd):
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError:
        raise

#====== transfer blca vcf file to plink file  ======
cmd = f'plink2 --vcf {vcf} --allow-extra-chr --chr 1-22 \
       --double-id --make-bed \
       --vcf-half-call h --out {pre}'
run(cmd)

#====== QC of the variants ===========================
#------ get allele frequency of the variants ---------
cmd = f'plink2 --bfile {pre} --freq --out {pre}'
run(cmd)

#------ get missingness of the genotype --------------
cmd = f'plink2 --bfile {pre} --missing --out {pre}'
run(cmd)

#------ check missing rate ---------------------------
# check the missing rate
vmiss_fn = f'{pre}.vmiss'
vmiss_df = pd.read_csv(vmiss_fn,sep='\t',header=0)
plt.figure()
ax = vmiss_df['F_MISS'].hist(bins=20)
ax.set_title('distribution of missingness in variant level')
plt.savefig(f'{pre}.vmiss.png',dpi=300)

smiss_fn = f'{pre}.smiss'
smiss_df = pd.read_csv(smiss_fn,sep='\t',header=0)
plt.figure()
ax = smiss_df['F_MISS'].hist()
ax.set_title('distribution of missingness in sample level')
plt.savefig(f'{pre}.smiss.png',dpi=300)

#------- plot the frequency ---------------------------
freq_fn = f'{pre}.afreq'
freq_df = pd.read_csv(freq_fn, sep='\t',header=0)
print('there are', freq_df.shape[0],'variants in total',
      'among them',freq_df.query('OBS_CT == 2').shape[0],'are singleton')
plt.figure()
ax = freq_df['OBS_CT'].hist(bins=50)
ax.set_title('distribution of observed allele count')
ax.set_xlabel('number of allele count')
ax.set_ylabel('number of SNPs')
ax.set_yscale('log')
plt.savefig(f'{pre}.AC_dist.png',dpi=300)

plt.figure()
ax = freq_df.query('OBS_CT >=20')['ALT_FREQS'].hist(bins=50)
ax.set_title('distribution of Allele Frequency')
ax.set_xlabel('Allele frequency')
ax.set_ylabel('number of SNPs')
plt.savefig(f'{pre}.AF_dist.png',dpi=300)

#---------- get the variants for downstream analysis ----------
keep_snp_fn = f'{pre}_keep_snps.tsv'
freq_df.query('OBS_CT > 20 and (@af_thred < ALT_FREQS < 1 - @af_thred)')['ID'].to_csv(keep_snp_fn,index=False,header=None)
# extract the variants for downstream analysis
rmSNP_pfile = f'{pre}.germ_rmSNP'
cmd = f'plink2 --bfile {pre} --extract {keep_snp_fn} --make-bed --out {rmSNP_pfile}'
run(cmd)





