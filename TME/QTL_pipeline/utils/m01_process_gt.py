import os
import subprocess
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('ggplot')

parser = argparse.ArgumentParser(description='process genotype for downstream QTL run')

parser.add_argument('-i','--vcf',action='store',dest='vcf',help='input vcf file')
parser.add_argument('--plink2',action='store',dest='plink2',help='plink2',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink2')
parser.add_argument('--plink1',action='store',dest='plink1',help='plink1',default='/hpc/grid/wip_drm_targetsciences/users/shangzhong/software/plink')
parser.add_argument('-p','--path',action='store',dest='path',help='path')
parser.add_argument('--af',type=float,action='store',dest='af',help='allele frequency threshold',default=0.01)

args = parser.parse_args()
vcf = args.vcf
path = args.path

plink2 = args.plink2
plink1 = args.plink1
af_thred = args.af

os.makedirs(path,exist_ok=True)
def run(cmd):
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError:
        raise

#====== transfer blca vcf file to plink file  ======
plink_path = f'{path}/plink'
os.makedirs(plink_path,exist_ok=True)

pfile = f'{plink_path}/germ38'
cmd = f'{plink2} --vcf {vcf} --allow-extra-chr --chr 1-22 \
       --double-id --make-bed \
       --vcf-half-call h --out {pfile}'
run(cmd)

#====== QC of the variants ===========================
#------ get allele frequency of the variants ---------
qc_path = f'{plink_path}/f01_qc'
os.makedirs(qc_path,exist_ok=True)
cmd = f'{plink2} --bfile {pfile} --freq --out {qc_path}/freq38'
run(cmd)

#------ get missingness of the genotype --------------
missing = f'{qc_path}/germ38'
cmd = f'{plink2} --bfile {pfile} --missing --out {missing}'
run(cmd)

#------ check missing rate ---------------------------
# check the missing rate
vmiss_fn = f'{qc_path}/germ38.vmiss'
vmiss_df = pd.read_csv(vmiss_fn,sep='\t',header=0)
plt.figure()
ax = vmiss_df['F_MISS'].hist(bins=20)
ax.set_title('distribution of missingness in variant level')
plt.savefig(f'{qc_path}/vmiss.png',dpi=300)

smiss_fn = f'{qc_path}/germ38.smiss'
smiss_df = pd.read_csv(smiss_fn,sep='\t',header=0)
plt.figure()
ax = smiss_df['F_MISS'].hist()
ax.set_title('distribution of missingness in sample level')
plt.savefig(f'{qc_path}/smiss.png',dpi=300)

#------- plot the frequency ---------------------------
freq_fn = f'{qc_path}/freq38.afreq'
freq_df = pd.read_csv(freq_fn, sep='\t',header=0)
print('there are', freq_df.shape[0],'variants in total',
      'among them',freq_df.query('OBS_CT == 2').shape[0],'are singleton')
plt.figure()
ax = freq_df['OBS_CT'].hist(bins=50)
ax.set_title('distribution of observed allele count')
ax.set_xlabel('number of allele count')
ax.set_ylabel('number of SNPs')
ax.set_yscale('log')
plt.savefig(f'{qc_path}/AC_dist.png',dpi=300)

plt.figure()
ax = freq_df.query('OBS_CT >=20')['ALT_FREQS'].hist(bins=50)
ax.set_title('distribution of Allele Frequency')
ax.set_xlabel('Allele frequency')
ax.set_ylabel('number of SNPs')
plt.savefig(f'{qc_path}/AF_dist.png',dpi=300)

#---------- get the variants for downstream analysis ----------
keep_snp_fn = f'{qc_path}/keep_snps38.tsv'
freq_df.query('OBS_CT > 20 and (@af_thred < ALT_FREQS < 1 - @af_thred)')['ID'].to_csv(keep_snp_fn,index=False,header=None)
# extract the variants for downstream analysis
rmSNP_pfile = f'{plink_path}/germ_rmSNP38'
cmd = f'{plink2} --bfile {pfile} --extract {keep_snp_fn} --make-bed --out {rmSNP_pfile}'
run(cmd)





