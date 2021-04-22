'''
Created on 2021/01/13 by Shangzhong.Li@pfizer.com
this file run repeat expansion length detection using nanopore

'''
from yaml import load, Loader
import os
import sys
import pysam
import glob
import pandas as pd
from matplotlib import pyplot as plt
plt.style.use('ggplot')
from intervaltree import Interval, IntervalTree
import numpy as np
from pybedtools import BedTool
import time


start = time.time()
param_fn = sys.argv[1]
with open(param_fn, 'r') as in_f:
    param = load(in_f,Loader=Loader)

disease = '/media' + param['disease']
target = '/media' + param['target']
repeat = '/media' + param['repeat_info']
fast5_path = param['fast5_path']
fastq_path = param['fastq_path']
out_path = '/media' + param['out_path']
ref_fa = '/media' + param['ref_fa']

if not os.path.exists(out_path):
    os.mkdir(out_path)


raw_fq_path = fastq_path
fq_file = out_path + '/reads.fq.gz'

os.chdir(out_path)
#============= 0. merge fastq files ==================
raw_fq_files = sorted(glob.glob(raw_fq_path + '/*.fastq'))
cmd = 'ls {fq_p}/*.fastq | xargs cat | gzip - > {o}'.format(fq_p=raw_fq_path,o=fq_file)
os.system(cmd)

fq_path = os.path.dirname(fq_file)
fq_name = fq_file.split('/')[-1]
ref_path = os.path.dirname(ref_fa)
ref_name = ref_fa.split('/')[-1]
# #============== 1. run QC using nanopack
try:
    cmd = ('NanoPlot --fastq {fq} \
     -o f01_nanoplot').format(fq=fq_name)
    print('Step1: run QC using nanopack')
    print(cmd)
    os.system(cmd)
    print('Finish Step1: run QC using nanopack\n\n')
except:
    assert False,'nanoplot failed'

# #============== 2. map the reads to the reference genome
map_path = out_path + '/f02_minimap'
if not os.path.exists(map_path):
    os.mkdir(map_path)
bam = map_path + '/minimap.bam'
cmd = ('minimap2 -ax map-ont {ref} {fq} | samtools sort -T {p}/pre -o {bam} && samtools index {bam}').format(p=fq_path,ref=ref_fa,fq=fq_file,bam=bam)

print('Step2: map reads to the reference genome')
print(cmd)
os.system(cmd)
print('Finish Step2: mapping reads to the reference genome\n\n')

#============== 3. count reads mapping to each chromosome
fig_path = out_path + '/f03_plot'
if not os.path.exists(fig_path):
    os.mkdir(fig_path)
# count reads mapping to each chromosome
def chrome_count(bam):
    chrom_n_dict = {}
    chroms = ['chr'+ str(i) for i in range(1,23)] + ['chrX','chrY']
    for c in chroms:
        chrom_n_dict[c] = 0
    obj = pysam.AlignmentFile(bam, "rb",check_sq=False)
    for read in obj:
        try:
            ref = read.reference_name
            chrom_n_dict[ref] += 1
        except:
            pass
    return chrom_n_dict

def plot_chrom_count(bam, fig_path):
    
    chrom_n_dict = chrome_count(bam) 
    chrom_count_df = pd.DataFrame.from_dict(chrom_n_dict,orient='index',
            columns=['count']).sort_values('count',ascending=False)
    ax = chrom_count_df['count'].plot(kind='bar')
    ax.set_ylabel('count')
    _ = ax.set_title('read count at each chromosome')
    plt.savefig(fig_path + '/chrom_read_count.png')
    
print('Step3: count reads mapping to each chromosome')
plot_chrom_count(bam, fig_path)
print('Finish Step3: counting reads\n\n')

#=============== 4. check coverage at each covering loci ==================
def cov_at_loci(bam, cov_fn):
    bed = BedTool(bam)
    cov_df = bed.genome_coverage(bg=True).to_dataframe()

    # get coverage
    inter = {}
    for idx, row in cov_df.iterrows():
        chrom,s,e,n = row
        s = int(s); e = int(e)
        if chrom not in inter:
            inter[chrom] = [[s,e]]
        else:
            if s - inter[chrom][-1][1] < 1000:
                inter[chrom][-1][1] = e
            else:
                inter[chrom].append([s,e])
    # prepare coverage list
    with open(cov_fn,'w') as out_f:
        for k,v in inter.items(): # k is chrom, v is list of region pos
            for c in v:
                s,e = c
                cov = cov_df.query('chrom==@k and start >=@s and end <=@e')['name'].mean()
                out_f.write('\t'.join([k,str(s),str(e),str(cov)])+'\n')

    cov_df = pd.read_csv(cov_fn,sep='\t',header=0,names=['chr','s','e','cov'])
    cov_df = cov_df.sort_values('cov',ascending=False)
    cov_df.to_csv(cov_fn,sep='\t',index=False)

def plot_loc_cov(bam, target, fig_path, disease):
    cov_fn = fig_path + '/mean_cov.txt'
    cov_at_loci(bam, cov_fn)
    # plot coverage
    cov_df = pd.read_csv(cov_fn,sep='\t',header=0,names=['chr','s','e','cov'])

    t_chr,t_s,t_e = target.split(':')
    t_s = int(t_s)
    t_e = int(t_e)
    target_cov = cov_df.query('chr==@t_chr and s >=@t_s and e <=@t_e')['cov'].mean()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(np.array(cov_df['cov'].tolist()),bins=50)
    ax.set_yscale('log')
    ax.set_xlabel('average coverage')
    ax.set_ylabel('number of regions')
    ax.set_title('distribution of coverage for '+disease)
    plt.axvline(x=float(target_cov), color='b',linestyle='-')
    plt.savefig(fig_path + '/mean_cov_in_each_region.png')


print('Step4: check coverage at each covering loci')
plot_loc_cov(bam, target, fig_path, disease)
print('Finish Step4: check cov at each cover loci\n\n')

#================== 5. detect repeat length ===================
# 1. index the reads
fast5_idx_path = out_path + '/f04_repeat_count'
if not os.path.exists(fast5_idx_path):
    os.mkdir(fast5_idx_path)

cmd = ('python /opt/STRique/scripts/STRique.py index {p} \
       --out_prefix {p}/ > {f}/reads.fofn'\
       .format(p=fast5_path,f=fast5_idx_path))

print('Step 5.1 index the raw reads')
print(cmd)
os.system(cmd)


# 2. calculate the length of the repeats
repeat_fn = repeat.split('/')[-1]

cmd = ('samtools view {bam} |  \
       python /opt/STRique/scripts/STRique.py count \
       --out {out_path}/f04_repeat_count/hg38.strique.tsv \
       {out_path}/f04_repeat_count/reads.fofn \
       /opt/STRique/models/r9_4_450bps.model \
       {config} \
       ').format(bam=bam,out_path=out_path,config=repeat)

print('Step 5.2 count the repeats')
print(cmd)
os.system(cmd)
end = time.time()
print(end - start)