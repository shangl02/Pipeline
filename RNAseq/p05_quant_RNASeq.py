'''
Created on 2022/01/21 by Shangzhong.Li@pfizer.com
This run STAR and salmon
'''
import glob, os
from natsort import natsorted
import argparse

parser = argparse.ArgumentParser(description='run RNAseq quant')

parser.add_argument('-i','--input',action='store',dest='input',help='input path')
parser.add_argument('-o','--out',action='store',dest='output',help='out path')
parser.add_argument('-s','--species',action='store',dest='species',help='species')
parser.add_argument('-t','--thread',action='store',dest='thread',help='thread')
parser.add_argument('-p','--paired',action='store',dest='paired',help='paired')
salmon = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/container/salmon_latest.sif'
star = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/container/star_2.7.1a.sif'

args       = parser.parse_args()
fq_path    = args.input
out_path   = args.output
org        = args.species
thread     = args.thread
paired     = args.paired

if org == 'human':
    star_idx = '/hpc/grid/wip_drm_ib/resources/star-index-v27/dna_pa_ercc.hs.h38.ens98/'
    gtf = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/hs/full_ercc.hs.h38.ens98.gtf'
    rna_fa = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/hs/full_ercc.hs.h38.ens98.fa'
elif org == 'mouse':
    star_idx = '/hpc/grid/wip_drm_ib/resources/star-index-v27/dna_pa_ercc.mm.m38.ens98/'
    gtf = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/mm/full_ercc.mm.m38.ens98.gtf'
    rna_fa = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/mm/full_ercc.mm.m38.ens98.fa'
elif org == 'rat':
    star_idx = '/hpc/grid/wip_drm_ib/resources/star-index-v27/dna_pa_ercc.rn.rnor60.ens98'
    gtf = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/rn/full_ercc.rn.rnor60.ens98.gtf'
    rna_fa = '/hpc/grid/wip_drm_ib/resources/ensembl/release-98/derived/rn/full_ercc.rn.rnor60.ens98.fa'

# path = '/hpc/grid/wip_drm_targetsciences/projects/p042_IL26'
# fq_path = '/hpc/grid/dsrd-globalpath/220114EPHRQSLiver/'
bam_path = f'{out_path}/bam'
os.makedirs(bam_path, exist_ok=True)
count_path = f'{out_path}/count'
os.makedirs(count_path,exist_ok=True)


def run_single_end(fq_files):
    for fq in fq_files:
        prefix = fq.split('/')[-1].split('.')[0]
        out = f'{bam_path}/{prefix}_'
        out_bam = f'{out}Aligned.toTranscriptome.out.bam'

        salmon_out = '_'.join(out_bam.split('/')[-1].split('_')[:2])
        salmon_out = f'{count_path}/{salmon_out}'
        if not os.path.exists(out_bam):
            cmd = f'singularity run -B /:/media {star} STAR --genomeDir /media/{star_idx} \
                   --sjdbGTFfile /media/{gtf} --readFilesIn /media/{fq} --readFilesCommand zcat \
                   --runThreadN {thread} --limitBAMsortRAM 53687091200 --alignEndsType EndToEnd \
                   --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate \
                   --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.05 \
                   --outFilterScoreMinOverLread 0.90 --outFilterMatchNminOverLread 0.90 \
                   --alignIntronMax 1000000 --outFileNamePrefix /media/{out} && \
                   singularity run -B /:/media {salmon} salmon quant -t /media/{rna_fa}\
           -l ISR -a /media{out_bam} -p {thread} -o /media/{salmon_out} --incompatPrior 0 \
           -g /media/{gtf}'
            log = f'{salmon_out}.log'
            bsub = f'bsub -o {log} -q medium -M 63687091200 -n {thread} \
                    -R \"select[hname!=amrndhl1296]\" \
                     -R \"span[ptile={thread}]\" \"{cmd}\"'
            os.system(bsub)
            # print(bsub)


def run_paired_end(fst_files, snd_files):
    for fst,snd in zip(fst_files,snd_files):
        prefix = fst.split('/')[-1].split('.')[0]
        out = f'{bam_path}/{prefix}_'
        out_bam = f'{out}Aligned.toTranscriptome.out.bam'

        salmon_out = '_'.join(out_bam.split('/')[-1].split('_')[:2])
        salmon_out = f'{count_path}/{salmon_out}'
        if not os.path.exists(out_bam):
            cmd = f'singularity run -B /:/media {star} STAR --genomeDir /media/{star_idx} \
                   --sjdbGTFfile /media/{gtf} \
                   --readFilesIn /media/{fst} /media/{snd} \
                   --readFilesCommand zcat \
                   --runThreadN {thread} --limitBAMsortRAM 53687091200 --alignEndsType EndToEnd \
                   --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate \
                   --alignSJDBoverhangMin 1 --outFilterMismatchNoverLmax 0.05 \
                   --outFilterScoreMinOverLread 0.90 --outFilterMatchNminOverLread 0.90 \
                   --alignIntronMax 1000000 --outFileNamePrefix /media/{out} && \
                   singularity run -B /:/media {salmon} salmon quant -t /media/{rna_fa}\
           -l ISR -a /media{out_bam} -p {thread} -o /media/{salmon_out} --incompatPrior 0 \
           -g /media/{gtf}'
            log = f'{salmon_out}.log'
            bsub = f'bsub -o {log} -q medium -M 63687091200 -n {thread} \
                    -R \"select[hname!=amrndhl1296]\" \
                     -R \"span[ptile={thread}]\" \"{cmd}\"'
            os.system(bsub)
            # print(bsub)


if paired == 'Y': 
    fq_files = natsorted(glob.glob(fq_path + '/*.fastq.gz'))
    fst_files = [fq_files[i] for i in range(0,len(fq_files),2)]
    snd_files = [fq_files[i] for i in range(1,len(fq_files),2)]
    run_paired_end(fst_files, snd_files)
else:
    fq_files = natsorted(glob.glob(fq_path + '/*.fastq.gz'))
    run_single_end(fq_files)