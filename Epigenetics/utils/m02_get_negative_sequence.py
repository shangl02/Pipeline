import os, sys
import pysam
import pandas as pd
import argparse
from seqdataloader import *
from seqdataloader.labelgen import *


parser = argparse.ArgumentParser(description='get negative seqeunce from peak file')


parser.add_argument('-r','--ref_fa',action='store',dest='ref_fa',help='reference genome fa file')
parser.add_argument('-p','--peak_fn',action='store',dest='peak_fn',help='peak file name')
parser.add_argument('-c','--chr_size',action='store',dest='chr_size_fn',help='file inidicate chromosome size')
parser.add_argument('--path',action='store',dest='path',help='path to save model')
parser.add_argument('--neg',action='store',dest='neg_code',help='code to extract negative code')
parser.add_argument('-f','--fold',action='store',dest='fold',help='fold of cross validation',type=int,default=1)

args = parser.parse_args()

ref_fa = args.ref_fa
peak_fn = args.peak_fn
model_path = args.path
neg_code = args.neg_code
chrom_size_fn = args.chr_size_fn
task_fn = f'{model_path}/task.txt'
fold = args.fold

def seqloader(model_path, peak_fn, task_fn, chrom_size_fn, i):
    test_chroms = [['chr1'],['chr2','chr19'],['chr3','chr20'],['chr6','chr13','chr22'],
                 ['chr5','chr16','chrY'],['chr4','chr15','chr21'],['chr7','chr14','chr18'],
                 ['chr11','chr17','chrX'],['chr9','chr12'],['chr8','chr10']]
    fold_path = f'{model_path}/fold{i}'
    os.makedirs(fold_path,exist_ok=True)
    tmp_path = f'{model_path}/tmp'
    os.makedirs(tmp_path,exist_ok=True)
    train_params={
        'task_list':task_fn,
        'outf':f'{model_path}/fold{i}/train.inputs.bed.gz',
        'output_type':'gzip',
        'chrom_sizes':chrom_size_fn,
        'chroms_to_exclude':test_chroms[i],
        "store_positives_only":False,
        'bin_stride':200,
        'left_flank':300,
        'right_flank':300,
        'bin_size':400,
        'chrom_threads':6,
        'task_threads':4,
        'allow_ambiguous':True,
        'labeling_approach':'peak_summit_in_bin_classification',
        'temp_dir':tmp_path
        }
    
    test_params={
        'task_list':task_fn,
        'outf':f'{model_path}/fold{i}/test.inputs.bed.gz',
        'output_type':'gzip',
        'chrom_sizes':chrom_size_fn,
        'chroms_to_keep':test_chroms[i],
        "store_positives_only":False,
        'bin_stride':150,
        'left_flank':300,
        'right_flank':300,
        'bin_size':400,
        'chrom_threads':6,
        'task_threads':4,
        'allow_ambiguous':True,
        'labeling_approach':'peak_summit_in_bin_classification',
        'temp_dir':tmp_path
        }
    # create task_fn
    cell = model_path.split('/')[-1]
    with open(task_fn,'w') as f:
        lines = [['task', 'narrowPeak'],[cell,f'{model_path}/peak.bed']]
        for line in lines:
            f.write('\t'.join(line) + '\n')
    genomewide_labels(train_params)
    genomewide_labels(test_params)



def get_seqs(bed, fasta, ref_fa):
    '''

    '''
    ref = pysam.FastaFile(ref_fa)

    df_bed = pd.read_csv(bed, sep='\t', header=None)
    with open(fasta,'w') as fa_file:
        counter = 0
        for index,row in df_bed.iterrows():
            seq=ref.fetch(row[0],int(row[1]),int(row[2]))
            head = f'>{row[0]}_{row[1]}_{row[2]}'
            fa_file.write(head + '\n')
            fa_file.write(seq.upper() + '\n')
            counter += 1


def get_neg_seq(model_path, ref_fa, i, dtype):
    '''
    get negative sequence
    * dtype: train or test
    '''
    ratio = 1.0 # neg to pos
    fold_path = f'{model_path}/fold{i}'
    # 1.1 prepare the file for train/test
    bed = f'{fold_path}/{dtype}_pos.bed'
    all_bed = f'{fold_path}/{dtype}_pos_all.bed'
    cmd = f'cut -f 1,2,3 {bed} | sed -e "s/$/\t1.0/" > {all_bed}'
    # print(cmd)
    os.system(cmd)
    # 1.2 add negative
    data_input = f'{fold_path}/{dtype}.inputs.bed.gz'
    cmd = f'zcat {data_input} | tail -n +2 | awk \'NF==4\' | \
                awk \'$4 != 1.0 \' >> {all_bed}'
    os.system(cmd)
    # print(cmd)
    # 1.3 get gc frequency
    gc_freq = f'{fold_path}/{dtype}_pos_gc.txt'
    cmd = f'python {neg_code} -b {all_bed} --ratio_neg_to_pos {ratio} \
                -o {gc_freq} --ref_fasta {ref_fa} --gc'
    print(cmd)
    os.system(cmd)
    # 1.4 get get dinnuc_freq
    final_bed = f'{fold_path}/{dtype}_final.bed'
    cmd = f'python {neg_code} -b {all_bed} --ratio_neg_to_pos {ratio} \
             -o {final_bed} --ref_fasta {ref_fa} \
             --gc --dinuc_freq {gc_freq}'
    print(cmd)
    os.system(cmd)
    # 1.5 
    with open(bed) as f:
        pos_len =  sum(1 for line in f if line.rstrip())
    with open(f'{final_bed}.0') as f:
        total_len =  sum(1 for line in f if line.rstrip())
    neg_len = total_len - pos_len
    final_pos_bed = f'{fold_path}/{dtype}_final_pos.bed'
    cmd = (f'head -n {pos_len} {final_bed}.0 > {final_pos_bed}')
    os.system(cmd)
    final_neg_bed = f'{fold_path}/{dtype}_final_neg.bed'    
    cmd = (f'tail -n {neg_len} {final_bed}.0 > {final_neg_bed}')
    os.system(cmd)

    get_seqs(final_pos_bed, f'{fold_path}/{dtype}_final_pos.fasta',ref_fa)
    get_seqs(final_neg_bed, f'{fold_path}/{dtype}_final_neg.fasta',ref_fa)


# loop for each fold
for i in range(fold):
    seqloader(model_path, peak_fn, task_fn, chrom_size_fn, i)
    get_neg_seq(model_path, ref_fa, i, 'train')
    get_neg_seq(model_path, ref_fa, i, 'test')

