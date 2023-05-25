import os
from Bio import SeqIO
import pandas as pd
import glob
import argparse

parser = argparse.ArgumentParser(description='get positive seqeunce from peak file')

parser.add_argument('-r','--ref_fa',action='store',dest='ref_fa',help='reference genome fa file')
parser.add_argument('-p','--peak_fn',action='store',dest='peak',help='peak file')
parser.add_argument('-d','--dist',action='store',dest='dist',help='distance up and downstream of peak center',default=500)
parser.add_argument('--path',action='store',dest='path',help='path to save model')
parser.add_argument('-f','--fold',action='store',dest='fold',help='number of total folds', default=10)
parser.add_argument('--header',action='store',dest='header',help='if file has header', default='no')

args = parser.parse_args()
ref_fa = args.ref_fa
peak_fn = args.peak
distance = int(args.dist)
model_path = args.path
fold = int(args.fold)
header = args.header


def get_peak_fa(df, fa_dict, out_fa, out_bed):
    '''get sequence of the peaks'''
    chroms = [f'chr{c}' for c in list(range(1,23))] + ['chrX','chrY']
    
    with open(out_fa,'w') as fa, open(out_bed,'w') as bed:
        for c in chroms:
            chrom_seq = fa_dict[c].seq
            sub_df = df[df[0].values == c]
            for idx, row in sub_df.iterrows():
                start = row['lo'] - 1
                end = row['hi']
                sequence = chrom_seq[start:end]
                if 'N' in sequence: continue
                fa.write(f'>{c}_{start}_{end}\n{sequence}\n')
                # col_len = len(row)
                bed.write('\t'.join([str(row[0]),str(row['lo']),str(row['hi'])] + \
                                   [str(i) for i in row[3:9]] + ['500'])+'\n')
                
                
def get_positive_fa(fa, bed, model_path, fold):
    folds = [['chr1'],['chr2','chr19'],['chr3','chr20'],['chr6','chr13','chr22'],
             ['chr5','chr16','chrY'],['chr4','chr15','chr21'],['chr7','chr14','chr18'],
             ['chr11','chr17','chrX'],['chr9','chr12'],['chr8','chr10']]
    n_fold = min(len(folds), fold)
    for i in range(n_fold):
        fold_path = f'{model_path}/fold{i}'
        os.makedirs(fold_path,exist_ok=True)
        train_fa = f'{fold_path}/train_pos.fa'
        test_fa = f'{fold_path}/test_pos.fa'
        train_bed = f'{fold_path}/train_pos.bed'
        test_bed = f'{fold_path}/test_pos.bed'
        test_chrom = folds[i]
        # output bed for train and test
        df = pd.read_csv(bed,sep='\t',header=None)
        train_df = df[~df[0].isin(test_chrom)]
        test_df = df[df[0].isin(test_chrom)]
        train_df.to_csv(train_bed,sep='\t',index=False,header=False)
        test_df.to_csv(test_bed,sep='\t',index=False,header=False)
        # output to fa file
        with open(train_fa, 'w') as train_f, open(test_fa, 'w') as test_f:
            for record in SeqIO.parse(fa,'fasta'):
                if record.id.split('_')[0] in test_chrom:
                    SeqIO.write(record, test_f,'fasta')
                else:
                    SeqIO.write(record,train_f,'fasta')
        
if __name__ == "__main__":
    # get peak center and extend in both direction
    if header == 'no':
        peak_df = pd.read_csv(peak_fn, sep='\t',header=None,compression='gzip')
    elif header == 'yes':
        peak_df = pd.read_csv(peak_fn, sep='\t',header=0,compression='gzip')
        peak_df.columns = range(peak_df.shape[1])
    peak_df['center'] = ((peak_df[1] + peak_df[2]) / 2).astype('int')
    peak_df['lo'] = peak_df['center'] - distance
    peak_df['hi'] = peak_df['center'] + distance


    index = SeqIO.index(ref_fa, 'fasta')
    os.makedirs(model_path, exist_ok=True)


    peak_fa = f'{model_path}/peak.fa'
    peak_bed = f'{model_path}/peak.bed'
    get_peak_fa(peak_df, index, peak_fa, peak_bed)
    get_positive_fa(peak_fa, peak_bed, model_path, fold)
