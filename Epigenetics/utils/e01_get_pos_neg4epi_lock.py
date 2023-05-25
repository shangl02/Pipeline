import argparse
import os

parser = argparse.ArgumentParser(description='get positive and negative sequence')

parser.add_argument('-p','--path',action='store',dest='path',help='path of origin train and test files')

args = parser.parse_args()
path = args.path


# merge positive and negative
train_pos_fn = f'{path}/train_final_pos.fasta'
test_pos_fn = f'{path}/test_final_pos.fasta'
cmd = f'cat {train_pos_fn} {test_pos_fn} > {path}/pos.fasta'
os.system(cmd)

train_neg_fn = f'{path}/train_final_neg.fasta'
test_neg_fn = f'{path}/test_final_neg.fasta'
cmd = f'cat {train_neg_fn} {test_neg_fn} > {path}/neg.fasta'
os.system(cmd)

train_pos_bed = f'{path}/train_final_pos.bed'
test_pos_bed = f'{path}/test_final_pos.bed'
cmd = f'cat {train_pos_bed} {test_pos_bed} > {path}/pos.bed'
os.system(cmd)

train_neg_bed = f'{path}/train_final_neg.bed'
test_neg_bed = f'{path}/test_final_neg.bed'
cmd = f'cat {train_neg_bed} {test_neg_bed} > {path}/neg.bed'
os.system(cmd)

# merge label
train_pos_bed = f'{path}/train_final_pos.bed'
test_pos_bed = f'{path}/test_final_pos.bed'
cmd = f'cut -f 4 {train_pos_bed} > {path}/pos.label'
os.system(cmd)
cmd = f'cut -f 4 {test_pos_bed} >> {path}/pos.label'
os.system(cmd)
train_neg_bed = f'{path}/train_final_neg.bed'
test_neg_bed = f'{path}/test_final_neg.bed'
cmd = f'cut -f 4 {train_neg_bed} > {path}/neg.label'
os.system(cmd)
cmd = f'cut -f 4 {test_neg_bed} >> {path}/neg.label'
os.system(cmd)


