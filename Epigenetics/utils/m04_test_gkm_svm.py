import os, sys

i = sys.argv[1]
path = '/hpc/grid/wip_drm_targetsciences/projects/EpiMiner/atac/p01_gkm_svm'
code_path = '~/Code/Scripts/HPC_Scripts/epigenetics/lsgkm/src'


predict_code = f'{code_path}/gkmpredict'


fold_path = f'{path}/ACardiomyocyte/fold{i}'

model = f'{fold_path}/gkm.model.txt'
test_pos_fn = f'{fold_path}/test_final_pos.fasta'
test_neg_fn = f'{fold_path}/test_final_neg.fasta'
test_pos_predict = f'{fold_path}/predict_test_pos.txt'
test_neg_predict = f'{fold_path}/predict_test_neg.txt'

test_cmd = f'{predict_code} {test_pos_fn} {model} {test_pos_predict}'
print(test_cmd)
os.system(test_cmd)

test_cmd = f'{predict_code} {test_neg_fn} {model} {test_neg_predict}'
print(test_cmd)
os.system(test_cmd)

# print(train_cmd)