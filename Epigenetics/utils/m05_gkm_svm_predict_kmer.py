import os, sys

i = sys.argv[1]
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm'
code_path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/lsgkm/src'
kmer_fa = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/gwas/kmer11.fa'

predict_code = f'{code_path}/gkmpredict'


fold_path = f'{path}/model/fold{i}'

model = f'{fold_path}/gkm.model.txt'
test_kmer_predict = f'{fold_path}/predict_kmer.txt'

test_cmd = f'{predict_code} {kmer_fa} {model} {test_kmer_predict}'
print(test_cmd)
os.system(test_cmd)
