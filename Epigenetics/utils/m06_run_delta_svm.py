import os, sys

i = sys.argv[1]
path = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm'
code_path = '/home/lis262/Code/Scripts/HPC_Scripts/epigenetics/deltasvm_script'


delta_code = f'{code_path}/deltasvm.pl'
fold_path = f'{path}/model/fold{i}'

risk_fa = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/gwas/Lung_GWAS_20210629.risk_50.fa'
norm_fa = '/hpc/grid/wip_drm_targetsciences/users/shangzhong/epigenetics/gkm_svm/gwas/Lung_GWAS_20210629.norm_50.fa'

kmer_predict = f'{fold_path}/predict_kmer.txt'
delta_out = f'{fold_path}/predict_delta_svm.txt'

delta_svm_cmd = f'perl {delta_code} {norm_fa} {risk_fa} {kmer_predict} {delta_out}'
# print(delta_svm_cmd)
os.system(delta_svm_cmd)
