import os
import argparse

parser = argparse.ArgumentParser(description='run matrixQTL')

parser.add_argument('-p','--path',action='store',dest='path',help='work path of input files for matrixQTL')


args = parser.parse_args()
work_path = args.path

matrixQTL_path = f'{work_path}/matrixQTL'
matrixQTL = '~/Code/Pipeline/TME/QTL_pipeline/utils/m06_matrixQTL.R'
matrixQTL = f'{os.path.dirname(os.path.realpath(__file__))}/m06_matrixQTL.R'

# 1. run matrixQTL
log = f'{work_path}/qtl_log.txt'
out_fn = 'qtl.txt'
cmd = f'bsub -R \"select[hname!=amrndhl1296]\" \
    -M 100457280 -q medium -o {log} \"Rscript {matrixQTL} \
      {matrixQTL_path} {out_fn} \
     && gzip {matrixQTL_path}/{out_fn}\"'
print(cmd)
# os.system(cmd)
