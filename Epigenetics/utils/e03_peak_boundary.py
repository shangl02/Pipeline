'''
Created by Shangzhong.Li@pfizer.com on 2022/06/01
'''

import argparse
import os
import json
import sys

code_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(description='predict peak boundary')


parser.add_argument('-p','--path',action='store',dest='path',help='path of origin files of train and test files')

args = parser.parse_args()
path = args.path

# prepare configure
if not os.path.exists(path):
    raise 'path does not exist'


pos_fa = f'{path}/pos.fasta'
neg_fa = f'{path}/neg.fasta'
label_fn = f'{path}/pos.label'
# get number of sequence
with open(pos_fa) as f:
    n = sum(1 for line in f if line.rstrip()) / 2

with open(label_fn,'w') as f:
    f.write('\n'.join(['1.0']*int(n))+'\n')
# create configure file for gradcam
config_fn = f'{path}/object_config.json'
config = {
  "pos_forward_path": pos_fa,
  "index": "all",
  "model_name": "scEpiLock",
  "n_class": 1,
  "model_path": f'{path}/model.pt',
  "target_category": [0],
  "subset": "False",
  "plot": "True",
  "result_path": f'{path}/cam'
}


with open(config_fn,'w') as f:
    json.dump(config, f)
epilock_gradCam_code = f'{code_path}/scEpiLock/grad_cam_module/main.py'

# grad cam module
cmd = f'python {epilock_gradCam_code} {config_fn}'
os.system(cmd)

