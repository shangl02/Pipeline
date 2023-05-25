'''
Created by Shangzhong.Li@pfizer.com on 2022/06/01
'''

import argparse
import os
import json
import sys

code_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(description='train the model')


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

# with open(label_fn,'w') as f:
#     f.write('\n'.join(['1.0']*int(n))+'\n')
# create configure file for training
config_fn = f'{path}/learn_config.json'
config = {
  "pos_forward_path": pos_fa,
  "encode_path": f'{code_path}/scEpiLock/encode.fa',
  "neg_path": neg_fa,
  "label_path": label_fn,
  "encode_n": 0,
  "neg_n": int(n),
  "model_name": "scEpiLock",
  "load_trans_model": "False",
  "num_epochs": 25,
  "batch_size": 64,
  "learning_rate": 0.0005,
  "weight_decay": 0,
  "use_pos_weight": "False",
  "cell_cluster": "all",
  "n_class": 111,
  "n_class_trans":111,
  "model_path": f'{path}/model.pt',
  "output_evaluation_data_path": f'{path}/',
  "trans_model_path": '',
  "subset": "False",
  "balance": "False"
}
with open(config_fn,'w') as f:
    json.dump(config, f)
epilock_train_code = f'{code_path}/scEpiLock/deep_learning_module/main.py'
#===== trian model ==========
# train the model
cmd = f'python {epilock_train_code} -e train -s main -c {config_fn}'
# os.system(cmd)
# test the model
cmd = f'python {epilock_train_code} -e test -s main -c {config_fn}'
os.system(cmd)


