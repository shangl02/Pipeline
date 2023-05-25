import os, sys
import argparse

code_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(description='train gkm svm')

parser.add_argument('-p','--data_path',action='store',dest='data',help='data path with train and test data')


args = parser.parse_args()
data_path = args.data

lsgkm_path = f'{code_path}/lsgkm/src'
train_code = f'{lsgkm_path}/gkmtrain'
predict_code = f'{lsgkm_path}/gkmpredict'

t = 16 # can only be 1, 4,or 16
def train(data_path, train_code):
    for i in range(10):
        fold_path = f'{data_path}/fold{i}'

        train_pos_fn = f'{fold_path}/train_final_pos.fasta'
        train_neg_fn = f'{fold_path}/train_final_neg.fasta'

        outpre = f'{fold_path}/gkm'
        train_cmd = f'{train_code} {train_pos_fn} {train_neg_fn} {outpre} -T {t} -m 30000'
        log = f'{fold_path}/log_train.txt'
        bsub = f'bsub -n {t} -R \"span[ptile={t}]\" -M 35457280 -o {log} \
                 -q long \"{train_cmd}\"'
        os.system(bsub)
        # print(train_cmd)


if __name__ == "__main__":
    train(data_path, train_code)