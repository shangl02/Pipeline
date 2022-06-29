import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils.class_weight import compute_class_weight
from utils.utils import Utils

class PreProcessor():
    # TODO: change the encode_path naming to neg_path
    def __init__(self, pos_fasta_path, encode_path, neg_path, label_path, encode_n, neg_n,
                 cell_cluster, subset):
        self.pos_fasta_path = pos_fasta_path
        self.encode_path = encode_path
        self.neg_path = neg_path
        self.label_path = label_path
        self.encode_n = encode_n
        self.neg_n = neg_n
        self.cell_cluster = cell_cluster
        self.subset = subset

    def concate_data(self):
        pos_fasta = pd.read_csv(self.pos_fasta_path,sep=">chr*",header=None, engine='python').values[1::2][:,0] # (n,) array, each element is string, dtype=object
        encode_fasta = pd.read_csv(self.encode_path,sep=">chr*",header=None, engine='python').values[1::2][:,0]
        neg_fasta = pd.read_csv(self.neg_path,sep=">chr*",header=None, engine='python').values[1::2][:,0]
        label = pd.read_csv(self.label_path, delimiter = ",", header=None)

        Utils.print_separator("finish pd.read_csv")
        Utils.print_separator("original shape")
        print("pos_size") # TODO Delte these print -YG
        print(pos_fasta.shape) #(114538,) / brain 221062
        print("encode_size")
        print(encode_fasta.shape) #(1127227,)/ brain 489777
        print("neg_size")
        print(neg_fasta.shape) #/ brain 1098300
        print("label")
        print(label.shape) #(114538, 1) /brain 221062

        np.random.seed(202101190)
        encode_index = np.random.choice(encode_fasta.shape[0], size=self.encode_n, replace=False)
        encode_fasta = encode_fasta[encode_index]

        np.random.seed(202101190)
        neg_index = np.random.choice(neg_fasta.shape[0], size=self.neg_n, replace=False)
        neg_fasta = neg_fasta[neg_index]


        data = np.concatenate([pos_fasta, encode_fasta, neg_fasta])
        label = label.to_numpy()

        label = label[:, self.cell_cluster] ## select the cluster for label

        encode_label = np.zeros((self.encode_n, len(self.cell_cluster)))
        neg_label = np.zeros((self.neg_n, len(self.cell_cluster)))

        label = np.concatenate([label, encode_label, neg_label])


        print("-----------shape right after encode/neg concate------------")
        print(data.shape) #(124538,) /brain 241062
        print(label.shape) #(124538, 1) / brain (241062, 6)

        if self.subset == "True":
            data, label = self._subset(data, label)
        print("-----------------shape after subset -------------")
        print(data.shape)
        print(label.shape)


        print("-----------original pos/neg ratio------------")
        # brain 0.24385499314593737
        print(np.count_nonzero(label)/(label.shape[0]*label.shape[1]-np.count_nonzero(label)))
        pos_weight = []
        for i in range(0,label.shape[1]):
            num_pos = np.count_nonzero(label[:, i])
            num_neg = label.shape[0] - num_pos
            pos_weight.append(float(num_neg)/num_pos)
        print(pos_weight)


        return data, label, pos_weight

    def split_train_test(self, data, label, test_size = 0.1):
        data_train_temp, data_test, label_train_temp, label_test = train_test_split(data, label, test_size=test_size, random_state=12)
        data_train, data_eval, label_train, label_eval = train_test_split(data_train_temp, label_train_temp, test_size=test_size, random_state=12)

        print("-----------train size concate------------")
        print(data_train.shape) #(173775,)
        print(label_train.shape) #(173775, 1)

        print("-----------eval size concate------------")
        print(data_eval.shape) #(19309,)
        print(label_eval.shape) #(19309, 1)

        print("-----------test size concate------------")
        print(data_test.shape) #(21454,)
        print(label_test.shape) #(21454,) no encode: (11454, 8)

        test_size = data_test.shape[0]

        return data_train, data_eval, data_test, label_train, label_eval, label_test, test_size

    def _subset(self, data, label):
        size = int(np.floor(data.shape[0] * 0.2))
        np.random.seed(202101190)
        index = np.random.choice(data.shape[0], size=size, replace=False)

        data = data[index]
        label = label[index, :]

        return data, label

    # TODO 04/12 add the balance option



class PreProcessorTrans():
    def __init__(self, pos_forward_path, neg_forward_path, balance, subset):
        self.pos_forward_path = pos_forward_path
        self.neg_forward_path = neg_forward_path
        self.balance = balance
        self.subset = subset

    def concate_data(self):
        pos_fasta = pd.read_csv(self.pos_forward_path, sep=">chr*",
                                header=None, engine='python').values[1::2][:, 0]  # (n,) array, each element is string, dtype=object
        neg_fasta = pd.read_csv(self.neg_forward_path, sep=">chr*",
                                header=None, engine='python').values[1::2][:, 0]

        print("-----------finish pd read_csv------------")
        print("-----------original shape------------")
        print("positive")
        print(pos_fasta.shape)
        print("negative")
        print(neg_fasta.shape)

        if self.balance == "True":
            pos_fasta, neg_fasta = self._balance(pos_fasta, neg_fasta)

        print("-----------positve/negative ratio after balance------------")
        print(pos_fasta.shape[0]/neg_fasta.shape[0])

        data = np.concatenate([pos_fasta, neg_fasta])

        n_pos = pos_fasta.shape[0]
        n_neg = neg_fasta.shape[0]

        pos_label = np.ones(n_pos)
        neg_label = np.zeros(n_neg)

        label = np.concatenate([pos_label, neg_label]).reshape((n_pos+n_neg), 1) # nx1

        print("-----------shape right after join pos and neg fasta------------")
        print(data.shape)  # (1681932,)
        print(label.shape)  # (1681932, 1)

        if self.subset == "True":
            data, label = self._subset(data, label)

        print("-----------------shape after subset -------------")
        print(data.shape)
        print(label.shape)

        print("-----------pos/neg ratio------------")
        print(np.count_nonzero(label) / (label.shape[0] * label.shape[1] - np.count_nonzero(label)))
        pos_weight = []
        for i in range(label.shape[1]):
            num_pos = np.count_nonzero(label[:, i])
            num_neg = label.shape[0] - num_pos
            pos_weight.append(float(num_neg) / num_pos)
        print(pos_weight)

        return data, label, pos_weight

    def split_train_test(self, data, label, test_size = 0.1):
        data_train_temp, data_test, label_train_temp, label_test = train_test_split(data, label, test_size=test_size, random_state=12)
        data_train, data_eval, label_train, label_eval = train_test_split(data_train_temp, label_train_temp, test_size=test_size, random_state=12)
        test_size = data_test.shape[0]

        return data_train, data_eval, data_test, label_train, label_eval, label_test, test_size

    def _subset(self, data, label):
        size = int(np.floor(data.shape[0] * 0.1))
        np.random.seed(202101190)
        index = np.random.choice(data.shape[0], size=size, replace=False)

        data = data[index]
        label = label[index, :]

        return data, label

    def _balance(self, pos, neg):
        pos_size = pos.shape[0]
        neg_size = neg.shape[0]

        if pos_size > neg_size:
            ratio = neg_size/pos_size
            size = int(np.floor(pos.shape[0] * ratio))
            np.random.seed(202101190)
            index = np.random.choice(pos.shape[0], size=size, replace=False)
            pos = pos[index]
        else:
            ratio = pos_size / neg_size
            size = int(np.floor(neg.shape[0] * ratio))
            np.random.seed(202101190)
            index = np.random.choice(neg.shape[0], size=size, replace=False)
            neg = neg[index]

        return pos, neg