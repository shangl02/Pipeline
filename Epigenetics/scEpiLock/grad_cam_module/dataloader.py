import numpy as np
import torch
import pandas as pd

class DataLoader:

    def __init__(self, fasta_path, index, subset):
        """

        :param fasta_path: path to the fasta file
        :param index: index of the peak to be evaluated
        """
        self.pos_fasta_path = fasta_path
        self.index = index
        self.subset = subset

    # return array of forward seq and reverse complement seq
    def get_seq(self):
        if self.index == "all":
            data = pd.read_csv(self.pos_fasta_path, sep=">chr*",
                              header=None, engine='python').values[1::2][:, 0]
        else:
            data = pd.read_csv(self.pos_fasta_path, sep=">chr*",
                              header=None, engine='python').values[1::2][:, 0][self.index]
            # change to np.array
            data = np.array([data])

        if self.subset == "True":
            size = int(np.floor(data.shape[0] * 0.01))
            np.random.seed(202101190)
            index = np.random.choice(data.shape[0], size=size, replace=False)
            data = data[index]


        # print(seq)
        # reverse complement
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        temp = []

        for seq in data:
            complement_seq = ''
            for base in seq:
                if base not in complement.keys():base='N'
                complement_seq = complement[base] + complement_seq
            temp.append(complement_seq)

        temp = np.array(temp, dtype=object)

        print("----------------- check the data and reverse complement shape -----------------") # TODO delete this print (YG)

        print(data.shape) #(114538,)
        print(temp.shape) #(114538,)
        return data, temp

    #TODO: update this one for array as input
    def seq_one_hot_tensor(self, data):
        temp_list = []
        for seq in data:
            temp = np.zeros((len(seq), 4))
            row_index = 0
            for base in seq:
                if base == 'A':
                    temp[row_index, 0] = 1
                elif base == 'T':
                    temp[row_index, 1] = 1
                elif base == 'G':
                    temp[row_index, 2] = 1
                elif base == 'C':
                    temp[row_index, 3] = 1
                row_index += 1
            temp_list.append(temp)

        data_tensor = torch.tensor(temp_list).float().permute(0, 2, 1)
        print("----------------- check the data tensor shape -----------------")
        print(data_tensor.shape)

        return data_tensor
