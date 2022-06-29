import numpy as np
import torch
from torch.utils.data import Dataset
from datetime import datetime


class DataLoader(Dataset):

    # The __init__ function is run once when instantiating the Dataset object.
    def __init__(self, data):
        # (n,) array, each element is string, dtype=object
        self.data = data # fasta of forward, no chr title, 1d np.array, shape is n
        print("-----------------shape before add RC -------------")
        print(self.data.shape)

        # add reverse complement
        temp = []
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N' : 'N'}
        print("reverse complement start time: ", datetime.now())
        for seq in self.data:
            complement_seq = ''
            for base in seq: ## need to check here, what is seq's shape?? why do they have seq[0] here
                complement_seq = complement[base] + complement_seq

            temp.append(complement_seq)# 0301 indented, TODO: check why before it could still train, even x and y have different n

        print("reverse complement end time: ", datetime.now())
        temp = np.array(temp, dtype=object)
        self.data = np.append(self.data, temp, axis=0)
        print("-----------------shape after subset and add RC-------------")
        print(self.data.shape)

    # The __len__ function returns the number of samples in our dataset.
    def __len__(self):
        return self.data.shape[0] ## check

    # The __getitem__ function loads and returns a sample from the dataset at the given index idx.
    def __getitem__(self, index):
        seq = self.data[index]
        row_index = 0
        temp = np.zeros((len(seq), 4))
        for base in seq: ## seq[0]??
            if base == 'A':
                temp[row_index, 0] = 1
            elif base == 'T':
                temp[row_index, 1] = 1
            elif base == 'G':
                temp[row_index, 2] = 1
            elif base == 'C':
                temp[row_index, 3] = 1
            row_index += 1

        X = torch.tensor(temp).float().permute(1,0) # change the dim to 4, 1000
        # y = torch.tensor(self.label[index]).float()

        return X