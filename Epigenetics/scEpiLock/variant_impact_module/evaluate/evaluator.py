import torch
import os
import numpy as np
import pandas as pd
from tqdm import tqdm


class Evaluator:
    def __init__(self, model, n_class, model_path, wt_loader, mutant_loader, result_path):
        self.model = model.get_model(n_class)
        self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.model.load_state_dict(torch.load(model_path,map_location=self.device))
        self.n_class = n_class
        self.wt_loader = wt_loader
        self.mutant_loader = mutant_loader
        self.result_path = result_path

    def evaluate(self):
        wt_loader = torch.utils.data.DataLoader(self.wt_loader, batch_size=1, shuffle=False)
        mutant_loader = torch.utils.data.DataLoader(self.mutant_loader, batch_size=1, shuffle=False)

        wt_yhat = self._predict(wt_loader)
        mutant_yhat = self._predict(mutant_loader)
        diff_yhat = np.subtract(wt_yhat, mutant_yhat)

        pd.DataFrame(wt_yhat).to_csv(os.path.join(self.result_path, "wt_yhat.csv"))
        pd.DataFrame(mutant_yhat).to_csv(os.path.join(self.result_path, "mutant_yhat.csv"))
        pd.DataFrame(diff_yhat).to_csv(os.path.join(self.result_path, "diff_yhat.csv"))

    def _predict(self, data_loader):
        self.model.eval()

        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.model.to(device)
        print(self.model.type)
        with torch.no_grad():

            y_hat_arr = np.empty((0, self.n_class))

            for X in tqdm(data_loader):
                X = X.to(device)
                y_hat = self.model (X.float())

                y_hat_arr = np.concatenate((y_hat_arr, y_hat.cpu().numpy()))

            print("------------------ data shape --------------------")
            print(y_hat_arr.shape)
            print(type(y_hat_arr)) #<class 'numpy.ndarray'>
            # print(y_hat_arr)
            y_hat_arr = 1 / (1 + np.exp(-y_hat_arr))  # 0320 evening, add sigmoid
            print(y_hat_arr)
            print("------------------ train/test prediction done --------------------")

            return y_hat_arr