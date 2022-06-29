import torch
import torch.nn as nn
from tqdm import tqdm
import copy
import matplotlib.pyplot as plt
import os
import numpy as np

class Trainer:
    def __init__(self, registered_model, train_data, eval_data, model_path, num_epochs, batch_size,
                 learning_rate, weight_decay, use_pos_weight, pos_weight, plot_path, n_class, n_class_trans,
                 load_trans_model, trans_model_path):

        '''
        :param model: the model
        :param train_data: complete train data, include x and y
        :param model_path: path the save the trained model
        '''

        self.registered_model = registered_model
        self.train_data = train_data
        self.eval_data = eval_data
        self._model_path = model_path
        self._num_epochs = num_epochs
        self._batch_size = batch_size
        self._learning_rate = learning_rate
        self._weight_decay = weight_decay
        self.use_pos_weight = use_pos_weight
        self.pos_weight = pos_weight
        self.plot_path = plot_path
        self.n_class = n_class
        self.n_class_trans = n_class_trans
        self.load_trans_model = load_trans_model
        self.trans_model_path = trans_model_path


    def train(self):
        """ train """
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

        # check whether to do transfer learning
        if self.load_trans_model == "True":

            model = self.registered_model.get_model(self.n_class_trans)
            # load the previously saved weight
            model.load_state_dict(torch.load(self.trans_model_path))
            # change the output layer dimension
            num_ftrs = model.Linear1.out_features
            print("-----------output for the linear1 layer from preload model-------")
            print(num_ftrs)
            model.Linear2 = nn.Linear(num_ftrs, self.n_class).to(device)

            print("----------- childern node count --------------")
            ct = 0  # total there are 9 chid for complexDQ

            for child in model.children():
                print(" child", ct, "is:")
                print(child)
                ct += 1
            print("total number of children is: ", ct)

        else:
            print("-----------no preload model loaded, load new model-------")
            model = self.registered_model.get_model(self.n_class)

        model = model.to(device)

        print("----------Current Model's state_dict-----------")
        for param_tensor in model.state_dict():
            print(param_tensor, "\t", model.state_dict()[param_tensor].size())
        print("----------Current Model's weight-----------")

        train_loss_hist = []
        eval_loss_hist = []

        train_acc_hist = []
        eval_acc_hist = []

        with open(os.path.join(self.plot_path, "train_vali_record.txt"), 'w') as file:
            file.write("Epoch \t Data \t Loss \t Acc \n")
            file.close()

        best_model_wts = copy.deepcopy(model.state_dict())
        best_loss = 1000

        if self.use_pos_weight == "False":
            print("not use pos weight in loss")
            criterion = nn.BCEWithLogitsLoss()
        elif self.use_pos_weight == "True":
            print("use pos weight in loss")
            criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([self.pos_weight]).to(device)) # TODO, compute weight
        optimizer = torch.optim.Adam(model.parameters(), lr=self._learning_rate,
                                     weight_decay=self._weight_decay, amsgrad=True)
        train_loader = torch.utils.data.DataLoader(self.train_data, batch_size=self._batch_size, shuffle=True)
        eval_loader = torch.utils.data.DataLoader(self.eval_data, batch_size=self._batch_size, shuffle=True)

        for epoch in range(self._num_epochs):
            print('Epoch {}/{}'.format(epoch, self._num_epochs - 1))
            print('-' * 10)

            # Train
            train_loss = 0.0
            train_acc = 0.0

            model.train()  # set to train mode, use dropout and batchnorm ##

            for train_step, (X, y) in enumerate(train_loader):

                X = X.to(device)
                y = y.to(device)

                # Forward pass: Compute predicted y by passing x to the model
                y_pred_logit = model(X.float()) # predicted value

                loss = criterion(y_pred_logit.float(), y.float())

                optimizer.zero_grad() # clears old gradients from the last step
                loss.backward() # for each parameter, calculate d(loss)/d(weight)
                optimizer.step() # update weights, causes the optimizer to take a step based on the gradients of the parameters

                # statistics
                y_pred_prob = 1 / (1 + np.exp(-y_pred_logit.cpu().detach().numpy())) # calculate prob from logit
                y_pred = y_pred_prob.round() ## 0./1.
                y = y.cpu().detach().numpy()

                train_loss += loss.item() * X.size(0) #loss.item() has be reducted by batch size
                if epoch == 0 and train_step == 0:
                    initial_loss = loss.item()
                    initial_acc = np.sum(y_pred == y)/float(self._batch_size*self.n_class)

                train_acc += np.sum(y_pred == y)

            train_loss = train_loss/len(train_loader.dataset)
            train_acc = train_acc/ float(len(train_loader.dataset)*self.n_class)
            print('Epoch {} {} Loss: {:.4f} Acc: {:.4f}'.format(epoch, 'train', train_loss, train_acc))

            with open(os.path.join(self.plot_path, "train_vali_record.txt"), 'a') as file:
                if epoch == 0:
                    file.write("{} \t {} \t {:.4f} \t {:.4f} \n".format(0, 'train', initial_loss, initial_acc))
                    train_loss_hist.append(initial_loss)
                    train_acc_hist.append(initial_acc)
                file.write("{} \t {} \t {:.4f} \t {:.4f} \n".format(epoch+1, 'train', train_loss, train_acc))
                file.close()

            train_loss_hist.append(train_loss)
            train_acc_hist.append(train_acc)

            # Evaluation
            eval_loss = 0
            eval_acc = 0

            model.eval() # added to not use drop out and batch norm for validation
            with torch.no_grad(): # disable gradiant calculation
                for X, y in tqdm(eval_loader):
                    optimizer.zero_grad() # make sure training and eval has minimum diff
                    X, y = X.to(device), y.to(device)
                    y_pred_logit = model(X.float())
                    eva_loss = criterion(y_pred_logit.float(), y.float())

                    y_pred_prob = 1 / (1 + np.exp(-y_pred_logit.cpu().detach().numpy()))  # calculate prob from logit
                    y_pred = y_pred_prob.round() ## 0./1.
                    y = y.cpu().detach().numpy()

                    # statistics
                    eval_loss += eva_loss.item() * X.size(0)
                    eval_acc += np.sum(y_pred == y)

                eval_loss = eval_loss / len(eval_loader.dataset)
                eval_acc = eval_acc / float(len(eval_loader.dataset)*self.n_class)

                print('{} Loss: {:.4f} Acc: {:.4f}'.format("validation", eval_loss, eval_acc))
                with open(os.path.join(self.plot_path, "train_vali_record.txt"), 'a') as file:
                    if epoch == 0:
                        # add the fake validation initial loss and acc for the plotting purpose
                        file.write("{} \t {} \t {:.4f} \t {:.4f} \n".format(0, 'validation', initial_loss, initial_acc))
                        eval_loss_hist.append(initial_loss)
                        eval_acc_hist.append(initial_acc)
                    file.write("{} \t {} \t{:.4f} \t {:.4f} \n".format(epoch+1, 'validation', eval_loss, eval_acc))
                    file.close()

            eval_loss_hist.append(eval_loss)
            eval_acc_hist.append(eval_acc)

            if eval_loss < best_loss:
                best_loss = eval_loss
                best_model_wts = copy.deepcopy(model.state_dict())

        torch.save(best_model_wts, self._model_path)
        self._plot_metrics(train_loss_hist, eval_loss_hist, "loss_history.pdf",'loss')
        self._plot_metrics(train_acc_hist, eval_acc_hist, "acc_history.pdf",'accuracy')

    def _plot_metrics(self, train_val, eval_val, plot_name, ylabel):
        plt.figure()
        plt.plot(train_val, color='green', label="train")
        plt.title('train loss')
        plt.ylabel(ylabel)
        plt.xlabel('epoch')

        plt.plot(eval_val, color='navy', label="validation")
        plt.title('metrics vs epoch')
        plt.ylabel(ylabel)
        plt.xlabel('epoch')
        plt.legend(loc='upper left')
        plt.savefig(os.path.join(self.plot_path, plot_name))


