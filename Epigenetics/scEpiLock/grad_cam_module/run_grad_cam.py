from cam.grad_cam import GradCAM
from model_register import Model_Register
from utils.utils import Utils
import torch
import numpy as np


class RunCam:
    def __init__(self, model_name, n_class, model_path, target_category, input):
        """
        :param model_name: name of model, eg. "DeepSea"
        :param n_class: number of cell cluster
        :param model_path: path of the pre-trained model weight
        :param target_category: target cell cluster label
        :param input: input tensor, shape (n_peak, 4, 1000)
        """
        self.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.model_name = model_name
        self.model = Model_Register(model_name).get_model(n_class)
        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.target_category = target_category
        self.input = input

    def compute_cam(self):
        Utils.print_separator("loaded model")

        for param_tensor in self.model.state_dict():
            print(param_tensor, "\t", self.model.state_dict()[param_tensor].size())

        if self.device == torch.device('cpu'):
            use_cude = False
        elif self.device == torch.device('cuda:0'):
            use_cude = True

        if self.model_name == "scEpiLock":
            target_layer = self.model.Conv4
        else:
            Utils.print_separator("No Model Retrieved!")

        cam = GradCAM(model=self.model, target_layer=target_layer, use_cuda=use_cude)

        grayscale_cam_list = []
        for i, x in enumerate(self.input):
            x = x.unsqueeze(0)
            grayscale_cam = cam(input_tensor=x, target_category=self.target_category)
            grayscale_cam_list.append(grayscale_cam)

        Utils.print_separator("grayscale_cam_list shape")
        print(len(grayscale_cam_list))
        grayscale_cam_array = np.asarray(grayscale_cam_list)
        print(grayscale_cam_array.shape)

        return(grayscale_cam_array)
