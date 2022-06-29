#import cv2
import numpy as np
import torch
from cam.activation_and_gradiant import ActivationsAndGradients

class BaseCAM:

    def __init__(self, model, target_layer, use_cuda=False):
        self.model = model.eval()
        self.taget_layer = target_layer
        self.cuda = use_cuda
        if self.cuda:
            self.model = model.cuda()

        self.activations_and_grads = ActivationsAndGradients(self.model,
                                                             target_layer)

    # make prediction
    def forward(self, input_img):
        return self.model(input_img)

    # place holder of computing the weight of cam
    # Grad-CAM: average over length, compute alpha
    def get_cam_weights(self,
                        input_tensor,
                        target_category,
                        activations,
                        grads):
        raise Exception("Not Implemented")

    # get loss of the specific target category
    def get_loss(self, output, target_category):
        return output[:, target_category]

    # main function
    def __call__(self, input_tensor, target_category=None):
        if self.cuda:
            input_tensor = input_tensor.cuda()

        # get the predicted result for input tensor
        output = self.activations_and_grads(input_tensor)

        # if not specified, target_category is the one with highest predicted prob
        if target_category is None:
            target_category = np.argmax(output.cpu().data.numpy())

        self.model.zero_grad()
        # TODO loss here is just the predicted result for given target_category?
        loss = self.get_loss(output, target_category)
        loss.backward(retain_graph=True)

        # output of last conv layer
        activations = self.activations_and_grads.activations[-1].cpu().data.numpy()[0, :]

        # gradient of output to output of last conv layer
        grads = self.activations_and_grads.gradients[-1].cpu().data.numpy()[0, :]

        # get weight (alpha)
        weights = self.get_cam_weights(input_tensor, target_category, activations, grads)

        # place holder of cam, shape is same as output of last conv (activation)
        cam = np.zeros(activations.shape[1:], dtype=np.float32)

        # weighted sum of alpha and A
        for i, w in enumerate(weights):
            # CHANGED activations[i, :, :] to activations[i, :]
            cam += w * activations[i, :]

        #print("-----------original cam ---------")
        # relu
        cam = np.maximum(cam, 0)
        # print(cam)
        # print(np.mean(cam))
        # CHANGED input_tensor.shape[2:][::-1] to input_tensor.shape[2]
        # print(input_tensor.shape)
        # cam = cv2.resize(cam, dsize=(input_tensor.shape[2],1))
        # print("-----------resize cam ---------")
        # # print(cam.shape)
        # # print(np.min(cam))
        # # print(np.max(cam))
        # # cam = cam - np.min(cam)
        # # cam = cam / np.max(cam)
        return cam


