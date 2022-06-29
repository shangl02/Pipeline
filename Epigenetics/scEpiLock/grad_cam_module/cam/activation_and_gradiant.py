class ActivationsAndGradients:
    """ Class for extracting activations and
    registering gradients from targeted intermediate layers """

    def __init__(self, model, target_layer):
        self.model = model
        self.gradients = []
        self.activations = []

        # register save_activation function as a callable of target layer
        target_layer.register_forward_hook(self.save_activation)
        target_layer.register_backward_hook(self.save_gradient)

    # save the output of the last conv
    def save_activation(self, module, input, output):
        self.activations.append(output)

    # save the gradient of last-layer_output/last-conv_output
    # will be used in loss.backward(retain_graph=True)!!
    def save_gradient(self, module, grad_input, grad_output):
        # Gradients are computed in reverse order

        self.gradients = [grad_output[0]] + self.gradients

    def __call__(self, x):
        self.gradients = []
        self.activations = []
        return self.model(x)