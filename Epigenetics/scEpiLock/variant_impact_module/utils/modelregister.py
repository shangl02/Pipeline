from model.scEpiLock import scEpiLock

class ModelRegister():

    def __init__(self, model_name):
        self.model_name = model_name

    def get_model(self, n_class):

        if self.model_name == "scEpiLock":
            print("scEpiLock")
            return scEpiLock(n_class)

        else:
            print("--------No Model Retrieved!-------")
            return
