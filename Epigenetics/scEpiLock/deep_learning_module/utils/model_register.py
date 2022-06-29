from model.scEpiLock import scEpiLock
from utils.utils import Utils

class Model_Register():

    def __init__(self, model_name):
        self.model_name = model_name

    def get_model(self, n_class):

        if self.model_name == "scEpiLock":
            Utils.print_separator("scEpiLock")
            return scEpiLock(n_class)

        else:
            Utils.print_separator("No Model Retrieved!")
            return
