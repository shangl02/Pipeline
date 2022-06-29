from config.config import Config
from config.config_trans import ConfigTrans
import json
import os

class Utils:

    @staticmethod
    def print_separator(input_s):
        print("-------------" + input_s + "-------------")

    @staticmethod
    def update_cell_cluster(cell_cluster, n_class):
        print(n_class)
        if cell_cluster == "all":
            cell_cluster = list(range(0, n_class, 1))
        return cell_cluster

    @staticmethod
    def create_path(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def read_json(file_path, task):
        """ Read json to dict. """
        with open(file_path, 'rt') as f:
            if task == "main":
                return json.load(f, object_hook=lambda a: Config(**a))
            elif task == "transfer_learning":
                return json.load(f, object_hook=lambda a: ConfigTrans(**a))

