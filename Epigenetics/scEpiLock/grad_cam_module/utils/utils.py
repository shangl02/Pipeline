from config.config import Config
import json

class Utils:
    def print_separator(input_s):
        print("-------------" + input_s + "-------------")

    @staticmethod
    def read_json(file_path):
        """ Read json to dict. """
        with open(file_path, 'rt') as f:
            return json.load(f, object_hook=lambda a: Config(**a))




