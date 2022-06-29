import argparse
import torch
from data.pre_processor import PreProcessor, PreProcessorTrans
from data.data_loader import DataLoader
from train.trainer import Trainer
from assess.assessor import Assessor
from utils.utils import Utils
from utils.model_register import Model_Register
from datetime import datetime

if __name__ == '__main__':
    # 1\ Parses the command line arguments and returns as a simple namespace.

    parser = argparse.ArgumentParser(description='main.py')
    parser.add_argument('-e', '--exe_mode', default='train', help='The execution mode.(train/test)')
    parser.add_argument('-c', '--config', default='./config/config_1.json', help='The config file of experiment.')
    parser.add_argument('-s', '--task', default='main', help='The task of model (transfer_learning/main)')

    args = parser.parse_args()

    # 2\ Get the configer.
    config = Utils.read_json(args.config, args.task)
    registered_model = Model_Register(config.model_name)

    # if config.cell_cluster == "all":
    #     config.cell_cluster = list(range(0, config.n_class, 1))

    config.cell_cluster = Utils.update_cell_cluster(config.cell_cluster, config.n_class)

    #print("create the output folder")
    Utils.print_separator("create the output folder")
    # if not os.path.exists(config.output_evaluation_data_path):
    #     os.makedirs(config.output_evaluation_data_path)
    Utils.create_path(config.output_evaluation_data_path)

    # 3\ Load, process and split the data set

    if args.task == "main":
        processor = PreProcessor(config.pos_forward_path, config.encode_path, config.neg_path,
                                 config.label_path, config.encode_n, config.neg_n, config.cell_cluster,
                                 config.subset)
        data, label, pos_weight = processor.concate_data()

        data_train, data_eval, data_test, label_train, \
        label_eval, label_test, test_size = processor.split_train_test(data, label)
    elif args.task == "transfer_learning":
        processor = PreProcessorTrans(config.pos_forward_path, config.neg_forward_path,
                                      config.balance, config.subset)
        data, label, pos_weight = processor.concate_data()

        data_train, data_eval, data_test, label_train, \
        label_eval, label_test, test_size = processor.split_train_test(data, label)

    # 4\ Selecting the execution mode.
    if args.exe_mode == 'train':
        Utils.print_separator("train data loader start")
        train_data_loader = DataLoader(data_train, label_train)
        Utils.print_separator("train data loader finish")

        Utils.print_separator("eval data loader start")
        eval_data_loader = DataLoader(data_eval, label_eval)
        Utils.print_separator("eval data loader finish")

        trainer = Trainer(registered_model, train_data_loader, eval_data_loader, config.model_path,
                          config.num_epochs, config.batch_size, config.learning_rate,
                          config.weight_decay, config.use_pos_weight, pos_weight,
                          config.output_evaluation_data_path, config.n_class, config.n_class_trans,
                          config.load_trans_model, config.trans_model_path)

        print("train start time: ", datetime.now())
        trainer.train() # include save model to destination
        print("train end time: ", datetime.now())

    elif args.exe_mode == 'test':
        Utils.print_separator("train data loader start")
        train_data_loader = DataLoader(data_train, label_train)
        Utils.print_separator("train data loader finish")

        Utils.print_separator("test data loader start")
        test_data_loader = DataLoader(data_test, label_test)
        Utils.print_separator("test data loader finish")

        Utils.print_separator("all data loader start")

        all_data_loader = DataLoader(data, label)

        Utils.print_separator("all data loader finish")

        assessor = Assessor(registered_model, train_data_loader, test_data_loader, all_data_loader, config.model_path,
                            config.output_evaluation_data_path, config.batch_size, test_size,
                            config.n_class)
        assessor.assess()

    else:
        print('No mode named {}.'.format(args.exe_mode))