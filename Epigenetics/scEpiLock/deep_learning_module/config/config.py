class Config():

    def __init__(self,
                 pos_forward_path,
                 encode_path,
                 neg_path,
                 label_path,
                 encode_n,
                 neg_n,
                 # X_train_data_path, # TODO: question, why these four does not have property functions
                 # X_test_data_path,
                 # y_train_data_path,
                 # y_test_data_path,
                 model_name,
                 load_trans_model,
                 num_epochs,
                 batch_size,
                 learning_rate,
                 weight_decay,
                 use_pos_weight,
                 model_path,
                 output_evaluation_data_path,
                 trans_model_path,
                 subset,
                 cell_cluster,
                 n_class,
                 n_class_trans,
                 balance):

        self.pos_forward_path = pos_forward_path
        self.encode_path = encode_path
        self.label_path = label_path
        self.neg_path = neg_path
        self.encode_n = encode_n
        self.neg_n = neg_n
        # self.X_train_data_path = X_train_data_path
        # self.X_test_data_path = X_test_data_path
        # self.y_train_data_path = y_train_data_path
        # self.y_test_data_path = y_test_data_path
        self._model_name = model_name
        self.load_trans_model = load_trans_model
        self._num_epochs = num_epochs
        self._batch_size = batch_size
        self._learning_rate = learning_rate
        self._weight_decay = weight_decay
        self.use_pos_weight = use_pos_weight
        self._model_path = model_path
        self._output_evaluation_data_path = output_evaluation_data_path
        self.trans_model_path = trans_model_path
        self._subset = subset
        self.cell_cluster = cell_cluster
        self.n_class = n_class
        self.n_class_trans = n_class_trans
        self.balance = balance

    @property
    def input_train_data_path(self):
        return self._input_train_data_path

    # @input_train_data_path.setter
    # def input_train_data_path(self, input_train_data_path):
    #     self._input_train_data_path = input_train_data_path

    @property
    def input_test_data_path(self):
        return self._input_test_data_path

    # @input_test_data_path.setter
    # def input_test_data_path(self, input_test_data_path):
    #     self._input_test_data_path = input_test_data_path

    @property
    def output_evaluation_data_path(self):
        return self._output_evaluation_data_path

    # @output_evaluation_data_path.setter
    # def output_evaluation_data_path(self, output_evaluation_data_path):
    #     self._output_evaluation_data_path = output_evaluation_data_path

    @property
    def model_name(self):
        return self._model_name

    @property
    def num_epochs(self):
        return self._num_epochs

    @property
    def batch_size(self):
        return self._batch_size

    @property
    def learning_rate(self):
        return self._learning_rate

    @property
    def weight_decay(self):
        return self._weight_decay

    # @model_name.setter
    # def model_name(self, model_name):
    #     self._model_name = model_name

    @property
    def model_path(self):
        return self._model_path

    # @model_path.setter
    # def model_path(self, model_path):
    #     self._model_path = model_path

    @property
    def subset(self):
        return self._subset



    # @subset.setter
    # def subset(self, subset):
    #     self._subset = subset