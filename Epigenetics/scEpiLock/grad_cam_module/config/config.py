class Config:
    def __init__(self,
                 pos_forward_path,
                 index,
                 model_name,
                 n_class,
                 model_path,
                 target_category,
                 subset,
                 plot,
                 result_path):
        self.pos_forward_path = pos_forward_path
        self.index = index
        self.model_name = model_name
        self.n_class = n_class
        self.model_path = model_path
        self.target_category = target_category
        self.subset = subset
        self.plot = plot
        self.result_path = result_path

