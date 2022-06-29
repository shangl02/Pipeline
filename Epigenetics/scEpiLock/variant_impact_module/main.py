import sys
import os
from utils.utils import Utils
from data.data_loader import DataLoader
from data.generate_fasta import FaGenerator
from utils.modelregister import ModelRegister
from evaluate.evaluator import Evaluator

if __name__== '__main__':
    # 1\ Parses the command line arguments and returns as a simple namespace.
    config = Utils.read_json(sys.argv[1])
    print("create the output folder")
    if not os.path.exists(config.result_path):
        os.makedirs(config.result_path)

    registered_model = ModelRegister(config.model_name)

    # 2\ Read bed and snp file, and generate the fasta
    fa_generator = FaGenerator(config.peak_snp_bed, config.ref_fasta, config.result_path)
    combined = fa_generator.read_bed()
    snp_pos = fa_generator.compute_pos_index(combined)
    wt_fa, mutant_fa = fa_generator.get_fasta(snp_pos, combined)
    # 3\ construct data loader
    wt_loader = DataLoader(wt_fa)
    mutant_loader = DataLoader(mutant_fa)

    # 4\ compute y_hat (and delta_y)
    evaluator = Evaluator(registered_model, config.n_class, config.model_weight_path, wt_loader, mutant_loader, config.result_path)
    evaluator.evaluate()