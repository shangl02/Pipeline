from utils.utils import Utils
import sys
import numpy as np
import os
from dataloader import DataLoader
from run_grad_cam import RunCam
from utils.plot import Plot

if __name__== '__main__':
    # 1\ Parses the command line arguments and returns as a simple namespace.
    config = Utils.read_json(sys.argv[1])
    print("create the output folder")
    if not os.path.exists(config.result_path):
        os.makedirs(config.result_path)

    # 2\ Load and process seq
    data = DataLoader(config.pos_forward_path, config.index, config.subset)
    seq, complement_seq = data.get_seq()
    data_tensor = data.seq_one_hot_tensor(seq)
    complement_data_tensor = data.seq_one_hot_tensor(complement_seq)

    # 3\ run Grad-CAM
    for i in config.target_category:
        grad_cam = RunCam(config.model_name, config.n_class, config.model_path,
                          i, data_tensor)
        cam = grad_cam.compute_cam()
        target_name = "cluster_"+str(i)
        np.savetxt(os.path.join(config.result_path, target_name+"_forward_cam.csv"), cam, delimiter=",")

        grad_cam_complement = RunCam(config.model_name, config.n_class, config.model_path,
                          i, complement_data_tensor)
        complement_cam = grad_cam_complement.compute_cam()
        np.savetxt(os.path.join(config.result_path, target_name+"_complement_cam.csv"), complement_cam, delimiter=",")

        # 4\ plot the computed cam
        if config.plot == "True":
            Plot(cam, config.result_path, target_name+"_forward").plot_cam(1)
            Plot(complement_cam, config.result_path, target_name+"_reverse_complement").plot_cam(2)
