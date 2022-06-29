import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
import os

class Plot:
    def __init__(self, cam, result_path, plot_name):
        self.cam = cam[0]
        self.result_path = result_path
        self.plot_name = plot_name

    def plot_cam(self, i):
        x = np.linspace(0, 999, len(self.cam))
        y = self.cam
        print(len(x))
        print(len(y))
        spl = InterpolatedUnivariateSpline(x, y)
        plt.figure(i)
        plt.plot(x, y, 'ro', ms=5)
        xs = np.linspace(0, 999, 1000)
        plt.plot(xs, spl(xs), 'g', lw=3, alpha=0.7)

        plot_name = self.plot_name + "_" + "cam.pdf"
        plt.savefig(os.path.join(self.result_path, plot_name))
