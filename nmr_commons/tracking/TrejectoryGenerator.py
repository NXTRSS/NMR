import numpy as np
from os import listdir, walk
from os.path import isfile, join
import re
import operator
import matplotlib.pyplot as plt

from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.tracking.TrajectoryList import TrajectoryList
from nmr_environment import Settings
from nmr_commons.utils import flatten

class TrajectoryGenerator:
    """Class for training and generating artificial list of Trajectories. This is a main class in data augmentation task
    for NMR protein purpose."""

    ####
    # Constructor
    def __init__(self, number_of_levels=None, path_list=None, curving_lim=0.12, k0=5, ksi0=0.0002):
        self.levels_number = number_of_levels
        self.lam = 0 #ratio of moving trajectories vs all trajectories in trajectories lists
        self.epsilon = 0
        self.mi = [0, 0]
        self.ni = 0
        self.r = 0
        self.k0 = k0
        self.ksi0 = ksi0
        self.ksi = 0
        self.Ibeta = [[0, 0], [0, 0]] # covariance matrix of noise. Noise calculated only on staying trajectories
        self.sigma = 0
        self.r = []
        self.path_list = path_list
        self.list_of_trajectory_list = []
        self.angle_dist = None
        self.angle_dist_params = [0, 0]
        self.gamma = 0
        self.curving_lim = curving_lim
        self.chem_shift_gen = None

    def train(self, verbose=False):
        if self.list_of_trajectory_list==[]:
            self.list_of_trajectory_list = self.get_list_of_trajectory_list(verbose)
        self.lam = self.calculate_lambda()
        lam = self.lam

        self.Ibeta = self.calculate_noise_cov()
        self.calculate_velocity()
        self.calculate_angle_noise()
        if verbose:
            print('Lambda: {:.4f}, k0: {:.4f}, ksi0: {:.4f}'.format(lam, self.k0, self.ksi0))
        # self.k0 = 5
        # self.ksi0 = 0.0002

    def get_list_of_trajectory_list(self, verbose=False):
        if self.path_list is None:
            self.path_list = get_path_list()
            if verbose:
                print(self.path_list)
        list_of_trajectory_list = []
        for path in self.path_list:
            trajectory_list = TrajectoryList.read_csv(path)
            list_of_trajectory_list.append(trajectory_list)
        return list_of_trajectory_list

    def calculate_lambda(self):
        n = 0
        N = 0
        for trajectory_list in self.list_of_trajectory_list:
            n += len(trajectory_list.get_moving_trajectories())
            N += len(trajectory_list)
        return float(n) / N

    def calculate_noise_cov(self):
        distances = list()
        for trajectory_list in self.list_of_trajectory_list:
            staying = trajectory_list.get_staying_trajectories()
            distances += trajectory_list.get_distances_min_max_norm(interval=staying)
        distances = flatten(distances)
        var = np.var(distances, axis=0)

        return [[var[0], 0], [0, var[1]]]

def get_path_list(general_path=None):
    if general_path is None:
        path_list = Settings.TRAJECTORIES_PATH
    else:
        path_list = {}
        for root, dirs, files in walk(general_path):
            for file in files:
                if isfile(join(root, file)) and file.endswith('.csv') and 'Tracking' in file:
                    last_directory = root.split('/')[-1]
                    last_directory_number = int(re.findall(r'\d+', last_directory)[0])
                    path_list[last_directory_number] = join(root, file)
        path_list = [path for number, path in sorted(path_list.items(), key=operator.itemgetter(0))]
    return path_list
