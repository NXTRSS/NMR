import numpy as np

from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.tracking.TrajectoryList import TrajectoryList

class TrajectoryGenerator:
    """Class for training and generating artificial list of Trajectories. This is a main class in data augmentation task
    for NMR protein purpose."""

    ####
    # Constructor
    def __init__(self, number_of_levels=None, path_list=None, curving_lim=0.12):
        self.levels_number = number_of_levels
        self.lam = 0
        self.epsilon = 0
        self.mi = [0, 0]
        self.ni = 0
        self.r = 0
        self.k0 = 0
        self.ksi0 = 0
        self.ksi = 0
        self.Ibeta = [[0, 0], [0, 0]]
        self.sigma = 0
        self.r = []
        self.path_list = path_list
        self.angle_dist = None
        self.angle_dist_params = [0, 0]
        self.gamma = 0
        self.curving_lim = curving_lim
        self.chem_shift_gen = None


