from nmr_commons.tracking.TrajectoryList import TrajectoryList
from nmr_commons.peak_picking.PeakList import PeakList
import numpy as np


class PeakTrajectoriesGenerator:

    #
    # Trajectory generation functions
    #

    # TODO
    # main function that generated list of trajectories
    def generate_peak_trajectories(self, peak_list):
        result = TrajectoryList()
        return result

    # TODO
    # n_j ~ Poi(n|N_j \lambda)
    def get_number_of_moving_peaks(self):
        pass

    # TODO
    # pkt. 2
    def get_list_of_moving_peaks(self, peak_list):
        list_of_moving_peaks = PeakList()
        return [list_of_moving_peaks, peak_list.get_list_excluding_peaks(list_of_moving_peaks)]

    # TODO
    # pkt. 3
    def get_velocity_prior_parameters(self):
        pass

    # TODO
    # pkt. 5, 6 i 7
    def get_sample_velocities(self, num_of_samples):
        rs = np.array(([]))
        ths = np.array(([]))
        ws = np.array(([]))
        return rs, ths, ws

    # TODO
    # e ~ N(e|0,b^2I)
    def get_position_noise(self, num_of_samples):
        pass

    # TODO
    # n_it ~ N(\eta, 0, \sigma^2)
    def get_velocity_noise(self, num_of_samples):
        pass

    #
    # Parameter estimation functions
    #

    # TODO
    # estimate parameters from data
    def build_generator(self, list_of_trajectory_list):
        pass