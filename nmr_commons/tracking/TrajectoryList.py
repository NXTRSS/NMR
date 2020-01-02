import numpy as np

from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.peak_picking.PeakList import PeakList

class TrajectoryList:
    """List of Trajectory objects. If Length of trajectories will be different then class add None in the end of shorter
    trajectories then the longest"""
    STYLE_TRUE_PEAK = 'bo'
    STYLE_ARTIFACT = 'mo'
    STYLE_CONN = 'r'
    STYLE_GAP_CONN = 'r--'
    STYLE_GAP_FILLER = 'cD'

    ####
    # Constructor
    def __init__(self, list_of_trajectory=None, norm_list_of_trajectories=None, path=None, directory=None):
        self.list = list_of_trajectory if list_of_trajectory else []
        if self.list != []:
            self.normalize_trajecotories_length()
        self.norm_list = norm_list_of_trajectories if norm_list_of_trajectories else []
        if self.norm_list != []:
            self.normalize_trajecotories_length_norm_list()
        self.path = path
        self.dir = directory
        self.number_of_levels = None
        self.mean_x = None
        self.mean_y = None
        self.min_x = None
        self.min_y = None
        self.max_x = None
        self.max_y = None

    def __iter__(self):
        return self.list.__iter__()

    def __getitem__(self, item):
        return self.list[item]

    def __len__(self):
        return len(self.list)

    def normalize_trajectories_length(self):
        max_length = max([len(trajectory) for trajectory in self])
        self.number_of_levels = max_length

        for trajectory in self:
            if len(trajectory) != max_length:
                for i in range(len(trajectory), max_length):
                    trajectory.add_point(None)

    def normalize_trajectories_length_norm_list(self):
        max_length = max([len(trajectory) for trajectory in self.norm_list])

        for trajectory in self.norm_list:
            if len(trajectory) != max_length:
                for i in range(len(trajectory), max_length):
                    trajectory.add_point(None)


    def calculate_min_max(self):
        max_x_trajectories = [max(trajectory.get_all_x(with_none_points=False)) for trajectory in self]
        min_x_trajectories = [min(trajectory.get_all_x(with_none_points=False)) for trajectory in self]
        max_y_trajectories = [max(trajectory.get_all_y(with_none_points=False)) for trajectory in self]
        min_y_trajectories = [min(trajectory.get_all_y(with_none_points=False)) for trajectory in self]
        self.max_x = max(max_x_trajectories)
        self.min_x = min(min_x_trajectories)
        self.max_y = max(max_y_trajectories)
        self.min_y = min(min_y_trajectories)

    def calculate_mean_of_points(self):
        trajectory_list_x = [trajectory.get_all_x(with_none_points=False) for trajectory in self.list]
        mean_x_trajectories = [elem for trajectory in trajectory_list_x for elem in trajectory]
        self.mean_x = np.mean(mean_x_trajectories)

        trajectory_list_y = [trajectory.get_all_y(with_none_points=False) for trajectory in self.list]
        mean_y_trajectories = [elem for trajectory in trajectory_list_y for elem in trajectory]
        self.mean_y = np.mean(mean_y_trajectories)

    def normalize(self):
        if self.min_x is None or self.max_x is None or self.min_y is None or self.max_y is None:
            self.calculate_min_max()
        # if self.mean_x is None or self.mean_y is None:
        #     self.mean_of_points()

        for trajectory in self.list:
            norm_trajectory = Trajectory()
            for point in trajectory:
                if point is None:
                    norm_trajectory.add_point(None)
                else:
                    x = (point[0] - self.min_x) / (self.max_x - self.min_x)
                    y = (point[1] - self.min_y) / (self.max_y - self.min_y)
                    norm_trajectory.add_point((x, y))
            self.norm_list.append(norm_trajectory)

    def add_trajectory(self, trajectory):
        self.append(trajectory)
        if len(trajectory) != self.number_of_levels:
            self.normalize_trajectories_length()

    def add_trajectory_norm_list(self, trajectory):
        self.norm_list.append(trajectory)
        if len(trajectory) != self.number_of_levels:
            self.normalize_trajectories_length_norm_list()

    def add_trajectories(self, trajectories_list):
        for trajectory in trajectories_list:
            self.add_trajectory(trajectory)

    def add_trajectories_norm_list(self, trajectories_list):
        for trajectory in trajectories_list:
            self.add_trajectory_norm_list(trajectory)

    def get_trajectory(self, number):
        return self[number]

    def get_trajectory_norm_list(self, number):
        return self.norm_list[number]

    def split_into_peak_lists(self):

        num_of_peak_list = self.number_of_levels
        peak_lists = [0]*num_of_peak_list

        # Initialize peak list objects
        for x in range(0, num_of_peak_list):
            peak_lists[x] = PeakList(axes_names=['w1', 'w2'])

        # Iterate over spectra
        for x in range(0, num_of_peak_list):

            # Iterate over trajectories
            for trajectory in self:

                # If trajectory has peak in x'th spectrum add peak to x'th peak list
                if trajectory[x] is not None:
                    peak_lists[x].add_peak(trajectory[x])

        return peak_lists