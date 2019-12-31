import numpy as np

from nmr_commons.tracking.Trajectory import Trajectory

class TrajectoryList:
    """List of Trajectory objects"""
    STYLE_TRUE_PEAK = 'bo'
    STYLE_ARTIFACT = 'mo'
    STYLE_CONN = 'r'
    STYLE_GAP_CONN = 'r--'
    STYLE_GAP_FILLER = 'cD'

    ####
    # Constructor
    def __init__(self, list_of_trajectory=None, norm_list_of_trajectories=None, path=None, directory=None):
        self.list = list_of_trajectory if list_of_trajectory else []
        self.norm_list = norm_list_of_trajectories if norm_list_of_trajectories else []
        self.path = path
        self.dir = directory
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

    def min_max(self):
        max_x_trajectories = [max(trajectory.get_all_x(with_none_points=False)) for trajectory in self.list]
        min_x_trajectories = [min(trajectory.get_all_x(with_none_points=False)) for trajectory in self.list]
        max_y_trajectories = [max(trajectory.get_all_y(with_none_points=False)) for trajectory in self.list]
        min_y_trajectories = [min(trajectory.get_all_y(with_none_points=False)) for trajectory in self.list]
        self.max_x = max(max_x_trajectories)
        self.min_x = min(min_x_trajectories)
        self.max_y = max(max_y_trajectories)
        self.min_y = min(min_y_trajectories)

    def mean_of_points(self):
        trajectory_list_x = [trajectory.get_all_x(with_none_points=False) for trajectory in self.list]
        mean_x_trajectories = [elem for trajectory in trajectory_list_x for elem in trajectory]
        self.mean_x = np.mean(mean_x_trajectories)

        trajectory_list_y = [trajectory.get_all_y(with_none_points=False) for trajectory in self.list]
        mean_y_trajectories = [elem for trajectory in trajectory_list_y for elem in trajectory]
        self.mean_y = np.mean(mean_y_trajectories)

    def normalization(self):
        if self.min_x is None or self.max_x is None or self.min_y is None or self.max_y is None:
            self.min_max()
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
