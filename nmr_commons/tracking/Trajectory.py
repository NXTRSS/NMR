import csv
import copy


class Trajectory:
    """Trajectory class where trajectory is a list of point and each point is a tuple of two points in space"""
    ####
    # Constructor
    def __init__(self, items=None, norm_items=None):
        # self.list = items if items else []    #this was before, I worried about mutable parameters
        self.list = list(items) if items else []
        # self.norm_list = norm_items if norm_items else []     #this was before, I worried about mutable parameters
        self.norm_list = copy.deepcopy(norm_items) if norm_items else []
        self.flag = None  # used in trajectories visualization,
        # True = trajectory created from real peak, False = trajectory created from artifact

    def __iter__(self):
        return iter(self.list)

    def __getitem__(self, item):
        return self.list[item]

    def __len__(self):
        return len(self.list)

    def add_point(self, point, position=None):
        if position is not None:
            self.list[position] = tuple(point)
        else:
            self.list.append(point)

    def add_points(self, points, start_position=None):
        if start_position is None:
            self.list = self.list+tuple(points)
        else:
            self.list = self.list[:start_position-1]+points

    def return_element(self, position):
        return self.list[position]          #I'm not sure if this is necessary - getitem doesn't do the same?

    def remove_point(self, position):
        self.list[position] = None

    def remove_points(self, position_list):
        for position in position_list:
            self.remove_point(position)

    def print_trajectory(self, position=None):
        if position is None:
            print(self.list)
        else:
            print(self.list[position])

    def to_csv(self, path):
        with open(path) as output:
            writer = csv.writer(output, lineterminator="\n", delimiter=';')
            writer.writerow(self.list)
