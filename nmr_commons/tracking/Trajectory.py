import csv
import copy


class Trajectory:
    """Trajectory class where trajectory is a list of points and each point is a tuple of two numbers.
    This trajectories are in the 2D space"""

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
            self.list = self.list + tuple(points)
        else:
            self.list = self.list[:start_position - 1] + points

    def return_element(self, position):
        return self.list[position]  # I'm not sure if this is necessary - getitem doesn't do the same?

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

    def get_next_point(self, point):
        index = None
        for i in range(point + 1, len(self.list)):
            if self.list[i] is not None:
                index = i
                break
        return index

    def get_distances_old(self):
        distance = list()
        for i in range(len(self.list) - 1):
            if self.list[i] is not None:
                k = i + 1
                while k != len(self.list) - 1 and self.list[k] is None:
                    k += 1
                if k == len(self.list) - 1 and self.list[k] is None:
                    return distance
                else:
                    for j in range(i, k):
                        distance.append([(self.list[i][0] - self.list[k][0]) / (k - i),
                                         (self.list[i][1] - self.list[k][1]) / (k - i)])
        return distance

    # def get_distances(self):
    #     distances = [None]*(len(self)-1)
    #     for i, point in enumerate(self[:-1]):
    #         if self[i] is not None:
    #             k = i + 1
    #             while k != (len(self)-1) and self[k] is None:
    #                 k += 1
    #             if k == len(self) - 1 and self[k] is None:
    #                 return distances
    #             else:
    #                 for j in range(i, k):
    #                     distances[j] =
