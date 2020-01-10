import csv
import copy
import numpy as np


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

    def __eq__(self, other):
        return self[:] == other[:]

    def get_x(self, position):
        if self[position] is not None:
            return self[position][0]
        return None

    def get_y(self, position):
        if self[position] is not None:
            return self[position][1]
        return None

    def get_vector(self, position):
        for i, point in enumerate(self[position:-1]):
            if point is not None and self[position + i + 1] is not None:
                return position + i
        return None

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
        if self[:position] == [None] * len(self[:position]) and self[position + 1:] == [None] * len(
                self[position + 1:]):
            raise ValueError('You can not remove all elements from trajectory')
        else:
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

    def get_next_point(self, position):
        index = None
        for i in range(position + 1, len(self)):
            if self[i] is not None:
                index = i
                break
        return index

    def get_distances(self):
        distances = list()
        for i, point in enumerate(self[:-1]):
            if self[i] is not None:
                k = i + 1
                if self[k] is not None:
                    distances.append(((point[0] - self.get_x(k)),
                                    (point[1] - self.get_y(k))))
                else:
                    while k != (len(self) - 1) and self[k] is None:
                        k += 1
                    if k == len(self) - 1 and self[k] is None:
                        return distances
                    else:
                        for j in range(i, k):
                            distances.append(((point[0] - self.get_x(k)) / (k - i),
                                            (point[1] - self.get_y(k)) / (k - i)))
        return distances

    def get_distances_sep(self):
        distances = list()
        for i, point in enumerate(self[:-1]):
            if self.get_vector(i) == i:
                distances.append(tuple((point[0] - self.get_x(i + 1), point[1] - self.get_y(i + 1))))
        return distances

    @staticmethod
    def unit_vector(vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def get_angles(self):
        angles = []
        for i in range(len(self) - 2):
            idx_first_vec = self.get_vector(i)
            if idx_first_vec is not None:
                a = self[idx_first_vec]
                b = self[idx_first_vec + 1]

                idx_second_vec = self.get_vector(idx_first_vec + 1)
                if idx_second_vec is not None:
                    c = self[idx_second_vec]
                    d = self[idx_second_vec + 1]
                    v1 = self.unit_vector(tuple([x1 - x2 for x1, x2 in zip(b, a)]))
                    v2 = self.unit_vector(tuple([x1 - x2 for x1, x2 in zip(d, c)]))
                    angles.append(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))
        return angles

    def get_all_x(self, with_none_points=True):
        if with_none_points:
            return [point[0] if point is not None else None for point in self]
        else:
            return [point[0] for point in self if point is not None]

    def get_all_y(self, with_none_points=True):
        if with_none_points:
            return [point[1] if point is not None else None for point in self]
        else:
            return [point[1] for point in self if point is not None]

    def get_mean(self):
        xn = self.get_all_x(with_none_points=False)
        yn = self.get_all_y(with_none_points=False)
        return [np.mean(xn), np.mean(yn)]

    def can_merge(self, trj):
        min_idx = min(len(self), len(trj))
        for i in range(min_idx):
            if (self[i] is not None) and (trj[i] is not None):
                return False
        return True

    def merge(self, trj):
        for i in range(min(len(self), len(trj))):
            if trj[i] is not None:
                self.add_point(trj[i], i)

    def get_first_last(self):
        not_none_idxs = [i for i in range(len(self)) if self[i] is not None]
        first_idx = not_none_idxs[0]
        last_idx = not_none_idxs[-1]
        first_point = self[first_idx]
        last_point = self[last_idx]

        return first_idx, last_idx, first_point, last_point

    def get_not_none_count(self):
        return sum(point is not None for point in self)

    def get_path_ratio(self, x_div=35.0, y_div=6.0):
        ratio = 1
        sum_dist = 0
        first_last_dist = 0
        if self.get_not_none_count() > 2:
            (first, last, first_point, last_point) = self.get_first_last()
            first_last_dist = np.sqrt(((first_point[0] - last_point[0]) / x_div) ** 2 +
                                      ((first_point[1] - last_point[1]) / y_div) ** 2)
            if first_last_dist > 0:
                for i in range(first + 1, last + 1):
                    if self[i] is not None:
                        last_point = self[i]
                        sum_dist = sum_dist + np.sqrt(((first_point[0] - last_point[0]) / x_div) ** 2
                                                      + ((first_point[1] - last_point[1]) / y_div) ** 2)
                        first_point = last_point
                ratio = float(sum_dist) / float(first_last_dist)
        return ratio, sum_dist, first_last_dist

    def get_stats(self, x_div, y_div):

        v_mod_mean = 0
        v_mod_std = 0
        angle_mean = 0
        angle_std = 0
        v_s = []
        angles = []
        if (self.get_not_none_count() > 2):
            (first, last, first_point, last_point) = self.get_first_last()
            div = 1.0
            for k in range(first + 1, last + 1):
                if self.list[k] is not None:
                    last_point = self.list[k];
                    dist = np.sqrt(((first_point[0] - last_point[0]) / x_div) ** 2 +
                                ((first_point[1] - last_point[1]) / y_div) ** 2) / div
                    angle = np.abs(np.arctan2(((first_point[1] - last_point[1]) / y_div),
                                              ((first_point[0] - last_point[0]) / x_div)))
                    v_s.append(dist)
                    angles.append(angle)
                    first_point = last_point
                    div = 1.0
                else:
                    div += 1.0

            v_mod_mean = np.mean(v_s)
            v_mod_std = np.std(v_s)
            angle_mean = np.mean(angles)
            angle_std = np.std(angles)
        return v_mod_mean, v_mod_std, angle_mean, angle_std, v_s, angles

    def is_correct(self, x_div=35.0, y_div=6.0):
        v_s = []
        angles = []
        if self.get_not_none_count() > 2:
            (first, last, first_point, last_point) = self.get_first_last()
            div = 1.0
            for k in range(first + 1, last + 1):
                if self[k] is not None:
                    last_point = self[k]
                    d = np.sqrt(((first_point[0] - last_point[0]) / x_div) ** 2 +
                                ((first_point[1] - last_point[1]) / y_div) ** 2) / div
                    en = np.abs(np.arctan2(((first_point[1] - last_point[1]) / y_div),
                                           ((first_point[0] - last_point[0]) / x_div)))
                    if len(v_s) >= 2:
                        if d > 0.1:
                            if np.abs(v_s[-2] - d) > 1:  # or (np.abs(engles[-2]-en)>0.8)):
                                return False

                    v_s.append(d)
                    angles.append(en)
                    first_point = last_point
                    div = 1.0
                else:
                    div += 1.0
        return True

    def init_nones(self, K):
        for k in range(K):
            self.add_point(None)

    @staticmethod
    def read_from_csv(row):
        line = Trajectory()
        for elem in row:
            if elem != '':
                corrected_elem = eval(elem)
                line.add_point(corrected_elem)
            else:
                line.add_point(None)

        return line
