import numpy as np
import csv
import os
import glob
import numbers
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from collections import OrderedDict


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
    MOVING_LIM_PARAMETER = 0.0006

    ####
    # Constructor
    def __init__(self, list_of_trajectories=None, norm_list_of_trajectories=None, path=None, directory=None):
        self.number_of_levels = None
        self.list = list_of_trajectories if list_of_trajectories else []
        if self.list:
            self.number_of_levels, self.list = normalize_trajectories_length(list_of_trajectories)
        self.norm_list = norm_list_of_trajectories if norm_list_of_trajectories else []
        if self.norm_list:
            _, self.norm_list = normalize_trajectories_length(self.norm_list)
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

    # def __getitem__(self, item):
    #     return self.list[item]

    def __getitem__(self, index):
        cls = type(self)
        if isinstance(index, slice):
            return cls(list_of_trajectories=self.list[index])
        elif isinstance(index, list) or isinstance(index, range):
            if isinstance(index, range):
                index = list(index)
            return cls(list_of_trajectories=[self[idx] for idx in index])
        elif isinstance(index, numbers.Integral):
            return self.list[index]
        else:
            msg = '{cls.__name__} indices must be integers'
            raise TypeError(msg.format(cls=cls))

    def __len__(self):
        return len(self.list)

    def calculate_min_max(self):
        max_x_trajectories = [max(trajectory.get_all_x(with_none_points=False)) for trajectory in self]
        min_x_trajectories = [min(trajectory.get_all_x(with_none_points=False)) for trajectory in self]
        max_y_trajectories = [max(trajectory.get_all_y(with_none_points=False)) for trajectory in self]
        min_y_trajectories = [min(trajectory.get_all_y(with_none_points=False)) for trajectory in self]
        self.max_x = max(max_x_trajectories)
        self.min_x = min(min_x_trajectories)
        self.max_y = max(max_y_trajectories)
        self.min_y = min(min_y_trajectories)

    def get_trajectories_bounding_box(self):
        bounding_box_list = [None]*len(self)
        for i, trajectory in enumerate(self):
            x_min, y_min = trajectory.get_min()
            x_max, y_max = trajectory.get_max()
            bounding_box_list[i] = (x_min, x_max, y_min, y_max)
        return bounding_box_list

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
            self.number_of_levels, self.list = normalize_trajectories_length(self.list)

    def add_trajectory_norm_list(self, trajectory):
        self.norm_list.append(trajectory)
        if len(trajectory) != self.number_of_levels:
            _, self.norm_list = normalize_trajectories_length(self.norm_list)

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

    def get_distances_old(self, interval=None, normalization=None, mode=1):
        distances = []
        if interval is None:
            interval = range(len(self.list))

        if normalization is None:
            if mode == 1:  # euclidean distances
                for k in interval:
                    trajectory = self.get_trajectory(k)
                    trajectory_dist = trajectory.get_distances()
                    one = []
                    for i in range(len(trajectory_dist))[:-1]:
                        one.append(np.sqrt((trajectory_dist[i + 1][0] - trajectory_dist[i][0]) ** 2 +
                                           (trajectory_dist[i + 1][1] - trajectory_dist[i][1]) ** 2))
                    distances.append(one)
            else:
                self.calculate_min_max()
                for i in interval:
                    trajectory = self.get_trajectory(i)
                    value = [tuple([a[0] / (self.max_x - self.min_x), a[1] / (self.max_y - self.min_y)])
                             for a in trajectory.get_distances_sep()]
                    distances.append(value)

        else:
            if mode == 1:  # euclidean distances
                for k in interval:
                    trajectory = self.norm_list[k]
                    trajectory_dist = trajectory.get_distances()
                    one = []
                    for i in range(len(trajectory_dist))[:-1]:
                        one.append(np.sqrt((trajectory_dist[i + 1][0] - trajectory_dist[i][0]) ** 2 +
                                           (trajectory_dist[i + 1][1] - trajectory_dist[i][1]) ** 2))
                    distances.append(one)
            else:
                for i in interval:
                    trajectory = self.get_trajectory(i)
                    value = [tuple([a[0] / (self.max_x - self.min_x), a[1] / (self.max_y - self.min_y)])
                             for a in trajectory.get_distances_sep()]
                    distances.append(value)

        return distances

    def get_distances_euclidean(self, interval=None, normalization=None):
        distances = []
        if interval is None:
            interval = range(len(self)) #take all trajectories in trajectories list

        if normalization is None:
            for trajectory in self[interval]:
                trajectory_dist = trajectory.get_distances()
                one_trajectory_list = []
                for i in range(len(trajectory_dist))[:-1]:
                    one_trajectory_list.append(np.sqrt((trajectory_dist[i + 1][0] - trajectory_dist[i][0]) ** 2 +
                                                       (trajectory_dist[i + 1][1] - trajectory_dist[i][1]) ** 2))
                distances.append(one_trajectory_list)

        else:
            if not self.norm_list:
                self.normalize()
            for trajectory in [self.norm_list[idx] for idx in interval]:
                trajectory_dist = trajectory.get_distances()
                one_trajectory_list = []
                for i in range(len(trajectory_dist))[:-1]:
                    one_trajectory_list.append(np.sqrt((trajectory_dist[i + 1][0] - trajectory_dist[i][0]) ** 2 +
                                                       (trajectory_dist[i + 1][1] - trajectory_dist[i][1]) ** 2))
                distances.append(one_trajectory_list)

        return distances

    def get_distances_min_max_norm(self, interval=None, normalization=None):
        distances = []
        if interval is None:
            interval = range(len(self)) #take all trajectories in trajectories list

        if self.min_x is None or self.max_x is None or self.min_y is None or self.max_y is None:
            self.calculate_min_max()
        for trajectory in self[interval]:
            value = [tuple([a[0] / (self.max_x - self.min_x), a[1] / (self.max_y - self.min_y)])
                     for a in trajectory.get_distances_sep()]
            distances.append(value)

        return distances

    def get_mean_angles_of_trajectory(self, interval=None):
        angles = []
        if interval is None:
            interval = range(len(self))
        for trajectory in self[interval]:
            angles_in_trajectory = trajectory.get_angles()
            if angles_in_trajectory != []:
                mean = np.mean(angles_in_trajectory)
            else:
                mean = None
            angles += [mean]
        return angles

    def get_distances_first_point_from_mean(self, interval=None):
        distances = []
        if interval is None:
            interval = range(len(self))
        for trajectory in self[interval]:
            [xm, ym] = trajectory.get_mean()
            x0 = trajectory.get_all_x(with_none_points=False)[0]
            y0 = trajectory.get_all_y(with_none_points=False)[0]
            # xl = [x for x in trajectory.get_x() if x is not None][-1]
            # yl = [y for y in trajectory.get_y() if y is not None][-1]
            # a = np.sqrt((x0-xm)**2+(y0-ym)**2)
            # b = np.sqrt((xl-xm)**2+(yl-ym)**2)
            distances.append([abs(x0-xm), abs(y0-ym)])

        return distances

    def get_distances_from_start_to_end(self, interval=None):
        distances = []
        if interval is None:
            interval = range(len(self))
        for trajectory in self[interval]:
            x_all = trajectory.get_all_x(with_none_points=False)
            x0 = x_all[0]
            xl = x_all[-1]

            y_all = trajectory.get_all_y(with_none_points=False)
            y0 = y_all[0]
            yl = y_all[-1]
            # a = np.sqrt((x0-xm)**2+(y0-ym)**2)
            # b = np.sqrt((xl-xm)**2+(yl-ym)**2)
            distances.append([abs(x0 - xl), abs(y0 - yl)])

        return distances

    def create_gap(self, gap_probability):
        for trajectory in self:
            for point in range(len(trajectory)):
                if np.random.uniform() <= gap_probability:
                    trajectory.remove_point(point)

    def trajectories_moving_staying_idx(self, directory=None, directory_path=None):
        lim = self.MOVING_LIM_PARAMETER * self.number_of_levels

        moving = []
        staying = []
        if self.max_x or self.max_y is None:
            self.calculate_min_max()
        if directory is None:
            for k, distances in enumerate(self.get_distances_from_start_to_end()):
                if np.sqrt((distances[0]/(self.max_x-self.min_x))**2 + (distances[1]/(self.max_y-self.min_y))**2) > lim:
                    moving.append(k)
                else:
                    staying.append(k)

        else:
            distances = self.get_distances_first_point_from_mean()
            distances_e = [(dist[0] / self.max_x[0]) ** 2 + (dist[1] / self.max_y[0]) ** 2 for dist in distances]
            idx = [i[0] for i in sorted(enumerate(distances_e), key=lambda x:x[1])]

            directory = self.dir if directory_path is None else directory_path
            os.chdir(directory)
            files = [file_name for file_name in glob.glob("*.txt")]
            file_name = open(files[0], 'r')
            n = file_name.read()

            moving = idx[0:n]
            staying = idx[n+1:]
        return [moving, staying]
    
    def get_moving_trajectories(self, directory=None):
        return self.trajectories_moving_staying_idx(directory)[0]

    def get_staying_trajectories(self, directory=None):
        return self.trajectories_moving_staying_idx(directory)[1]
    
    def get_crossing_trajectories_idx(self):
        crossing = []
        moving = self.get_moving_trajectories()
        bounding_boxes = TrajectoryList.get_trajectories_bounding_box(self[moving])

        for i, first_idx in enumerate(moving):
            trajectory1 = self[first_idx]
            bb1 = bounding_boxes[i]
            for j, second_idx in enumerate(moving[i:]):
                trajectory2 = self[second_idx]
                bb2 = bounding_boxes[j]

                if bb1[0] > bb2[1] or bb1[1] < bb2[0] or bb1[2] > bb2[3] or bb1[3] < bb2[2]:
                    continue  # no crossing for sure

                for k in range(len(trajectory1) - 1):
                    for l in range(len(trajectory2) - 1):
                        a_start_in = trajectory1[k]
                        a_next = trajectory1.get_next_point(k)

                        b_start_in = trajectory2[l]
                        b_next = trajectory2.get_next_point(l)

                        if a_next is not None and b_next is not None:
                            a_end_in = trajectory1[a_next]
                            b_end_in = trajectory2[b_next]

                            if cross_check_of_two_points(a_start_in, a_end_in, b_start_in, b_end_in):
                                crossing.append((moving[first_idx], moving[second_idx]))

        return crossing

    def split_into_peak_lists(self):

        num_of_peak_list = self.number_of_levels
        peak_lists = [0] * num_of_peak_list

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

    def visualize_3d(self, show=None, level_interval=None, trajectories_interval=None):
        if level_interval is None:
            level_interval = list(range(self.number_of_levels))
        if trajectories_interval is None:
            trajectories_interval = range(len(self))
        fig1 = plt.figure()
        ax = fig1.add_subplot(1, 1, 1, projection='3d')
        for trajectory in self[trajectories_interval]:
            not_none_point = None
            for idx in level_interval:
                point = trajectory[idx]
                if point is not None:

                    ax.scatter(point[0], point[1], idx)
                    if not_none_point is None:
                        not_none_point = (point, idx)
                    else:
                        previous_point, previous_idx = not_none_point
                        ax.plot([point[0], previous_point[0]],
                                [point[1], previous_point[1]], [idx, previous_idx], alpha=0.3, color='red')
                        not_none_point = (point, idx)

        if show is not None:
            plt.show()
        return fig1

    def visualize_2d(self, show=None, title=None, level_interval=None, trajectories_interval=None, save_path=None, formats=None):
        if level_interval is None:
            level_interval = range(self.number_of_levels)

        if trajectories_interval is None:
            trajectories_interval = range(len(self))

        fig = plt.gcf()
        plot_queue = OrderedDict([
            (TrajectoryList.STYLE_GAP_FILLER, []),
            (TrajectoryList.STYLE_ARTIFACT, []),
            (TrajectoryList.STYLE_TRUE_PEAK, []),
            (TrajectoryList.STYLE_CONN, []),
            (TrajectoryList.STYLE_GAP_CONN, [])
        ])

        for i, trajectory in enumerate(self[trajectories_interval]):
            plot_queue, start_point = add_trajectory_to_plot_queue(trajectory, level_interval, plot_queue)
            if start_point is not None:
                plt.text(start_point[0], start_point[1], str(i), color='green', size=10)

        for style, points in plot_queue.items():
            if len(points) > 0:
                if style in [TrajectoryList.STYLE_CONN, TrajectoryList.STYLE_GAP_CONN]:
                    for elem in points:
                        xs, ys = zip(*elem)
                        plt.plot(xs, ys, style, alpha=0.8)
                else:
                    xs, ys = zip(*points)
                    plt.plot(xs, ys, style)

        if title is not None:
            plt.title(title)

        if show is not None:
            plt.show(fig)

        if save_path is not None:
            fig.set_size_inches(2*11.7, 2*8.27)
            if formats is None:
                plt.savefig(save_path, dpi=200, format='pdf')
            else:
                for fmt in formats:
                    plt.savefig(save_path + '.' + fmt, dpi=200, format=fmt)
            plt.close()
        return fig

    def save_to_csv(self, path):
        with open(path, 'w') as output:
            writer = csv.writer(output, delimiter=';')
            for trajectory in self.list:
                writer.writerow(trajectory)

def read_csv(path):
    csv_trajectory_list = TrajectoryList()
    csv_trajectory_list.path = path
    csv_trajectory_list.dir = os.path.dirname(path)
    with open(path) as f:
        reader = csv.reader(f, delimiter=';')
        for row in reader:
            csv_trajectory_list.add_trajectory(Trajectory.read_from_csv(row))

    return csv_trajectory_list

def normalize_trajectories_length(list_of_trajectories):
    max_length = max([len(trajectory) for trajectory in list_of_trajectories])

    for trajectory in list_of_trajectories:
        if len(trajectory) != max_length:
            for i in range(len(trajectory), max_length):
                trajectory.add_point(None)
    return max_length, list_of_trajectories

def cross_check_of_two_points(a_start, a_end, b_start, b_end):
    if a_start is None or a_end is None or b_start is None or b_end is None:
        return 0

    crossing = 0

    a_end = np.array([a_end[0], a_end[1]])
    a_beg = np.array([a_start[0], a_start[1]])

    b_end = np.array([b_end[0], b_end[1]])
    b_beg = np.array([b_start[0], b_start[1]])

    point_temp = a_end - a_beg
    beg_cross = np.cross((b_beg - a_beg), point_temp)
    end_cross = np.cross((b_end - a_beg), point_temp)
    if beg_cross != 0 and end_cross != 0 and np.sign(beg_cross) != np.sign(end_cross):
        temp2 = b_end - b_beg
        beg2_cross = np.cross((a_beg - b_beg), temp2)
        end2_cross = np.cross((a_end - b_beg), temp2)
        if beg2_cross != 0 and end2_cross != 0 and np.sign(beg2_cross) != np.sign(end2_cross):
            crossing = 1

    return crossing

def add_trajectory_to_plot_queue(trajectory, level_interval, plot_queue):
    point_style = TrajectoryList.STYLE_TRUE_PEAK if trajectory.flag else TrajectoryList.STYLE_ARTIFACT
    valid_points = [(idx, point) for idx, point in enumerate(trajectory) if idx in level_interval and point is not None]
    start_point = valid_points[0][1] if len(valid_points) > 0 else None

    not_none_point = None
    for idx, point in valid_points:
        plot_queue[point_style].append(point)
        # k = trajectory.get_next_point(idx)

        if not_none_point is None:
            not_none_point = (point, idx)
        else:
            previous_point, previous_idx = not_none_point
            if previous_idx + 1 == idx:
                plot_queue[TrajectoryList.STYLE_CONN].append([previous_point, point])
            else:
                plot_queue[TrajectoryList.STYLE_GAP_CONN].append([previous_point, point])
                for i in range(idx - previous_idx - 1):
                    middle_point_x = previous_point[0] + (i + 1) / (idx - previous_idx) * (point[0] - previous_point[0])
                    middle_point_y = previous_point[1] + (i + 1) / (idx - previous_idx) * (point[1] - previous_point[1])
                    plot_queue[TrajectoryList.STYLE_GAP_FILLER].append((middle_point_x, middle_point_y))
            not_none_point = (point, idx)

    return (plot_queue, start_point)



