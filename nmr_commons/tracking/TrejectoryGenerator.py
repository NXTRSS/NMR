import numpy as np
from os import listdir, walk
from os.path import isfile, join
import re
import operator
import matplotlib.pyplot as plt
import numbers
from scipy.stats import rv_continuous, gamma, cauchy, norm, kstest
import math
import warnings
import numpy.linalg as la
import scipy

from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.tracking.TrajectoryList import TrajectoryList
from nmr_environment import Settings
from nmr_commons.utils import flatten
from chemical_shift_assignment.components.ChemicalShiftGenerator import ChemicalShiftGenerator
from nmr_commons.database.BMRBDatabase import BMRBDatabase
from nmr_commons.peak_picking.PeakList import PeakList
from spectrum_generator.components.PeakListGenerator import PeakListGenerator as PLGen, PeakListGenParams
from nmr_commons.spectrum.Spectrum import Spectrum


class TrajectoryGenerator:
    """Class for training and generating artificial list of Trajectories. This is a main class in data augmentation task
    for NMR protein purpose."""

    GAMMA_DISTRIBUTION_PARAMETER_FOR_VELOCITY = 2
    GAMMA_DISTRIBUTION_SHIFTING_PARAMETER_FOR_VELOCITY = 0
    MINIMAL_AVERAGE_ANGLE_IN_ONE_TRAJECTORY = 0.01 * np.pi * 2
    DROP_QUANTILE_FROM_ANGLES = 0.2
    DROP_QUANTILE_OF_CURVING_TRAJECTORIES_RATIO = 0.2

    ####
    # Constructor
    def __init__(self, number_of_levels=None, path_list=None, list_of_trajectory_list=None,
                 curving_lim=0.12, k0=5, ksi0=0.0002):
        self.levels_number = number_of_levels
        self.lam = 0  # ratio of moving trajectories vs all trajectories in trajectories lists
        self.epsilon = 0
        self.mi = [0, 0]
        self.ni = 0
        self.r = 0
        self.k0 = k0
        self.ksi0 = ksi0
        self.ksi = 0
        self.Ibeta = [[0, 0], [0, 0]]  # covariance matrix of noise. Noise calculated only on staying trajectories
        self.sigma = 0
        self.r = []
        self.path_list = path_list
        if list_of_trajectory_list is None:
            self.list_of_trajectory_list = []
        else:
            self.list_of_trajectory_list = list_of_trajectory_list
        self.angle_dist = None
        self.angle_dist_params = [0, 0]
        self.curving_ratio_dist = None
        self.curving_ratio_params = [0, 0]
        self.gamma = 0
        # self.curving_lim = curving_lim
        self.chem_shift_gen = None

    def __iter__(self):
        return self.list_of_trajectory_list.__iter__()

    def __getitem__(self, index):
        cls = type(self)
        if isinstance(index, slice):
            return cls(list_of_trajectory_list=self.list_of_trajectory_list[index])
        elif isinstance(index, list) or isinstance(index, range):
            if isinstance(index, range):
                index = list(index)
            return cls(list_of_trajectory_list=[self[idx] for idx in index])
        elif isinstance(index, numbers.Integral):
            return self.list_of_trajectory_list[index]
        else:
            msg = '{cls.__name__} indices must be integers'
            raise TypeError(msg.format(cls=cls))

    def __len__(self):
        return len(self.list_of_trajectory_list)

    def __eq__(self, other):
        return tuple(self) == tuple(other)

    def train(self, verbose=False):
        if not self.list_of_trajectory_list:
            self.list_of_trajectory_list = self.get_list_of_trajectory_list(verbose)
        self.lam = self.calculate_lambda()
        lam = self.lam

        self.Ibeta = self.calculate_noise_cov()
        self.sigma, self.k0, self.ksi0 = self.calculate_velocity_parameters()
        self.angle_dist, self.angle_dist_params = self.calculate_angle_noise_old()
        if verbose:
            print('Lambda: {:.4f}, k0: {:.4f}, ksi0: {:.4f}'.format(lam, self.k0, self.ksi0))
        # self.k0 = 5
        # self.ksi0 = 0.0002

    def generate(self, bmrb_id, peak_gen_params=None, extra_shifts_distr=None):
        """

        Args:
            bmrb_id: int, BMRB protein file ID
            peak_gen_params: PeakListGenParams object
            extra_shifts_distr: one of EmpiricalDistribution.{LAPLACE, GAUSS, UNIFORM, ALWAYS_MEDIAN, ALWAYS_MEAN}
        Returns:

        """
        if peak_gen_params is None:
            peak_gen_params = PeakListGenParams()
        else:
            peak_gen_params = peak_gen_params.copy()

        shifts = BMRBDatabase.parse_bmrb_file('{}bmr{}.str'.format(Settings.BMRB_ROOT, bmrb_id))

        temp_params = PeakListGenParams().incomplete(peak_gen_params.incompleteProb)
        peaks, _, ideal_num_peaks = PLGen.generate_peak_list(shifts, Spectrum.NHSQC, temp_params, verbose=False)
        peaks.set_responses(1)

        peak_gen_params.incomplete(1)
        if self.chem_shift_gen is None:
            self.chem_shift_gen = ChemicalShiftGenerator(Settings.CHEMICAL_SHIFT_GENERATOR_PATH)

        for amino_acid in shifts.sequence:
            if amino_acid.getOneLetterCode() == 'U':
                amino_acid.chemicalShifts = []
        extra_shifts = [self.chem_shift_gen.generateShiftsForSequence(shifts, extra_shifts_distr)
                        for _ in range(peak_gen_params.extraTrueMax)]
        extra_peaks, _, _ = PLGen.generate_peak_list(shifts, Spectrum.NHSQC, peak_gen_params, extra_shifts,
                                                     verbose=False)
        extra_peaks.set_responses(0)
        all_peaks = PeakList.get_concatenation(peaks, extra_peaks)
        print('\tPeaks: {:>4} ({:>6.2f}% of original)'.format(len(all_peaks), 100. * len(all_peaks) / ideal_num_peaks),
              end=' ')
        res = self.generate_trajectories(all_peaks)
        return res

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
        for trajectory_list in self:
            n += len(trajectory_list.get_moving_trajectories())
            N += len(trajectory_list)
        return float(n) / N

    def calculate_noise_cov(self):
        distances = list()
        for trajectory_list in self:
            staying = trajectory_list.get_staying_trajectories()
            distances += trajectory_list.get_distances_min_max_norm(interval=staying)
        distances = flatten(distances)
        var = np.var(distances, axis=0)

        return [[var[0], 0], [0, var[1]]]

    def calculate_velocity_parameters(self):
        ksi = []
        var = []
        for trajectory_idx, trajectory_list in enumerate(self):
            moving = trajectory_list.get_moving_trajectories()
            distances = trajectory_list.get_distances_euclidean(interval=moving)
            var += calculate_velocity_noise(distances)
            distances = flatten(distances)
            try:
                (alpha, loc, beta) = rv_continuous.fit(gamma, distances,
                                                       floc=self.GAMMA_DISTRIBUTION_SHIFTING_PARAMETER_FOR_VELOCITY,
                                                       f0=self.GAMMA_DISTRIBUTION_PARAMETER_FOR_VELOCITY)
                ksi.append(beta)
            except:
                print('Problem with fitting gamma distribution (velocity) in seq: %s' % self.path_list[trajectory_idx])

        sigma = np.mean([x for x in var if not math.isnan(x)])
        fit_alpha, fit_loc, fit_beta = rv_continuous.fit(gamma, ksi,
                                                         floc=self.GAMMA_DISTRIBUTION_SHIFTING_PARAMETER_FOR_VELOCITY)
        k0 = fit_alpha
        ksi0 = fit_beta
        return sigma, k0, ksi0

    def calculate_angle_noise_old(self, verbose=None):
        angles = []
        for trajectory_idx, trajectory_list in enumerate(self):
            trajectory_list.normalize()
            if trajectory_list.number_of_levels >= 5:
                moving = trajectory_list.get_moving_trajectories()
                levels = trajectory_list.number_of_levels

                for trajectory_id in moving:
                    norm_trajectory = trajectory_list.norm_list[trajectory_id]
                    start = 0
                    start_bool = 0
                    while start < levels - 4:
                        a = list(filter(None, norm_trajectory[start + 0:start + 4]))
                        if len(a) == 4:
                            ax, ay = zip(*a)
                            v1 = [np.mean(ax[0:2]), np.mean(ay[0:2])]
                            v2 = [np.mean(ax[2:4]), np.mean(ay[2:4])]
                            v = tuple([v2[0] - v1[0], v2[1] - v1[1]])
                            if (v[0] != 0 or v[1] != 0):
                                start_bool = 1
                                break
                        start += 1

                    end = levels
                    end_bool = 0
                    while end > 5:
                        b = list(filter(None, norm_trajectory[end - 4:end]))
                        if len(b) == 4:
                            bx, by = zip(*b)
                            u1 = [np.mean(bx[0:2]), np.mean(by[0:2])]
                            u2 = [np.mean(bx[2:4]), np.mean(by[2:4])]
                            u = tuple([u2[0] - u1[0], u2[1] - u1[1]])
                            if (u[0] != 0 or u[1] != 0):
                                end_bool = 1
                                break
                        end -= 1

                    if start_bool and end_bool:
                        angle = angle_between(v, u)

                        angles.append((angle) / (end - start - 4))

                        # plt.plot(ax, ay, 'bo')
                        # plt.plot([v1[0], v2[0]], [v1[1], v2[1]])
                        # plt.plot(v1[0], v1[1], 'go')
                        # plt.plot(bx, by, 'ro')
                        # plt.plot([u1[0], u2[0]], [u1[1], u2[1]], 'r')
                        # plt.plot(u1[0], u1[1], 'go')
                        # plt.show()
                        # plt.plot([0,v[0]],[0,v[1]],'r')
                        # plt.plot([0,u[0]], [0,u[1]], 'g')
                        # plt.show()
                        if ((angle) / (end - start - 4)) == 0:
                            # plt.plot(ax, ay, 'bo')
                            # plt.plot([v1[0], v2[0]], [v1[1], v2[1]])
                            # plt.plot(v1[0], v1[1], 'go')
                            # plt.plot(bx, by, 'ro')
                            # plt.plot([u1[0], u2[0]], [u1[1], u2[1]], 'r')
                            # plt.plot(u1[0], u1[1], 'go')
                            # plt.show()
                            # plt.plot([0,v[0]],[0,v[1]],'r')
                            # plt.plot([0,u[0]], [0,u[1]], 'g')
                            # plt.show()
                            if verbose:
                                print(trajectory_idx, trajectory_id)

        angles = [angle % np.pi if angle > 0 else -(-angle % np.pi) for angle in angles]

        angles_array = np.asarray(angles)

        # angles = map(abs, angles)

        # filtered_angles = filter(lambda x: x > self.curving_lim or x < -self.curving_lim, angles)

        # filtered_angles = map(abs, filtered_angles)
        filtered_angles = angles
        if verbose:
            print(sum(angles_array == 0))

        filtered_angles = list(filter(lambda x: x != 0, angles))
        filtered_angles = np.sort(np.asarray(filtered_angles))

        filtered_angles = filtered_angles[int(0.05 * len(filtered_angles)):int(0.95 * len(filtered_angles))]

        # filtered_angles_shift = (filtered_angles - np.mean(filtered_angles)) / (np.std(filtered_angles))
        #
        # cauchy_params = scipy.stats.cauchy.fit(filtered_angles, floc=0)
        #
        # filtered_angles_cauchy = (filtered_angles - cauchy_params[0])/cauchy_params[1]
        #
        # X = np.linspace(-0.5, 0.5, 1000)
        # Y = scipy.stats.cauchy.pdf(X, loc=cauchy_params[0], scale=cauchy_params[1])
        # plt.plot(X, Y, 'r--')
        #
        # test_resp = scipy.stats.kstest(filtered_angles_cauchy, 'cauchy')
        #
        # print (test_resp)
        # #
        # # test_resp = scipy.stats.kstest(scipy.stats.cauchy.rvs(size=100), 'cauchy')
        # #
        # # print(test_resp)
        #
        # gauss_params = scipy.stats.norm.fit(filtered_angles, loc=0)
        # filtered_angles_gauss = (filtered_angles - gauss_params[0])/gauss_params[1]
        #
        # X = np.linspace(-0.5, 0.5, 1000)
        # Y = scipy.stats.norm.pdf(X, loc=gauss_params[0], scale=gauss_params[1])
        # plt.plot(X, Y, 'g--')
        #
        # test_resp = scipy.stats.kstest(filtered_angles_gauss, 'norm')
        #
        # print (test_resp)
        #
        # test_resp = scipy.stats.kstest(scipy.stats.norm.rvs(size=100), 'norm')
        # print(test_resp)
        # plt.show()

        DISTRIBUTIONS = [cauchy]

        mark = []

        for i, distribution in enumerate(DISTRIBUTIONS):
            # print(distribution.name)
            # test_resp = scipy.stats.kstest(filtered_angles, 'alpha')
            name = distribution.name
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')

                params = distribution.fit(filtered_angles, floc=0)
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore')

                    test_resp = kstest((filtered_angles - params[0]) / params[1], name)
                except:
                    test_resp = [0, 0]
            except:
                test_resp = [0, 0]

            mark.append([i, test_resp[0], test_resp[1]])

        mark_array = np.array(mark)
        maks = np.max(mark_array[:, -1])
        mark_v = list(mark_array[:, -1])
        p = mark_v.index(maks)

        # print(mark[p], DISTRIBUTIONS[p])

        best_dist_params = DISTRIBUTIONS[p].fit(filtered_angles, floc=0)
        X = np.linspace(-0.5, 0.5, 1000)
        Y = DISTRIBUTIONS[p].pdf(X, loc=best_dist_params[0], scale=best_dist_params[1])
        plt.plot(X, Y, 'y--')

        # print('Angles will be draw from %s distribution with loc paramter %f and scale paramter %f\n'
        #       % (DISTRIBUTIONS[p].name, best_dist_params[0], best_dist_params[1]))
        # print('KS test p-value for %s distribution was %f' % (DISTRIBUTIONS[p].name, maks))

        angle_dist = DISTRIBUTIONS[p]
        angle_dist_params = best_dist_params

        n, bins, patches = plt.hist(filtered_angles, bins=50, density=True)
        plt.close()
        return angle_dist, angle_dist_params

    def calculate_angle_noise(self, verbose=False):
        angles = []
        ratio_curving_trajectories = []
        for trajectory_idx, trajectory_list in enumerate(self):
            trajectory_list.normalize()
            if trajectory_list.number_of_levels >= 5:
                moving = trajectory_list.get_moving_trajectories()
                levels = trajectory_list.number_of_levels
                trajectory_moving_number = len(moving)
                trajectory_curving_number = 0
                for trajectory_id in moving:
                    norm_trajectory = trajectory_list.norm_list[trajectory_id]
                    one_trajectory_angles = []
                    for idx in range(levels - 2):
                        if None not in norm_trajectory[idx:idx + 3]:
                            v = (norm_trajectory.get_x(idx + 1) - norm_trajectory.get_x(idx),
                                 norm_trajectory.get_y(idx + 1) - norm_trajectory.get_y(idx))
                            u = (norm_trajectory.get_x(idx + 2) - norm_trajectory.get_x(idx + 1),
                                 norm_trajectory.get_y(idx + 2) - norm_trajectory.get_y(idx + 1))
                            one_trajectory_angles.append(angle_between(v, u))
                    if np.mean(one_trajectory_angles) > self.MINIMAL_AVERAGE_ANGLE_IN_ONE_TRAJECTORY:
                        trajectory_curving_number += 1
                        angles.extend(one_trajectory_angles)
                        # TrajectoryList(list_of_trajectories=[norm_trajectory]).visualize_2d()
                        # plt.show()
                ratio_curving_trajectories.append(trajectory_curving_number / trajectory_moving_number)

        angles = [angle % np.pi if angle > 0 else -(-angle % np.pi) for angle in angles]

        filtered_angles = np.sort(np.asarray(angles))

        lower_quantile = self.DROP_QUANTILE_FROM_ANGLES / 2
        upper_quantile = 1 - lower_quantile

        filtered_angles = filtered_angles[
                          int(lower_quantile * len(filtered_angles)):int(upper_quantile * len(filtered_angles))]

        angle_dist, angle_dist_params = fit_distribution(filtered_angles, floc=0, objective_name='angles',
                                                         verbose=verbose)

        ratio_curving_trajectories = np.sort(np.array(ratio_curving_trajectories))
        ratio_curving_trajectories_filtered = ratio_curving_trajectories[
                                              int(lower_quantile * len(ratio_curving_trajectories)):int(
                                                  upper_quantile * len(ratio_curving_trajectories))]

        curving_ratio_dist, curving_ratio_dist_params = fit_distribution(ratio_curving_trajectories_filtered,
                                                                         objective_name='curving ratio',
                                                                         verbose=verbose)

        return angle_dist, angle_dist_params, curving_ratio_dist, curving_ratio_dist_params


def fit_distribution(sample_for_fitting, floc=None, objective_name='not known', distributions=None, verbose=False):
    if distributions is None:
        distributions = [cauchy, norm]

    distribution_fitting_log = [None] * len(distributions)

    for i, distribution in enumerate(distributions):
        name = distribution.name
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
            if floc is None:
                params = distribution.fit(sample_for_fitting)
            else:
                params = distribution.fit(sample_for_fitting, floc=floc)
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')

                test_resp = kstest((sample_for_fitting - params[0]) / params[1], name)
            except:
                print('KS-test for {} distribution failed, {} case'.format(name, objective_name))
                test_resp = [None, None]
        except:
            print('Fitting distribution {} for angles failed, {} case'.format(name, objective_name))
            test_resp = [None, None]
        if None not in test_resp and verbose:
            plt.hist(sample_for_fitting, density=True)
            X = np.linspace(-1.5, 1.5, 1000)
            Y = distribution.pdf(X, loc=params[0], scale=params[1])
            plt.plot(X, Y, 'y--')
            plt.title('Histogram of {} and fitted {} distribution'.format(objective_name, name))
            plt.show()

        distribution_fitting_log[i] = (test_resp[0], test_resp[1])

    distribution_array = np.array(distribution_fitting_log)
    max_p_value = np.max(distribution_array[:, -1])
    p_values = list(distribution_array[:, -1])
    p = p_values.index(max_p_value)

    if floc is None:
        best_dist_params = distributions[p].fit(sample_for_fitting)
    else:
        best_dist_params = distributions[p].fit(sample_for_fitting, floc=0)

    if verbose:
        plt.hist(sample_for_fitting, bins=10, density=True)
        X = np.linspace(-1.5, 1.5, 1000)
        Y = distributions[p].pdf(X, loc=best_dist_params[0], scale=best_dist_params[1])
        plt.plot(X, Y, 'y--')
        plt.title('Histogram of {} and fitted {} as best distribution'.format(objective_name, distributions[p].name))
        plt.show()

        n, bins, patches = plt.hist(sample_for_fitting, bins=50, density=True)
        plt.close()

        print('{} will be draw from {} distribution with loc parameter {} and scale parameter {}\n'
              .format(objective_name, distributions[p].name, best_dist_params[0], best_dist_params[1]))
        print('KS test p-value for {} distribution was {}'.format(distributions[p].name, max_p_value))

    return distributions[p], best_dist_params


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang1 = np.dot((1, 0), v1)
    sinang1 = la.norm(np.cross((1, 0), v1))

    cosang2 = np.dot((1, 0), v2)
    sinang2 = la.norm(np.cross((1, 0), v2))
    return np.arctan2(sinang1, cosang1) - np.arctan2(sinang2, cosang2)


def calculate_velocity_noise(distances):
    var = []
    for trajectory_distances in distances:
        trajectory_distances_sum = sum(trajectory_distances)
        norm_trajectory_distances = [x / trajectory_distances_sum for x in trajectory_distances]
        var.append(np.var(norm_trajectory_distances))
    return var


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
