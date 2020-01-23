import unittest
from numpy.testing import assert_almost_equal

from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.tracking.TrajectoryList import TrajectoryList
from nmr_commons.tracking.TrejectoryGenerator import TrajectoryGenerator, get_path_list
from nmr_environment import Settings


class TestTrajectoryMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.trajectory = Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5)])
        self.trajectory_complement_short = Trajectory([None, None, (3, 3), None, None])
        self.trajectory_complement_long = Trajectory([None, None, (3, 3), (4, 4), None, None, (10, 13)])
        self.trajectory_angle = Trajectory([(1, 1), (2, 2), None, (3, 3), (4, 2)])

    def test_getitem(self):
        self.assertEqual(self.trajectory[4], (8, 8))

    def test_eq_method(self):
        self.assertEqual(self.trajectory, Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5)]))

    def test_get_x_method_normal_case(self):
        self.assertEqual(self.trajectory.get_x(5), 9)

    def test_get_x_method_none_case(self):
        self.assertEqual(self.trajectory.get_x(3), None)

    def test_get_y_method_normal_case(self):
        self.assertEqual(self.trajectory.get_y(5), 10.5)

    def test_get_y_method_none_case(self):
        self.assertEqual(self.trajectory.get_x(3), None)

    def test_get_vector_method_normal_case(self):
        self.assertEqual(self.trajectory.get_vector(0), 0)

    def test_get_vector_method_none_between(self):
        self.assertEqual(self.trajectory.get_vector(1), 4)

    def test_get_vector_method_last_point(self):
        self.assertEqual(self.trajectory.get_vector(5), None)

    def test_next_point_method_normal_case(self):
        self.assertEqual(self.trajectory.get_next_point(0), 1)

    def test_next_point_method_next_none_present_point(self):
        self.assertEqual(self.trajectory.get_next_point(1), 4)

    def test_next_point_method_next_none_present_none(self):
        self.assertEqual(self.trajectory.get_next_point(2), 4)

    def test_next_point_method_last_point(self):
        self.assertEqual(self.trajectory.get_next_point(5), None)

    def test_get_distances(self):
        self.assertEqual(self.trajectory.get_distances(),
                         [(-1 * x[0], -1 * x[1]) for x in [(1, 1), (2, 2), (2, 2), (2, 2), (1, 2.5)]])

    def test_get_distances_sep(self):
        self.assertEqual(self.trajectory.get_distances_sep(),
                         [(-1 * x[0], -1 * x[1]) for x in [(1, 1), (1, 2.5)]])

    def test_get_angles(self):
        self.assertAlmostEqual(self.trajectory_angle.get_angles()[0], 1.57, places=2)

    def test_get_angles(self):
        self.assertAlmostEqual(self.trajectory.get_angles()[0], 0.4, places=2)

    def test_get_all_x(self):
        self.assertEqual(self.trajectory.get_all_x(), [1, 2, None, None, 8, 9])

    def test_get_all_x(self):
        self.assertEqual(self.trajectory.get_all_x(with_none_points=False), [1, 2, 8, 9])

    def test_get_all_y(self):
        self.assertEqual(self.trajectory.get_all_y(), [1, 2, None, None, 8, 10.5])

    def test_get_mean(self):
        self.assertEqual(self.trajectory.get_mean(), [5, 5.375])

    def test_can_marge_true_short(self):
        self.assertEqual(self.trajectory.can_merge(self.trajectory_complement_short), True)

    def test_can_marge_true_long(self):
        self.assertEqual(self.trajectory.can_merge(self.trajectory_complement_long), True)

    def test_can_marge_false(self):
        self.assertEqual(self.trajectory.can_merge(self.trajectory_angle), False)

    def test_marge_short(self):
        self.trajectory.merge(self.trajectory_complement_short)
        self.assertEqual(self.trajectory[:], [(1, 1), (2, 2), (3, 3), None, (8, 8), (9, 10.5)])

    def test_marge_long(self):
        self.trajectory.merge(self.trajectory_complement_long)
        self.assertEqual(self.trajectory[:], [(1, 1), (2, 2), (3, 3), (4, 4), (8, 8), (9, 10.5)])

    def test_get_first_last_normal_calse(self):
        self.assertEqual(self.trajectory.get_first_last(), (0, 5, (1, 1), (9, 10.5)))

    def test_get_first_last_none_in_the_beggining(self):
        self.assertEqual(self.trajectory_complement_long.get_first_last(), (2, 6, (3, 3), (10, 13)))

    def test_get_first_last_only_one_element(self):
        self.assertEqual(self.trajectory_complement_short.get_first_last(), (2, 2, (3, 3), (3, 3)))

    def test_get_not_none_count_normal(self):
        self.assertEqual(self.trajectory.get_not_none_count(), 4)

    def test_get_path_ratio_sum_dist(self):
        self.assertAlmostEqual(self.trajectory.get_path_ratio()[1], 1.6013, places=4)

    def test_get_path_ratio_first_last_dist(self):
        self.assertAlmostEqual(self.trajectory.get_path_ratio()[2], 1.5997, places=4)


class TestTrajectoryListMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.trajectory = Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5)])
        self.trajectory_complement_short = Trajectory([None, None, (3, 3), None, None])
        self.trajectory_complement_long = Trajectory([None, None, (3, 3), (4, 4), None, None, (10, 13)])
        self.trajectory_angle = Trajectory([(1, 1), (2, 2), None, (3, 3), (4, 2)])
        self.trajectory_list = TrajectoryList(list_of_trajectories=[self.trajectory,
                                                                    self.trajectory_complement_short,
                                                                    self.trajectory_complement_long,
                                                                    self.trajectory_angle])

        self.trajectory_for_crossing1 = Trajectory([(0, 0), (1, 1)][::-1])
        self.trajectory_for_crossing2 = Trajectory([(1, 0), (0, 1)])
        self.trajectory_for_crossing3 = Trajectory([(0, 0), None, (5, 5)])
        self.trajectory_for_crossing4 = Trajectory([(0, 0), (1, 0), (2, 1), (2, 2), (1, 3), (0, 3)][::-1])
        self.trajectory_for_crossing5 = Trajectory([(3, 0), (2, 0), (1, 1), (1, 2), (2, 3), (3, 3)])
        self.trajectory_list_for_crossing = TrajectoryList(list_of_trajectories=[self.trajectory_for_crossing1,
                                                                                 self.trajectory_for_crossing2,
                                                                                 self.trajectory_for_crossing3,
                                                                                 self.trajectory_for_crossing4,
                                                                                 self.trajectory_for_crossing5])

    def test_normalize_trajecotories_length(self):
        self.assertEqual([len(trajectory) for trajectory in self.trajectory_list], [7] * 4)
        self.assertEqual(self.trajectory_list.number_of_levels, 7)

    def test_len_method(self):
        self.assertEqual(len(self.trajectory_list), 4)

    def test_getitem(self):
        self.assertEqual(self.trajectory_list[3][4], (4, 2))

    def test_getitem_all(self):
        self.assertEqual(self.trajectory_list[2], Trajectory([None, None, (3, 3), (4, 4), None, None, (10, 13)]))

    def test_getitem_small_list(self):
        small_trajectory_list = self.trajectory_list[[0, 3]]
        self.assertEqual(small_trajectory_list[0],
                         Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5), None]))
        self.assertEqual(small_trajectory_list[1],
                         Trajectory([(1, 1), (2, 2), None, (3, 3), (4, 2), None, None]))

    def test_getitem_slice(self):
        small_trajectory_list = self.trajectory_list[:3:2]
        self.assertEqual(small_trajectory_list[0], Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5), None]))
        self.assertEqual(small_trajectory_list[1], Trajectory([None, None, (3, 3), (4, 4), None, None, (10, 13)]))

    def test_get_bounding_box_of_some_trajectories_idx_list(self):
        self.assertListEqual(self.trajectory_list[[0, 3]].get_trajectories_bounding_box(),
                             [(1, 9, 1.0, 10.5), (1, 4, 1, 3)])

    def test_get_bounding_box_of_some_trajectories_idx_slice(self):
        self.assertListEqual(self.trajectory_list[0:3:2].get_trajectories_bounding_box(),
                             [(1, 9, 1.0, 10.5), (3, 10, 3, 13)])

    def test_mean_of_points(self):
        self.trajectory_list.calculate_min_max()
        self.assertEqual(self.trajectory_list.min_x, 1)
        self.assertEqual(self.trajectory_list.max_x, 10)
        self.assertEqual(self.trajectory_list.min_y, 1)
        self.assertEqual(self.trajectory_list.max_y, 13)

    def test_mean_of_points(self):
        self.trajectory_list.calculate_mean_of_points()
        self.assertAlmostEqual(self.trajectory_list.mean_y, 4.375, places=3)

    def test_normalization(self):
        self.trajectory_list.normalize()
        self.assertAlmostEqual(self.trajectory_list.norm_list[0][1][0], 0.11, places=2)
        self.assertAlmostEqual(self.trajectory_list.norm_list[0][1][1], 0.083, places=3)
        self.assertAlmostEqual(self.trajectory_list.norm_list[0][2], None, places=3)
        self.assertAlmostEqual(self.trajectory_list.norm_list[0][5][0], 0.89, places=2)
        self.assertAlmostEqual(self.trajectory_list.norm_list[0][5][1], 0.79, places=2)
        self.assertEqual(len(self.trajectory_list.norm_list), 4)

    def test_get_trajectory(self):
        self.assertEqual(self.trajectory_list.get_trajectory(0)[:],
                         Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5), None])[:])

    def test_number_of_spectra_sequence(self):
        self.assertEqual(self.trajectory_list.number_of_levels, 7)

    def test_get_distances(self):
        self.assertListEqual(self.trajectory_list.get_distances_euclidean(interval=[0, 3]),
                             [[1.4142135623730951, 0.0, 0.0, 1.118033988749895],
                              [0.7071067811865476, 0.0, 1.5811388300841898]])

    def test_get_distances(self):
        self.assertListEqual(self.trajectory_list.get_distances_euclidean(),
                             [[1.4142135623730951, 0.0, 0.0, 1.118033988749895],
                              [],
                              [2.23606797749979, 0.0, 0.0],
                              [0.7071067811865476, 0.0, 1.5811388300841898]])

    def test_get_distances(self):
        self.assertListEqual(self.trajectory_list.get_distances_min_max_norm(interval=[0, 3]),
                             [[(-0.1111111111111111, -0.08333333333333333),
                               (-0.1111111111111111, -0.20833333333333334)],
                              [(-0.1111111111111111, -0.08333333333333333),
                               (-0.1111111111111111, 0.08333333333333333)]])

    def test_get_distances(self):
        self.assertListEqual(self.trajectory_list.get_mean_angles_of_trajectory(),
                             [0.4048917862850834, None, None, 1.5707963267948966])

    def test_get_distances(self):
        self.assertListEqual(self.trajectory_list.get_distances_first_point_from_mean(interval=[0, 3]),
                             [[4, 4.375], [1.5, 1]])

    def test_get_distances_from_start_to_end(self):
        self.assertListEqual(self.trajectory_list.get_distances_from_start_to_end(),
                             [[8, 9.5], [0, 0], [7, 10], [3, 1]])

    def test_trajectories_moving(self):
        self.assertListEqual(self.trajectory_list.trajectories_moving_staying_idx(), [[0, 2, 3], [1]])

    def test_get_crossing_trajectories(self):
        self.assertListEqual(self.trajectory_list_for_crossing.get_crossing_trajectories_idx(),
                             [(0, 1), (1, 2), (3, 4), (3, 4)])

    # def test_visualization_3d(self):
    #     self.trajectory_list_for_crossing.visualize_3d(show=True)

    # def test_visualization_2d(self):
    #     self.trajectory_list_for_crossing.visualize_2d(show=True)


class TestTrajectoryGeneratorMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.trajectory_generator = TrajectoryGenerator()
        self.trajectory_generator.list_of_trajectory_list = self.trajectory_generator.get_list_of_trajectory_list()

    def test_get_path_list(self):
        self.assertListEqual(get_path_list('/Users/zl449wx/projects/NMR-new/NMR/Datasets/Tracking'),
                             Settings.TRAJECTORIES_PATH)

    def test_get_list_of_trajectory_list(self):
        self.assertEqual(self.trajectory_generator.get_list_of_trajectory_list()[0][0],
                         Trajectory([(111.394, 3.571), (111.394, 3.571), (111.395, 3.572), (111.398, 3.572),
                                     (111.399, 3.573), (111.4, 3.573), (111.402, 3.574), (111.404, 3.574),
                                     (111.406, 3.575), (111.408, 3.575), (111.41, 3.575)]))
        self.assertEqual(self.trajectory_generator.get_list_of_trajectory_list()[-1][-1],
                         Trajectory([(113.205, 4.405), (113.197, 4.399), (113.198, 4.399), (113.199, 4.395),
                                     (113.185, 4.396), (113.19, 4.394), None, None, None]))

    def test_caltulation_of_lambda(self):
        self.assertAlmostEqual(self.trajectory_generator.calculate_lambda(), 0.3216, places=4)

    def test_caltulation_of_noise_covariance(self):
        self.assertAlmostEqual(self.trajectory_generator.calculate_noise_cov()[0][0], 4.892e-07, places=10)
        self.assertAlmostEqual(self.trajectory_generator.calculate_noise_cov()[1][1], 4.802e-07, places=10)


if __name__ == '__main__':
    unittest.main()
