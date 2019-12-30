import unittest
from nmr_commons.tracking.Trajectory import Trajectory
from nmr_commons.tracking.TrajectoryList import TrajectoryList


class TestTrajectoryMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.trajectory = Trajectory([(1, 1), (2, 2), None, None, (8, 8), (9, 10.5)])
        self.trajectory_complement_short = Trajectory([None, None, (3, 3), None, None])
        self.trajectory_complement_long = Trajectory([None, None, (3, 3), (4, 4), None, None, (10, 13)])
        self.trajectory_angle = Trajectory([(1, 1), (2, 2), None, (3, 3), (4, 2)])

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
        self.trajectorylist = TrajectoryList(list_of_trajectory=[self.trajectory,
                                                                 self.trajectory_complement_long,
                                                                 self.trajectory_complement_short,
                                                                 self.trajectory_angle])
    def test_len_method(self):
        self.assertEqual(len(self.trajectorylist), 4)

    def test_getitem(self):
        self.assertEqual(self.trajectorylist[3][4], (4, 2))



if __name__ == '__main__':
    unittest.main()