import unittest
from nmr_commons.tracking.Trajectory import Trajectory

class TestTrajectoryMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.trajectory = Trajectory([(1, 1), (2, 2), None, None, (4, 4), (5, 5)])

    def test_next_point_method_normal_case(self):
        self.assertEqual(self.trajectory.get_next_point(0), 1)

    def test_next_point_method_next_none_present_point(self):
        self.assertEqual(self.trajectory.get_next_point(1), 4)

    def test_next_point_method_next_none_present_none(self):
        self.assertEqual(self.trajectory.get_next_point(2), 4)

    def test_next_point_method_last_point(self):
        self.assertEqual(self.trajectory.get_next_point(5), None)

if __name__ == '__main__':
    unittest.main()