import os
import sys
import unittest

# Make sure we can import the /src.
sys.path.append(os.path.abspath("../../code"))

from src.soil_properties import SoilProperties


class TestSoilProperties(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestSoilProperties - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestSoilProperties - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Creates a test object with the default input parameters.
        :return: None.
        """
        self.test_obj = SoilProperties()
    # _end_def_

    def test_n(self):
        """
        Soil property 'n' should be > 1.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.n = 0.0
        # _end_with_
    # _end_def_

    def test_alpha(self):
        """
        Soil property 'alpha' should be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.alpha = 0.0
        # _end_with_
    # _end_def_

    def test_psi_sat(self):
        """
        Soil property 'psi_sat' should be <= 0.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.psi_sat = 1.0
        # _end_with_
    # _end_def_

    def test_epsilon(self):
        """
        Soil property 'epsilon' should be > 0.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.epsilon = 0.0
        # _end_with_
    # _end_def_

    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            _ = SoilProperties(n=0.0)
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
