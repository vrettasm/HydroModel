import unittest
from code.src.soil_properties import SoilProperties

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
            self.test_obj.param_n = 0.0
        # _end_with_
    # _end_def_

    def test_a0(self):
        """
        Soil property 'a0' should be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.param_a0 = 0.0
        # _end_with_
    # _end_def_

    def test_psi_sat(self):
        """
        Soil property 'psi_sat' should be <= 0.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.psiSat = 1.0
        # _end_with_
    # _end_def_

    def test_epsilon(self):
        """
        Soil property 'epsilon' should be > 0.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.param_epsilon = 0.0
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
