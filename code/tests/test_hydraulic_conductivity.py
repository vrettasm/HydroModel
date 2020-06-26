import unittest

from code.src.hydraulic_conductivity import HydraulicConductivity


class TestHydraulicConductivity(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestHydraulicConductivity - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestHydraulicConductivity - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Creates a test object with the default input parameters.
        :return: None.
        """
        self.test_obj = HydraulicConductivity()
    # _end_def_

    def test_sat_soil(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.sat_soil = 0.0
        # _end_with_
    # _end_def_

    def test_sat_saprolite(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.sat_saprolite = 0.0
        # _end_with_
    # _end_def_

    def test_sat_fresh_bedrock(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.sat_fresh_bedrock = 0.0
        # _end_with_
    # _end_def_

    def test_sigma_noise(self):
        """
        Noise model amplitude parameters can't be negative.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.sigma_noise = -1.0
        # _end_with_
    # _end_def_

    def test_lambda_exponent(self):
        """
        Noise model exponent parameters can't be negative.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.lambda_exponent = -1.0
        # _end_with_
    # _end_def_

    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            # Saturated values should be strictly positive.
            _ = HydraulicConductivity(sat_soil=0.0)
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
