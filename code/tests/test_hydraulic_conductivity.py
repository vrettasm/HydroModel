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

    def test_satSoil(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.satSoil = 0.0
        # _end_with_
    # _end_def_

    def test_satSaprolite(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.satSaprolite = 0.0
        # _end_with_
    # _end_def_

    def test_satFreshBedrock(self):
        """
        Saturated values must be strictly positive.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.satFreshBedrock = 0.0
        # _end_with_
    # _end_def_

    def test_sigmaAmp(self):
        """
        Noise model amplitude parameters can't be negative.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.sigmaAmp = -1.0
        # _end_with_
    # _end_def_

    def test_lambdaExp(self):
        """
        Noise model exponent parameters can't be negative.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.lambdaExp = -1.0
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
