import unittest
import numpy as np
from math import isclose
from code.src.root_density import RootDensity

class TestRootDensity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestRootDensity - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestRootDensity - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Add fields to the main object. These will be used to create
        different root density profiles with identical input values.

        :return: None.
        """

        # Select randomly a space discretization step [L: cm].
        self.dz = np.random.randint(1, 10)

        # Number of discrete cells forming the root zone.
        self.ln = 200

        # Test grid (0 - max_depth) [L: cm].
        self.z_grid = np.linspace(0.0, self.ln * self.dz, self.ln)
    # _end_def_

    def test_call(self):
        """
        Test the __call__() method with all the known root pdf models.
        :return: None
        """

        # Make a list with all the know models.
        root_models = ["uniform", "negative_exp", "gamma_pdf", "mixture"]

        # Test all the models one at a time.
        for model_i in root_models:

            # Print info.
            print(" Testing {0} model ... ".format(model_i))

            # Create a linear porosity object.
            test_obj = RootDensity(self.ln, self.dz, model_i)

            # Get the full root profile.
            root_pdf_1 = test_obj()

            # Check if the integrated density sums to one.
            self.assertTrue(isclose(np.sum(root_pdf_1) * self.dz, 1.0, rel_tol=1.0e-5))

            # Make sure the max root depths match.
            self.assertEqual(test_obj.max_depth_m, (self.ln * self.dz)/100.0)

            # Get the root profile on the test grid.
            root_pdf_2 = test_obj(self.z_grid)

            # Check if the densities are almost equal.
            self.assertTrue(isclose(np.sum(root_pdf_1) * self.dz,
                                    np.sum(root_pdf_2) * self.dz,
                                    rel_tol=1.0e-4))
        # _end_for_

    # _end_def_

    # _end_def_
    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            # The input model is unknown to the class.
            _ = RootDensity(self.ln, self.dz, "Gaussian")
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
