import unittest
from math import isclose

import numpy as np

from code.src.porosity import Porosity
from code.src.soil_properties import SoilProperties
from code.src.tree_roots import TreeRoots
from code.src.water_content import WaterContent


class TestTreeRoots(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestTreeRoots - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestTreeRoots - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Add fields to the main object. These will be used to create
        different root density profiles with identical input values.

        :return: None.
        """

        # Grid parameters [L: cm].
        self.dz, z_max = 5.0, 1500.0

        # Test vertical domain.
        self.z_grid = np.arange(0.0, z_max + self.dz, self.dz)

        # Layers.
        self.layers = (0, 50, 200, z_max)

        # Water content object with default parameters.
        self.theta = WaterContent()

        # Soil properties object with default parameters.
        self.soil = SoilProperties()

        # Number of discrete cells forming the root zone.
        self.ln = 200

        # Test grid (0 - max_depth) [L: cm].
        self.z_roots = np.linspace(0.0, self.ln * self.dz, self.ln)

        # Create a linear porosity object.
        self.porous = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Linear')
    # _end_def_

    def test_call(self):
        """
        Test the __call__() method with all the known tree root pdf models.
        :return: None
        """

        # Make a list with all the know models.
        root_models = ["uniform", "negative_exp", "gamma_pdf", "mixture"]

        # Test all the models one at a time.
        for model_i in root_models:

            # Print info.
            print(" Testing {0} model ... ".format(model_i))

            # Create a linear porosity object.
            test_obj = TreeRoots(self.ln, self.dz, model_i, self.porous)

            # Get the full root profile.
            root_pdf_1 = test_obj()

            # Check if the integrated density sums to one.
            self.assertTrue(isclose(np.sum(root_pdf_1) * self.dz, 1.0, rel_tol=1.0e-5))

            # Make sure the max root depths match.
            self.assertEqual(test_obj.max_root_depth, (self.ln * self.dz))

            # Get the root profile on the test grid.
            root_pdf_2 = test_obj(self.z_roots)

            # Check if the densities are almost equal.
            self.assertTrue(isclose(np.sum(root_pdf_1) * self.dz,
                                    np.sum(root_pdf_2) * self.dz,
                                    rel_tol=1.0e-4))
        # _end_for_

    # _end_def_

    def test_efficiency(self):
        # Create a linear porosity object.
        test_obj = TreeRoots(self.ln, self.dz, "uniform", self.porous)

        # Print object.
        print(test_obj)

        # Zero water version.
        theta_0 = np.zeros(self.ln)
        rho_theta0, _ = test_obj.efficiency(theta_0, self.z_roots)
        self.assertTrue(np.sum(rho_theta0) == 0.0)

        # Simple version.
        theta_1d = np.random.rand(self.ln)
        rho_theta1, _ = test_obj.efficiency(theta_1d, self.z_roots)

        # Check the input/output dimensions.
        self.assertEqual(theta_1d.shape, rho_theta1.shape)

        # Check if the integrated water efficiency sums to one.
        self.assertTrue(isclose(np.sum(rho_theta1) * self.dz, 1.0, rel_tol=1.0e-5))

        # Vectorized version.
        theta_2d = np.random.rand(4, self.ln)
        rho_theta2, _ = test_obj.efficiency(theta_2d, self.z_roots)

        # Check the input/output dimensions.
        self.assertEqual(theta_2d.shape, rho_theta2.shape)

        # Check if the integrated water efficiency sums to one.
        total_sum = np.sum(rho_theta2, axis=1) * self.dz
        for sum_i in total_sum:
            # Check each dimension in the vector.
            self.assertTrue(isclose(sum_i, 1.0, rel_tol=1.0e-5))
        # _end_if_
    # _end_def_

    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            # The input model is unknown to the class.
            _ = TreeRoots(self.ln, self.dz, "Gaussian", self.porous)
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
