import unittest
import numpy as np
from code.src.porosity import Porosity
from code.src.water_content import WaterContent
from code.src.soil_properties import SoilProperties


class TestPorosity(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestPorosity - START -")
        cls.ttime = 'noe'
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestPorosity - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Creates test objects with default input parameters.
        These will be used to create different porosity profiles.

        :return: None.
        """
        # Grid parameters [L: cm].
        dz, z_max = 5.0, 1500.0

        # Test vertical domain.
        self.z_grid = np.arange(0.0, z_max + dz, dz)

        # Layers.
        self.layers = (0, 50, 200, z_max)

        # Water content object with default parameters.
        self.theta = WaterContent()

        # Soil properties object with default parameters.
        self.soil = SoilProperties()
    # _end_def_

    def test_call_linear(self):
        """
        Test the 'Linear' Porosity profile.
        :return: None
        """
        # Create a linear porosity object.
        test_obj = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Linear')

        # Print the object.
        print(test_obj)

        # Get all the profiles.
        poros_0, field_cap_0, wilting_0 = test_obj()

        # Make sure the shapes are equal.
        self.assertEqual(self.z_grid.shape, poros_0.shape)
        self.assertEqual(poros_0.shape, field_cap_0.shape)
        self.assertEqual(field_cap_0.shape, wilting_0.shape)

        # Make sure the output sizes match.
        poros_1, field_cap_1, wilting_1 = test_obj(101.4)
        self.assertEqual(1, poros_1.size)
        self.assertEqual(poros_1.size, field_cap_1.size)
        self.assertEqual(field_cap_1.size, wilting_1.size)

        # Check Max limit.
        self.assertTrue(np.all(self.theta.max >= poros_0))

        # Check Min limit.
        self.assertTrue(np.all(self.theta.min <= poros_0))
    # _end_def_

    def test_call_constant(self):
        """
        Test the 'Constant' Porosity profile.
        :return: None
        """
        # Create a linear porosity object.
        test_obj = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Constant')

        # Print the object.
        print(test_obj)

        # Get all the profiles.
        poros_0, field_cap_0, wilting_0 = test_obj()

        # Make sure the shapes are equal.
        self.assertEqual(self.z_grid.shape, poros_0.shape)
        self.assertEqual(poros_0.shape, field_cap_0.shape)
        self.assertEqual(field_cap_0.shape, wilting_0.shape)

        # Make sure the output sizes match.
        poros_1, field_cap_1, wilting_1 = test_obj(101.4)
        self.assertEqual(1, poros_1.size)
        self.assertEqual(poros_1.size, field_cap_1.size)
        self.assertEqual(field_cap_1.size, wilting_1.size)

        # Check Max limit.
        self.assertTrue(np.all(self.theta.max >= poros_0))

        # Check Min limit.
        self.assertTrue(np.all(self.theta.min <= poros_0))
    # _end_def_

    def test_call_exp(self):
        """
        Test the 'Exponential' Porosity profile.
        :return: None
        """
        # Create a linear porosity object.
        test_obj = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Exponential')

        # Print the object.
        print(test_obj)

        # Get all the profiles.
        poros_0, field_cap_0, wilting_0 = test_obj()

        # Make sure the shapes are equal.
        self.assertEqual(self.z_grid.shape, poros_0.shape)
        self.assertEqual(poros_0.shape, field_cap_0.shape)
        self.assertEqual(field_cap_0.shape, wilting_0.shape)

        # Make sure the output sizes match.
        poros_1, field_cap_1, wilting_1 = test_obj(101.4)
        self.assertEqual(1, poros_1.size)
        self.assertEqual(poros_1.size, field_cap_1.size)
        self.assertEqual(field_cap_1.size, wilting_1.size)

        # Check Max limit.
        self.assertTrue(np.all(self.theta.max >= poros_0))

        # Check Min limit.
        self.assertTrue(np.all(self.theta.min <= poros_0))
    # _end_def_

    def test_call_stratified(self):
        """
        Test the 'Stratified' Porosity profile.
        :return: None
        """
        # Create a linear porosity object.
        test_obj = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Stratified')

        # Print the object.
        print(test_obj)

        # Get all the profiles.
        poros_0, field_cap_0, wilting_0 = test_obj()

        # Make sure the shapes are equal.
        self.assertEqual(self.z_grid.shape, poros_0.shape)
        self.assertEqual(poros_0.shape, field_cap_0.shape)
        self.assertEqual(field_cap_0.shape, wilting_0.shape)

        # Make sure the output sizes match.
        poros_1, field_cap_1, wilting_1 = test_obj(101.4)
        self.assertEqual(1, poros_1.size)
        self.assertEqual(poros_1.size, field_cap_1.size)
        self.assertEqual(field_cap_1.size, wilting_1.size)

        # Check Max limit.
        self.assertTrue(np.all(self.theta.max >= poros_0))

        # Check Min limit.
        self.assertTrue(np.all(self.theta.min <= poros_0))
    # _end_def_

    def test_call_noisy(self):
        """
        Test the 'Noisy' Porosity profile.
        :return: None
        """
        # Create a linear porosity object.
        test_obj = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Noisy')

        # Print the object.
        print(test_obj)

        # Get all the profiles.
        poros_0, field_cap_0, wilting_0 = test_obj()

        # Make sure the shapes are equal.
        self.assertEqual(self.z_grid.shape, poros_0.shape)
        self.assertEqual(poros_0.shape, field_cap_0.shape)
        self.assertEqual(field_cap_0.shape, wilting_0.shape)

        # Make sure the output sizes match.
        poros_1, field_cap_1, wilting_1 = test_obj(101.4)
        self.assertEqual(1, poros_1.size)
        self.assertEqual(poros_1.size, field_cap_1.size)
        self.assertEqual(field_cap_1.size, wilting_1.size)

        # Check Max limit.
        self.assertTrue(np.all(self.theta.max >= poros_0))

        # Check Min limit.
        self.assertTrue(np.all(self.theta.min <= poros_0))
    # _end_def_

    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            # Unknown value for the profile is given ...
            _ = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Unknown')
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
