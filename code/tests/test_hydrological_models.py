import os
import sys
import unittest
import numpy as np
from math import isclose

# Make sure we can import the /src.
sys.path.append(os.path.abspath("../../code"))

from src.porosity import Porosity
from src.water_content import WaterContent
from src.soil_properties import SoilProperties
from src.models.vrettas_fung import VrettasFung
from src.models.vanGenuchten import vanGenuchten
from src.models.hydrological_model import HydrologicalModel
from src.hydraulic_conductivity import HydraulicConductivity


class TestHydrologicalModels(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestHydrologicalModels - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestHydrologicalModels - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Creates test objects with default input parameters.
        These will be used to create different hydrological models.

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

        # Hydraulic conductivity object with default parameters.
        self.k_hc = HydraulicConductivity()

        # Create a linear porosity object.
        self.porous = Porosity(self.z_grid, self.layers, self.theta, self.soil, 'Noisy')
    # _end_def_

    def test_base_model(self):
        """
        Test parent (base) class model.

        :return: None
        """
        print(" Testing HydrologicalModel (parent class).")

        # Create a 'test-model' with default parameters.
        test_model = HydrologicalModel(self.soil, self.porous, self.k_hc, self.theta.res, self.dz)

        # Create a 'hypothetical' vol. water content.
        theta_z = np.linspace(self.theta.max, self.theta.min, self.z_grid.size)

        # Get the equivalent pressure head and effective saturation.
        psi_z, s_eff_z = test_model.pressure_head(theta_z, self.z_grid)

        # Check the dimensions.
        self.assertEqual(theta_z.shape, psi_z.shape)
        self.assertEqual(theta_z.shape, s_eff_z.shape)

        # Test 2D version.
        theta_2d = theta_z.repeat(4).reshape(theta_z.size, 4)

        # Get the equivalent pressure head and effective saturation.
        psi_2d, *_ = test_model.pressure_head(theta_2d, self.z_grid)

        # Check the dimensions.
        self.assertEqual(psi_2d.shape, theta_2d.shape)
    # _end_def_

    def test_vanG_model(self):
        """
        Test vanGenuchten class model.

        :return: None
        """
        print(" Testing vanGenuchten (child class).")

        # Create a 'test-model' with default parameters.
        test_model = vanGenuchten(self.soil, self.porous, self.k_hc, self.theta.res, self.dz)

        # Create a 'hypothetical' vol. water content.
        theta_1d = np.linspace(self.theta.max, self.theta.min, self.z_grid.size)

        # Get the equivalent pressure head and effective saturation.
        psi_1d, s_eff_1d = test_model.pressure_head(theta_1d, self.z_grid)

        # Get back the theta.
        theta_new, *_ = test_model(psi_1d, self.z_grid)

        # The new theta should be the same (within some error).
        self.assertTrue(isclose(np.abs(np.mean(theta_new - theta_1d)),
                                0.0, abs_tol=0.1))
        # Check the dimensions.
        self.assertEqual(theta_1d.shape, psi_1d.shape)
        self.assertEqual(theta_1d.shape, s_eff_1d.shape)

        # Test 2D version.
        psi_2d = psi_1d.repeat(4).reshape(psi_1d.size, 4)
        theta_2d, *_ = test_model(psi_2d, self.z_grid)

        # Check the dimensions.
        self.assertEqual(psi_2d.shape, theta_2d.shape)
    # _end_def_

    def test_VrettasFung_model(self):
        """
        Test Vrettas-Fung class model.

        :return: None
        """
        print(" Testing Vrettas-Fung (child class).")

        # Create a 'test-model' with default parameters.
        test_model = VrettasFung(self.soil, self.porous, self.k_hc, self.theta.res, self.dz)

        # Create a 'hypothetical' vol. water content.
        theta_1d = np.linspace(self.theta.max, self.theta.min, self.z_grid.size)

        # Get the equivalent pressure head and effective saturation.
        psi_1d, s_eff_1d = test_model.pressure_head(theta_1d, self.z_grid)

        # Create a random vector.
        n_rnd = np.random.randn(self.z_grid.size)

        # Get back the theta.
        theta_new, *_ = test_model(psi_1d, self.z_grid, {"n_rnd": n_rnd})

        # The new theta should be the same (within some error).
        self.assertTrue(isclose(np.abs(np.mean(theta_new-theta_1d)),
                                0.0, abs_tol=0.1))
        # Check the dimensions.
        self.assertEqual(theta_1d.shape, psi_1d.shape)
        self.assertEqual(theta_1d.shape, s_eff_1d.shape)

        # Test 2D version.
        psi_2d = psi_1d.repeat(4).reshape(psi_1d.size, 4)
        theta_2d, *_ = test_model(psi_2d, self.z_grid, {"n_rnd": n_rnd})

        # Check the dimensions.
        self.assertEqual(psi_2d.shape, theta_2d.shape)
    # _end_def_


if __name__ == '__main__':
    unittest.main()
