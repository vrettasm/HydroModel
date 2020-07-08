import os
import sys
import unittest
import numpy as np
from math import isclose
from numpy.random import default_rng

# Make sure we can import the /src.
sys.path.append(os.path.abspath("../../code"))

from src.utilities import logN_rnd, find_wtd


class TestUtilities(unittest.TestCase):

    def test_logN_rnd(self):
        # Log-Normal parameters.
        mu, sigma = 20.5, 1.0

        # Number of samples.
        n_samples = 25000

        # Get a random number generator.
        # Set the seed = 0 (for reproducibility).
        rng = default_rng(0)

        # Sample n_samples parameters.
        x = logN_rnd(mu, sigma, rng.standard_normal(n_samples))

        # Check if the result is close enough (~ 1.0e-3).
        self.assertTrue(isclose(np.mean(x), mu, rel_tol=1.0e-3))

        # Check the vector version (1D).
        num = 100
        mu_vec = np.linspace(10, 1, num)
        sigma_ = np.ones(mu_vec.shape)

        # Sample n_samples parameters.
        x_out = logN_rnd(mu_vec, sigma_, rng.standard_normal(num))

        # Check if the shapes match.
        self.assertEqual(mu_vec.shape, x_out.shape)

        # Check the vector version (2D).
        dim_d, dim_m = 100, 4
        mu_mat = np.random.rand(dim_d, dim_m)
        sigma_ = np.ones(mu_mat.shape)

        # Sample n_samples parameters.
        x_out = logN_rnd(mu_mat, sigma_, rng.standard_normal((dim_d, dim_m)))

        # Check if the shapes match.
        self.assertEqual(mu_mat.shape, x_out.shape)
    # _end_def_

    def test_find_wtd(self):
        # Extreme case no.1
        # There are no saturated cells in the domain vector.
        x1 = np.array([False, False, False, False, False, False, False, False, False, False])

        # The expected index is at the end of the domain (i_correct = 9).
        self.assertEqual(9, find_wtd(x1))

        # Extreme case no.2
        # All the cells are saturated in the domain vector.
        x2 = np.array([True, True, True, True, True, True, True, True, True, True])

        # The expected index is at the beginning of the domain (i_correct = 0).
        self.assertEqual(0, find_wtd(x2))

        # Typical case no.3
        # The saturated cells are all at the bottom of the domain vector.
        x3 = np.array([False, False, False, False, False, False, True, True, True, True])

        # The expected index is at the first "True",
        # from the end of the domain (i_correct = 6).
        self.assertEqual(6, find_wtd(x3))

        # Typical case no.4
        # Variably saturated cells in the domain vector.
        x3 = np.array([False, True, False, True, False, False, False, True, True, True])

        # The expected index is at the first "True",
        # from the end of the domain (i_correct = 7).
        self.assertEqual(7, find_wtd(x3))
    # _end_def_


if __name__ == '__main__':
    unittest.main()
