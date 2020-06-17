import unittest
import numpy as np
from math import isclose
from code.src.utilities import logN_rnd

class TestUtilities(unittest.TestCase):

    def test_logN_rnd(self):
        # Log-Normal parameters.
        mu, sigma = 20.5, 1.0

        # Number of samples.
        n_samples = 25000

        # Set the random seed (for reproducibility.
        np.random.seed(0)

        # Sample n_samples parameters.
        x = logN_rnd(mu, sigma, np.random.randn(n_samples))

        # Check if the result is close enough (~ 1.0e-3).
        self.assertTrue(isclose(np.mean(x), mu, rel_tol=1.0e-3))
    # _end_def_


if __name__ == '__main__':
    unittest.main()
