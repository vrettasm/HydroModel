import numpy as np
from numba import njit

@njit(fastmath=True)
def _local_fast(mx, vx, en):
    """
    Local (hidden) function optimized with numba.njit.
    The input signature is the same as the logN_rnd().
    """
    # Square mean values.
    mx_sq = mx ** 2

    # Parameter 'mu' for the LogNormal distribution.
    mu0 = np.log(mx_sq / np.sqrt(vx + mx_sq))

    # Parameter 'sigma' for the LogNormal distribution.
    sig = np.sqrt(np.log(vx / mx_sq + 1.0))

    return np.exp(mu0 + sig * en)
# _end_def_

def logN_rnd(mx, vx, en):
    """
    Log-Normal function. If X follows the log-normal distribution with parameters 'mu'
    and 'sigma': X ~ LogN(mu, sigma), then log(X) follows the normal distribution with
    mean 'mu' and standard deviation 'sigma': log(X) ~ N(mu, sigma).

    Example (from python console):

    >> x = logN_rnd(5.0, 1.0, np.random.randn(10000))

    >> x_mean = np.mean(x)

    >> ... 4.992025079986667 (which is almost 5.0)

    :param mx: desired mean (of the Normal distribution).

    :param vx: desired variance (of the Normal distribution).

    :param en: random variable N(0, 1)

    :return: a random variable drawn from LogNormal(mu, sigma)
    """

    # Ensure the input is at least 1-D.
    mx, vx, en = np.atleast_1d(mx, vx, en)

    # Avoid division by zero error.
    mx[mx == 0.0] = 1.0e-7

    # Random variable with (mu, sigma).
    return _local_fast(mx, vx, en)
# _end_def_

@njit
def find_wtd(y=None):
    """
    This function identifies  the location of  the water table  depth (wtd)
    and returns  its index value.  It starts the search  from the bottom of
    the input vector and continues going upwards until the cell in question
    is not saturated. It assumes that the "wtd" exists right below the first
    non-saturated cell.

    :param y: bool vector with the saturated cells. (dim_d,).

    :return: the index at the first "True", from the end of the vector.
    """

    # If no input, return None.
    if y is None:
        return None
    # _end_if_

    # We want the input vector to be 1-D.
    y = np.atleast_1d(y)

    # Get the length of the input vector.
    dim_d = y.shape[0]

    # Initialize return index. In the unlikely event that the whole state
    # vector is saturated, the wtd should point at the surface (1st cell).
    i = 0

    # Initialize the index at the bottom.
    #  - i.e. start the search backwards.
    for j, y_n in enumerate(y[::-1]):
        # Check if the cell is saturated.
        if not y_n:
            i = dim_d - j
            break
        # _end_if_
    # _end_for_

    # Make sure the index does not exceed the upper limit: return index "i"
    # should be in [0, dim_d-1].  This is because  we assume that the |wtd|
    # always exists inside the spatial domain.
    return int(np.minimum(i, dim_d-1))
# _end_def_
