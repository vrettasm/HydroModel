import numpy as np

def logN_rnd(mx, vx, en):
    """
    Log-Normal function. If X follows the log-normal distribution with parameters 'mu'
    and 'sigma': X ~ LogN(mu, sigma), then log(X) follows the normal distribution with
    mean 'mu' and standard deviation 'sigma': log(X) ~ N(mu, sigma).

    Example:
    >> x = logN_rnd(5.0, 1.0, np.random.randn(10000))
    >> x_mean = np.mean(x)
    >> ... 4.992025079986667 (almost 5.0)

    :param mx: desired mean (of the Normal distribution).

    :param vx: desired variance (of the Normal distribution).

    :param en: random variable N(0, 1)

    :return: a random variable drawn from LogNormal(mu, sigma)
    """

    # Avoid division by zero error.
    if np.any(mx == 0.0):
        mx[mx == 0.0] = 1.0e-7
    # _end_if_

    # Parameter 'mu' for the LogNormal distribution.
    mu0 = np.log((mx ** 2) / np.sqrt(vx + mx ** 2))

    # Parameter 'sigma' for the LogNormal distribution.
    sig = np.sqrt(np.log(vx / (mx ** 2) + 1.0))

    # Random variable with (mu, sigma).
    return np.exp(mu0 + sig * en)
# _end_def_

def find_wtd(y=None):
    """
    This function identifies  the location of  the water table  depth (wtd)
    and returns  its index value.  It starts the search  from the bottom of
    the input vector and continues going upwards until the cell in question
    is not saturated. It assumes that the "wtd" exists right below the first
    non-saturated cell.

    :param y: bool vector with the saturated cells. [dim_d x 1].

    :return: the index at the first "True", from the end of the vector.
    """

    # If no input, return None.
    if y is None:
        return None
    # _end_if_

    # Get the length of the input vector.
    dim_d = y.size

    # Initialize return index. In the unlikely event that the whole state
    # vector is saturated, the wtd should point at the surface (1st cell).
    i = 0

    # Initialize the index at the bottom.
    #  - i.e. start the search backwards.
    for j, y_n in enumerate(reversed(y)):
        # Check if the cell is saturated.
        if not y_n:
            i = dim_d - j
            break
        # _end_if_
    # _end_for_

    # Make sure the index does not exceed the upper limit: return index "i"
    # should be in [0, dim_d-1].  This is because  we assume that the |wtd|
    # always exists inside the spatial domain.
    return np.minimum(i, dim_d-1)
# _end_def_
