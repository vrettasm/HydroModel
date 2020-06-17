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

    :param en: random variable N(0,1)

    :return: a random variable drawn from LogNormal(mu, sigma)
    """

    # Parameter 'mu' for the LogNormal distribution.
    mu0 = np.log((mx ** 2) / np.sqrt(vx + mx ** 2))

    # Parameter 'sigma' for the LogNormal distribution.
    sig = np.sqrt(np.log(vx / (mx ** 2) + 1.0))

    # Random variable with (mu, sigma).
    return np.exp(mu0 + sig * en)
# _end_def_
