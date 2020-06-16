import numpy as np

def logN_muf(mi, vi):
    """
    Log-Normal functions.

    :param mi:

    :param vi:

    :return:
    """
    return np.log((mi ** 2) / np.sqrt(vi + mi ** 2))
# _end_def_

def logN_sig(mi, vi):
    """
    Log-Normal function

    :param mi:

    :param vi:

    :return:
    """
    return np.sqrt(np.log(vi / (mi ** 2) + 1.0))
# _end_def_

def logN_rnd(mu, sig, en):
    """
    Log-Normal function

    :param mu:

    :param sig:

    :param en:

    :return:
    """
    return np.exp(mu + sig * en)
# _end_def_
