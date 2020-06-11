import numpy as np

class SoilProperties(object):
    """
        This class represents Soil properties.
    """

    __slots__ = ("n", "a0", "psi_sat", "epsilon")

    # Default constructor.
    def __init__(self, n: float = 2.0, a0: float = 0.009, psi_sat: float = -100.0,
                 epsilon: float = 1.0e-7):
        """
        Soil properties default constructor.

        :param n: (float) is a measure of the pore-size distribution [-]
        :param a0: (float) is related to the inverse of the air entry suction [1/L : 1/cm]
        :param psi_sat: (float) suction at saturation [L: cm]
        :param epsilon: (float) tiny positive threshold value [-]
        """

        # Check for allowed range.
        if n > 1.0:
            # Pore size distribution parameter.
            self.n = n
        else:
            raise ValueError(" Soil property 'n0' should be > 1.")
        # _end_if_

        # Check for allowed range.
        if a0 > 0.0:
            self.a0 = a0
        else:
            raise ValueError(" Soil property 'a0' should be strictly positive.")
        # _end_if_

        # Pressure head (suction) should be only negative.
        # If zero value is given it means full saturation.
        self.psi_sat = np.minimum(psi_sat, 0.0)

        # This should be a small positive number.
        self.epsilon = np.maximum(epsilon, 1.0e-8)
    # _end_def_

    @property
    def param_n(self):
        return self.n
    # _end_def_

    @param_n.setter
    def param_n(self, new_value):
        # Accept only positive values.
        if new_value > 1.0:
            # Make the change.
            self.n = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'n': {},"
                             " should be > 1.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def param_m(self):
        return 1.0 - (1.0/self.n)
    # _end_def_

    @property
    def param_a0(self):
        return self.a0
    # _end_def_

    @param_a0.setter
    def param_a0(self, new_value):
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.a0 = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'a0': {},"
                             " should be strictly positive.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def psiSat(self):
        return self.psi_sat
    # _end_def_

    @psiSat.setter
    def psiSat(self, new_value):
        # Accept only negative values.
        if new_value <= 0.0:
            # Make the change.
            self.psi_sat = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'psi_sat': {},"
                             " should be <= 0.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def param_epsilon(self):
        return self.epsilon
    # _end_def_

    @param_epsilon.setter
    def param_epsilon(self, new_value):
        # Accept only negative values.
        if new_value > 0.0:
            # Make the change.
            self.epsilon = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'epsilon': {},"
                             " should be (tiny) positive.".format(new_value))
        # _end_if_
    # _end_def_

# _end_class_
