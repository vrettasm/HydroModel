import numpy as np

class SoilProperties(object):
    """
        This class represents Soil properties.
    """

    __slots__ = ("param_n", "param_a", "param_psi_sat", "param_epsilon")

    # Default constructor.
    def __init__(self, n: float = 2.0, alpha: float = 0.009, psi_sat: float = -100.0,
                 epsilon: float = 1.0e-7):
        """
        Soil properties default constructor.

        :param n: (float) is a measure of the pore-size distribution [-]
        :param alpha: (float) is related to the inverse of the air entry suction [1/L : 1/cm]
        :param psi_sat: (float) suction at saturation [L: cm]
        :param epsilon: (float) tiny positive threshold value [-]
        """

        # Check for allowed range.
        if n > 1.0:
            # Pore size distribution parameter.
            self.param_n = n
        else:
            raise ValueError(" Soil property 'n' should be > 1.")
        # _end_if_

        # Check for allowed range.
        if alpha > 0.0:
            self.param_a = alpha
        else:
            raise ValueError(" Soil property 'alpha' should be strictly positive.")
        # _end_if_

        # Pressure head (suction) should only be negative.
        # If zero value is given it means full saturation.
        self.param_psi_sat = np.minimum(psi_sat, 0.0)

        # This should be a small positive number.
        self.param_epsilon = np.maximum(epsilon, 1.0e-8)
    # _end_def_

    @property
    def n(self):
        return self.param_n
    # _end_def_

    @n.setter
    def n(self, new_value):
        # Accept only positive values.
        if new_value > 1.0:
            # Make the change.
            self.param_n = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'n': {},"
                             " should be > 1.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def m(self):
        return 1.0 - (1.0/self.param_n)
    # _end_def_

    @property
    def alpha(self):
        return self.param_a
    # _end_def_

    @alpha.setter
    def alpha(self, new_value):
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_a = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'alpha': {},"
                             " should be strictly positive.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def psi_sat(self):
        return self.param_psi_sat
    # _end_def_

    @psi_sat.setter
    def psi_sat(self, new_value):
        # Accept only negative values.
        if new_value <= 0.0:
            # Make the change.
            self.param_psi_sat = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'psi_sat': {},"
                             " should be <= 0.".format(new_value))
        # _end_if_
    # _end_def_

    @property
    def epsilon(self):
        return self.param_epsilon
    # _end_def_

    @epsilon.setter
    def epsilon(self, new_value):
        # Accept only negative values.
        if new_value > 0.0:
            # Make the change.
            self.param_epsilon = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" Soil property 'epsilon': {},"
                             " should be (tiny) positive.".format(new_value))
        # _end_if_
    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the SoilProperties object.
        This will include its id(), along with its fields values.

        :return: a string representation of a SoilProperties object.
        """
        return " SoilProperties Id({0}): n={1}, alpha={2}, psi_sat={3}," \
               " epsilon={4}".format(id(self), self.param_n, self.param_a,
                                     self.param_psi_sat, self.param_epsilon)
    # _end_def_

# _end_class_
