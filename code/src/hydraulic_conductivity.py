from numpy import maximum as np_max

class HydraulicConductivity(object):
    """
    This class represents the hydraulic conductivity of the hydrologic model.
    """

    __slots__ = ("param_sat_soil", "param_sat_saprolite", "param_sat_fresh_bedrock",
                 "param_sigma_noise", "param_lambda_exponent")

    # Default constructor.
    def __init__(self, sat_soil: float = 8.5, sat_saprolite: float = 3.2, sat_fresh_bedrock: float = 0.1,
                 sigma_noise: float = 2.0, lambda_exponent: float = 1.0):
        """
        Constructs an object that will hold all the hydraulic conductivity variables.

        Note: These values are defined at the top of the layer and can change by several
        orders of magnitude between the different soil layers.

        :param sat_soil: (float) saturation at the soil layer ([L/T: cm/0.5hrs]).

        :param sat_saprolite: (float) saturation at the saprolite layer ([L/T: cm/0.5hrs]).

        :param sat_fresh_bedrock: (float) saturation at the fresh bedrock layer ([L/T: cm/0.5hrs]).

        :param sigma_noise: (float) noise amplitude (dimensionless [-]).

        :param lambda_exponent: (float) power parameter (conceptual).
        """

        # Soil related variables.
        if sat_soil > 0.0:
            self.param_sat_soil = sat_soil
        else:
            raise ValueError(" {0}: The saturated value of the Soil layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__,
                                                                         sat_soil))
        # _end_if_

        if sat_saprolite > 0.0:
            self.param_sat_saprolite = sat_saprolite
        else:
            raise ValueError(" {0}: The saturated value of the Saprolite layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__,
                                                                         sat_saprolite))
        # _end_if_

        if sat_fresh_bedrock > 0.0:
            self.param_sat_fresh_bedrock = sat_fresh_bedrock
        else:
            raise ValueError(" {0}: The saturated value of the Fresh Bedrock layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__,
                                                                         sat_fresh_bedrock))
        # _end_if_

        # Noise related variables.
        self.param_sigma_noise = np_max(sigma_noise, 0.0)
        self.param_lambda_exponent = np_max(lambda_exponent, 0.0)
    # _end_def_

    @property
    def sat_soil(self):
        return self.param_sat_soil
    # _end_def_

    @sat_soil.setter
    def sat_soil(self, new_value):
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_soil = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" {0}: The saturated value of the Soil layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__, new_value))
        # _end_if_
    # _end_def_

    @property
    def sat_saprolite(self):
        return self.param_sat_saprolite
    # _end_def_

    @sat_saprolite.setter
    def sat_saprolite(self, new_value):
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_saprolite = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" {0}: The saturated value of the Saprolite layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__, new_value))
        # _end_if_
    # _end_def_

    @property
    def sat_fresh_bedrock(self):
        return self.param_sat_fresh_bedrock
    # _end_def_

    @sat_fresh_bedrock.setter
    def sat_fresh_bedrock(self, new_value):
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_fresh_bedrock = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" {0}: The saturated value of the Fresh Bedrock layer:"
                             " {1}, should be strictly positive.".format(self.__class__.__name__, new_value))
        # _end_if_
    # _end_def_

    @property
    def sigma_noise(self):
        return self.param_sigma_noise
    # _end_def_

    @sigma_noise.setter
    def sigma_noise(self, new_value):
        # Accept only positive values.
        if new_value >= 0.0:
            # Make the change.
            self.param_sigma_noise = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" {0}: The sigma amplitude value of the noise model:"
                             " {1}, should be non-negative.".format(self.__class__.__name__, new_value))
        # _end_if_
    # _end_def_

    @property
    def lambda_exponent(self):
        return self.param_lambda_exponent
    # _end_def_

    @lambda_exponent.setter
    def lambda_exponent(self, new_value):
        # Accept only positive values.
        if new_value >= 0.0:
            # Make the change.
            self.param_lambda_exponent = new_value
        else:
            # Raise an error with a message.
            raise ValueError(" {0}: The lambda exponent value of the noise model:"
                             " {1}, should be non-negative.".format(self.__class__.__name__, new_value))
        # _end_if_
    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the HydraulicConductivity object.
        This will include its id(), along with its fields values.

        :return: a string representation of a HydraulicConductivity object.
        """
        return " HydraulicConductivity Id({0}):"\
               " Sat-Soil={1}, Sat-Saprolite={2}, Sat-Fresh-Bedrock={3},"\
               " Sigma={4}, Lambda={5}".format(id(self), self.sat_soil, self.sat_saprolite,
                                               self.sat_fresh_bedrock, self.sigma_noise, self.lambda_exponent)
    # _end_def_

# _end_class_
