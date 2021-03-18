from numpy import maximum as np_max

class HydraulicConductivity(object):
    """
    This class represents the hydraulic conductivity of the hydrological model.
    """

    __slots__ = ("param_sat_soil", "param_sat_saprolite", "param_sigma_noise",
                 "param_sat_fresh_bedrock", "param_lambda_exponent")

    # Default constructor.
    def __init__(self, sat_soil: float = 8.5, sat_saprolite: float = 3.2,
                 sat_fresh_bedrock: float = 0.1, sigma_noise: float = 2.0,
                 lambda_exponent: float = 1.0):
        """
        Constructs an object that will hold all the hydraulic conductivity
        variables. These values are defined at the top of the layer and can
        change by several orders of magnitude between the different soil layers.

        Note: L is cm, T is 0.5hrs.

        :param sat_soil: (float) saturation at the soil layer [L/T].

        :param sat_saprolite: (float) saturation at the saprolite layer [L/T].

        :param sat_fresh_bedrock: (float) saturation at the fresh bedrock layer [L/T].

        :param sigma_noise: (float) noise amplitude (dimensionless [-]).

        :param lambda_exponent: (float) power parameter (conceptual).

        :raises ValueError: if there the input values are not in range.
        """

        # Soil related variables.
        if sat_soil > 0.0:
            self.param_sat_soil = sat_soil
        else:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Soil layer: "
                             f" {sat_soil} should be strictly positive.")
        # _end_if_

        if sat_saprolite > 0.0:
            self.param_sat_saprolite = sat_saprolite
        else:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Saprolite layer: "
                             f" {sat_saprolite} should be strictly positive.")
        # _end_if_

        if sat_fresh_bedrock > 0.0:
            self.param_sat_fresh_bedrock = sat_fresh_bedrock
        else:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Fresh Bedrock layer: "
                             f" {sat_fresh_bedrock} should be strictly positive.")
        # _end_if_

        # Noise related variables.
        self.param_sigma_noise = np_max(sigma_noise, 0.0)
        self.param_lambda_exponent = np_max(lambda_exponent, 0.0)
    # _end_def_

    @property
    def sat_soil(self):
        """
        Accessor method.

        :return: the saturated value of the soil layer.
        """
        return self.param_sat_soil
    # _end_def_

    @sat_soil.setter
    def sat_soil(self, new_value):
        """
        Accessor method.

        :param new_value: for the saturation of the soil layer.

        :return: None.
        """
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_soil = new_value
        else:
            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Soil layer: "
                             f" {new_value} should be strictly positive.")
        # _end_if_
    # _end_def_

    @property
    def sat_saprolite(self):
        """
        Accessor method.

        :return: the saturated value of the saprolite layer.
        """
        return self.param_sat_saprolite
    # _end_def_

    @sat_saprolite.setter
    def sat_saprolite(self, new_value):
        """
        Accessor method.

        :param new_value: for the saturation of the saprolite layer.

        :return: None.
        """
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_saprolite = new_value
        else:
            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Saprolite layer:"
                             f" {new_value} should be strictly positive.")
        # _end_if_
    # _end_def_

    @property
    def sat_fresh_bedrock(self):
        """
        Accessor method.

        :return: the saturated value of the fresh bedrock layer.
        """
        return self.param_sat_fresh_bedrock
    # _end_def_

    @sat_fresh_bedrock.setter
    def sat_fresh_bedrock(self, new_value):
        """
        Accessor method.

        :param new_value: for the saturation of the fresh bedrock layer.

        :return: None.
        """
        # Accept only positive values.
        if new_value > 0.0:
            # Make the change.
            self.param_sat_fresh_bedrock = new_value
        else:
            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The saturated value of the Fresh Bedrock layer:"
                             f" {new_value} should be strictly positive.")
        # _end_if_
    # _end_def_

    @property
    def sigma_noise(self):
        """
        Accessor method.

        :return: the sigma noise amplitude of the noise model.
        """
        return self.param_sigma_noise
    # _end_def_

    @sigma_noise.setter
    def sigma_noise(self, new_value):
        """
        Accessor method.

        :param new_value: of the sigma noise amplitude.

        :return: None.
        """
        # Accept only positive values.
        if new_value >= 0.0:
            # Make the change.
            self.param_sigma_noise = new_value
        else:
            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The sigma amplitude value of the noise model:"
                             f" {new_value} should be non-negative.")
        # _end_if_
    # _end_def_

    @property
    def lambda_exponent(self):
        """
        Accessor method.

        :return: the lambda exponent of the noise model.
        """
        return self.param_lambda_exponent
    # _end_def_

    @lambda_exponent.setter
    def lambda_exponent(self, new_value):
        """
        Accessor method.

        :param new_value: for the lambda exponent of the noise model.

        :return: None.
        """
        # Accept only positive values.
        if new_value >= 0.0:
            # Make the change.
            self.param_lambda_exponent = new_value
        else:
            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The lambda exponent value of the noise model:"
                             f" {new_value} should be non-negative.")
        # _end_if_
    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the object.
        This will include its id(), along with its fields values.

        :return: a string representation of a HydraulicConductivity.
        """
        # New line character.
        new_line = '\n'

        # Return the string.
        return f" HydraulicConductivity Id({id(self)}): {new_line}"\
               f" Sat-Soil={self.sat_soil} {new_line}"\
               f" Sat-Saprolite={self.sat_saprolite} {new_line}"\
               f" Sat-Fresh-Bedrock={self.sat_fresh_bedrock} {new_line}"\
               f" Sigma={self.sigma_noise} {new_line}"\
               f" Lambda={self.lambda_exponent}"
    # _end_def_

# _end_class_
