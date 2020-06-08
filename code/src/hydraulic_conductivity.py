from numpy import minimum as np_min
from numpy import maximum as np_max

class HydraulicConductivity(object):
    """
        This class represents Water content, or soil moisture content.
    """

    __slots__ = ["sat_soil", "sat_saprolite", "sat_fresh_bedrock", "sigma_noise", "lambda_exponent"]

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
        self.sat_soil = sat_soil
        self.sat_saprolite = sat_saprolite
        self.sat_fresh_bedrock = sat_fresh_bedrock

        # Noise related variables.
        self.sigma_noise = np_max(sigma_noise, 0.0)
        self.lambda_exponent = np_max(lambda_exponent, 0.0)
    # _end_def_


# _end_class_
