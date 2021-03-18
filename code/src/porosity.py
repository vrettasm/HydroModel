import numpy as np
from scipy.interpolate import interp1d

from .utilities import logN_rnd


class Porosity(object):
    """
    This class creates a porosity object. A porosity profile (field) is generated,
    as a function of the depth (z), along with the relevant field capacity and
    wilting point profiles (as functions of depth).

    Naturally the porosity in the soil decreases with depth because the free empty
    pores are smaller in size, therefore allowing for less water to be kept at any
    given depth. However, depending on the type (or mixture) of soils this profile
    can vary.
    """

    __slots__ = ("profile", "field_cap", "wilting_point", "f_interp", "p_model",
                 "p_layers")

    def __init__(self, z_grid, layers, theta, soil, p_model):
        """
        Constructs an object that will hold the porosity, field capacity and
        wilting point profiles, for the z-domain (increasing downwards).

        :param z_grid: depth values at which we want to generate the porosity
        values.

        :param layers: tuple with the underground layers.

        :param theta: water content object.

        :param soil: soil properties object.

        :param p_model: porosity profile type:
                        1) Constant (for modeling short depths)
                        2) Linear (typical choice for small depths)
                        3) Exponential (typical choice for medium depths)
                        4) Stratified (for modeling high depths)
                        5) Noisy (assumes random "fractures" in the underground)

        :raises ValueError: if there is no spatial input or wrong porosity profile.

        :raises RuntimeError: if the domain is not increasing.
        """
        # Ensure the input is 1-D.
        z_grid = np.atleast_1d(z_grid)

        # Make sure the z-grid domain is not empty.
        if not z_grid.size:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" Input array z_grid is empty.")
        # _end_if_

        # Make sure the z-grid domain is increasing.
        if np.any(np.diff(z_grid) <= 0.0):
            raise RuntimeError(f" {self.__class__.__name__}:"
                               f" Space domain z_grid is not increasing.")
        # _end_if_

        # Store the type of profile.
        self.p_model = p_model

        # Get the size of the spatial grid.
        len_z = z_grid.size

        # Extract the limits of the soil layers:
        # (1) Soil, (2) Saprolite, (3) Weathered Bedrock, (4) Fresh Bedrock.
        (_, l_sapr, l_wbed, l_fbed) = layers

        # Add the layers to the object.
        self.p_layers = layers

        # Find the indices of each underground layer.
        sap_layer_idx = (z_grid >= l_sapr) & (z_grid <= l_wbed)
        web_layer_idx = (z_grid >= l_wbed) & (z_grid <= l_fbed)

        # Create the profile according to the selected type.
        if str.upper(p_model) == "CONSTANT":

            # The porosity is set to its maximum value.
            q_sat = theta.max * np.ones(len_z)

        elif str.upper(p_model) == "LINEAR":

            # Increase linearly from max (top of the domain)
            # to min (bottom of the domain).
            q_sat = np.linspace(theta.max, theta.min, len_z)

        elif str.upper(p_model) == "EXPONENTIAL":

            # Compute the two parameters of the exponential function:
            # f(z) = p0*exp(-z*p1),
            # where f(0) = theta_max, and f(z_max) = theta_min
            p0 = theta.max
            p1 = np.log(p0 / theta.min) / z_grid[-1]

            # Fit the exponential function.
            q_sat = p0 * np.exp(-z_grid * p1)

        elif str.upper(p_model) == "STRATIFIED":

            # Initialize the field at maximum value.
            q_sat = theta.max * np.ones(len_z)

            if np.any(sap_layer_idx):
                # Fit a straight (decreasing) line in the saprolite layer.
                q_sat[sap_layer_idx] = np.linspace(theta.max, theta.mid, sap_layer_idx.sum())
            # _end_if_

            if np.any(web_layer_idx):
                # Compute the two parameters of the exponential function:
                # f(z) = p0*exp(-z*p1),
                p0 = theta.mid
                p1 = np.log(p0 / theta.min) / z_grid[-1]

                # Fit the exponential function.
                q_sat[web_layer_idx] = p0 * np.exp(-np.linspace(0, l_fbed, web_layer_idx.sum()) * p1)
            # _end_if_

        elif str.upper(p_model) == "NOISY":
            # Create a random number generator.
            rng = np.random.default_rng()

            # Variance Scale Factor (small positive number).
            _vsf = 5.0e-5

            # Initialize the field at maximum value.
            q_sat = theta.max * np.ones(len_z)

            if np.any(sap_layer_idx):
                # Number of cells in the saprolite layer.
                n_sap = sap_layer_idx.sum()

                # Mean (saprolite) values.
                mean_sap = np.linspace(theta.max, theta.mid, n_sap)

                # Update the profile in the saprolite zone.
                q_sat[sap_layer_idx] = logN_rnd(mean_sap, _vsf * np.min(mean_sap),
                                                rng.standard_normal(n_sap))
            # _end_if_

            if np.any(web_layer_idx):
                # Number of cells in the weathered bedrock layer.
                n_web = web_layer_idx.sum()

                # Compute the two parameters of the exponential function:
                p0 = theta.mid
                p1 = np.log(p0 / theta.min) / z_grid[-1]

                # Mean values.
                mean_web = p0 * np.exp(-np.linspace(0, l_fbed, n_web) * p1)

                # Update the profile in the weathered bedrock zone.
                q_sat[web_layer_idx] = logN_rnd(mean_web, _vsf * np.min(mean_web),
                                                rng.standard_normal(n_web))
            # _end_if_

        else:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" Wrong porosity profile type: {p_model}")
        # _end_if_

        # Safeguard: make sure this profile is [MIN <= q_sat <= MAX]
        q_sat = np.minimum(np.maximum(q_sat, theta.min), theta.max)

        # Make sure the profile is flatten before assigning it to the object.
        self.profile = q_sat.flatten()

        # Helper function.
        def fun_wrc(psi):
            return theta.res + (q_sat - theta.res) /\
                   (1.0 + (soil.alpha * psi) ** soil.n) ** soil.m
        # _end_def_

        # Field capacity over the whole z-domain.
        self.field_cap = np.maximum(fun_wrc(theta.flc), theta.res)

        # Wilting points over the whole z-domain.
        self.wilting_point = np.minimum(fun_wrc(theta.wlt), self.field_cap)

        # Interpolate the profiles on the grid and store the functions.
        self.f_interp = interp1d(z_grid, [self.profile, self.field_cap, self.wilting_point])
    # _end_def_

    def __call__(self, z_new=None):
        """
        A call of the porosity object will return a tuple containing all three
        profiles (i.e. porosity, field capacity, wilting point) for the specific
        input 'z_new'. If no new depth values 'z' are given, the call will return
        the full profiles.

        :param z_new: the depth(s) at which we want the new profiles.

        :return: the porosity profile at depth(s) 'z'.
        """

        # Check if we have given input.
        if np.any(z_new):
            # Return all three profiles at the new input 'z_new'.
            return self.f_interp(z_new)
        else:
            # Return all the profiles (over the whole z-domain).
            return self.profile, self.field_cap, self.wilting_point
        # _end_if_

    # _end_def_

    @property
    def layers(self):
        """
        Accessor of the underground layers.

        :return: Tuple with the (L0, L1, ...)
        """
        return self.p_layers
    # _end_def_

    def __str__(self):
        """
        Override to print a readable string presentation of the Porosity object.
        At the moment this will include its id(), along with its profile type.

        :return: a string representation of a Porosity object.
        """
        return f" Porosity Id({id(self)}): Type={self.p_model}"
    # _end_def_

# _end_class_
