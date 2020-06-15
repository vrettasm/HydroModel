import numpy as np


class Porosity(object):
    """
    This class creates a porosity object. A porosity profile (field) is generated,
    as a function of the depth (z), along with the relevant field capacity and
    wilting point profiles (as functions of depth).
    """

    __slots__ = ("profile", "field_cap", "wilting_point", "sap_layer_idx", "web_Layer_idx")

    def __init__(self, z_grid, layers, theta, soil, p_model):

        # Get the size of the spatial grid.
        len_z = z_grid.size

        # Extract the limits of the soil layers:
        # (1) Soil, (2) Saprolite, (3) Weathered Bedrock, (4) Fresh Bedrock.
        (_, l_sapr, l_wbed, l_fbed) = layers

        # Find the indices of each underground layer.
        self.sap_layer_idx = (z_grid >= l_sapr) & (z_grid <= l_wbed)
        self.web_Layer_idx = (z_grid >= l_wbed) & (z_grid <= l_fbed)

        # Create the profile according to the selected type.
        if str.upper(p_model) == "CONSTANT":

            # The porosity is set to its maximum value.
            q_sat = theta.max * np.ones(len_z, 1)

        elif str.upper(p_model) == "LINEAR":

            # Increase linearly from max (top of the domain) to min (bottom of the domain).
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

            if np.any(self.sap_layer_idx):
                # Fit a straight (decreasing) line in the saprolite layer.
                q_sat[self.sap_layer_idx] = np.linspace(theta.max, theta.mid,
                                                        self.sap_layer_idx.sum())
            # _end_if_

            if np.any(self.web_Layer_idx):
                # Compute the two parameters of the exponential function:
                # f(z) = p0*exp(-z*p1),
                p0 = theta.mid
                p1 = np.log(p0 / theta.min) / z_grid[-1]

                # Fit the exponential function.
                q_sat[self.web_Layer_idx] = p0 * np.exp(-np.linspace(0, l_fbed,
                                                                     self.web_Layer_idx.sum()) * p1)
            # _end_if_

        elif str.upper(p_model) == "NOISY":

            # Log-Normal functions.
            logN_muf = lambda m, v: np.log(m ** 2) / np.sqrt(v + m ** 2)
            logN_sig = lambda m, v: np.sqrt(np.log(v / (m ** 2) + 1.0))
            logN_rnd = lambda mu, sig, en: np.exp(mu + sig * en)

            # Variance Scale Factor (small positive number).
            _vsf = 5.0e-5

            # Initialize the field at maximum value.
            q_sat = theta.max * np.ones(len_z)

            if np.any(self.sap_layer_idx):
                # Number of cells in the saprolite layer.
                n_sap = self.sap_layer_idx.sum()

                # Mean (saprolite) values.
                mean_sap = np.linspace(theta.max, theta.mid, n_sap)

                # Noise model components (mean, variance).
                mu_sap = logN_muf(mean_sap, _vsf * np.min(mean_sap))
                sig_sap = logN_sig(mean_sap, _vsf * np.min(mean_sap))

                # Update the profile in the sparolite zone.
                q_sat[self.sap_layer_idx] = logN_rnd(mu_sap, sig_sap,
                                                     np.random.randn(n_sap, 1))
            # _end_if_

            if np.any(self.web_Layer_idx):
                # Number of cells in the weathered bedrock layer.
                n_web = self.web_Layer_idx.sum()

                # Compute the two parameters of the exponential function:
                p0 = theta.mid
                p1 = np.log(p0 / theta.min) / z_grid[-1]

                # Mean values.
                mean_web = p0 * np.exp(-np.linspace(0, l_fbed, n_web) * p1)

                # Noise model components (mean, variance).
                mu_web = logN_muf(mean_web, _vsf * np.min(mean_web))
                sig_web = logN_sig(mean_web, _vsf * np.min(mean_web))

                # Update the profile in the weathered bedrock zone.
                q_sat[self.web_Layer_idx] = logN_rnd(mu_web, sig_web,
                                                     np.random.randn(n_web, 1))
            # _end_if_

        else:
            raise ValueError(" Wrong porosity profile type: {0}".format(p_model))
        # _end_if_

        # Safeguard: make sure this profile is [MIN <= q_sat <= MAX]
        q_sat = np.minimum(np.maximum(q_sat, theta.min), theta.max)

        # Make sure the profile is flatten before assigning it to the object.
        self.profile = q_sat.flatten()

        # Extract soil parameters.
        n = soil.n
        m = soil.m
        a0 = soil.alpha

        # helper function.
        fun_wrc = lambda psi: theta.res + (q_sat - theta.res) / (1.0 + (a0 * psi) ** n) ** m

        # Field capacity over the whole z-domain.
        self.field_cap = np.maximum(fun_wrc(theta.flc), theta.res)

        # Wilting points over the whole z-domain.
        self.wilting_point = np.minimum(fun_wrc(theta.wlt), self.field_cap)
    # _end_def_

    def __call__(self, z=None):
        """
        TBD

        :param z:
        :return: the porosity profile at depth(s) 'z'.
        """

        # Return all the profiles (over the whole domain).
        if not z:
            return self.profile
        # _end_if_

        # Get the pre-computed porosity profile (over the whole domain).
        porosity_z = self.profile

        return porosity_z

    # _end_def_

    def __str__(self):
        pass
    # _end_def_

# _end_class_
