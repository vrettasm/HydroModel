import numpy as np
from scipy.interpolate import interp1d

from .hydrological_model import HydrologicalModel
from ..utilities import logN_rnd


class VrettasFung(HydrologicalModel):
    """
        This class represents a hydrologic model based on the Vrettas-Fung model.

        -- REFERENCE(S) --
        Vrettas, M. D., and I. Y. Fung (2015), "Toward a new parameterization of hydraulic conductivity
        in climate models: Simulation of rapid groundwater fluctuations in Northern California", Journal
        of Advances in Modeling Earth Systems, 07, doi:10.1002/2015MS000516.

        Vrettas, M. D., and I. Y. Fung (2017), "Sensitivity of transpiration to subsurface properties:
        Exploration with a 1-D model", Journal of Advances in Modelling Earth Systems, 9,
        doi:10.1002/2016MS000901.
        -- REFERENCE(S) --

    """

    def __init__(self, soil, porous, k_hc, theta_res, dz):
        # Call the constructor of the parent class.
        super().__init__(soil, porous, k_hc, theta_res, dz)
    # _end_def_

    def __call__(self, psi, z, *args):
        """
        A direct call to an object of this class will return the water content, along other
        related quantities, at a specific depth 'z', given the input pressure head (suction).

        :param psi: pressure head (suction) [dim_m x dim_d].

        :param z: depth values (increasing downwards) [dim_d x 1].

        :param args: in here we pass additional parameters for the noise model.

        :return: q (water content), K (unsaturated hydraulic conductivity),
        C (specific moisture capacity), Kbkg (background hydraulic conductivity),
        and q_inf_max (max infiltration capacity) [dim_d x dim_m].

        :raises ValueError: if there is a mismatch in the input dimensions.
        """

        # Ensure the input is 1-D.
        z, psi = np.atleast_1d(z, psi)

        # Get the vector size.
        dim_d, dim_m = psi.shape[0], None

        # Check if the input is 2D.
        if len(psi.shape) == 2:
            dim_m = psi.shape[1]
        # _end_if_

        # Check the input dimensions (of the vertical domain).
        if dim_d != z.shape[0]:
            raise ValueError(" {0}: Input size dimensions do not match:"
                             " {1} not equal to {2}.".format(self.__class__.__name__, dim_d, z.shape[0]))
        # _end_if_

        # Get the porosity field at 'z'.
        porous_z, *_ = self.porous(z)

        # Make sure the porosity is at least 1-D.
        porous_z = np.atleast_1d(porous_z)

        # Vectorized version.
        if dim_m:
            porous_z = porous_z.repeat(dim_m).reshape(dim_d, dim_m)
        # _end_if_

        # Initialise at None. This will cause an error
        # if the n_rnd is not given as input. (revisit)
        n_rnd = None

        # Extract additional parameters.
        if "n_rnd" in args[0]:
            n_rnd = np.atleast_1d(args[0]["n_rnd"])
        # _end_if_

        # Pre-compute constant parameters.
        delta_s = porous_z - self.theta_res

        # Check if there are saturated cells.
        id_sat = np.where(psi >= self.psi_sat)

        # Compute the volumetric moisture content in unsaturated cells.
        q = self.theta_res + delta_s * (1.0 + (self.alpha * np.abs(psi)) ** self.n) ** (-self.m)

        # Volumetric water content in saturated cells.
        q[id_sat] = porous_z[id_sat]

        # Compute the effective saturation (Se \in [0,1]).
        # (i.e. the "normalized" water content)
        s_eff = (q - self.theta_res) / delta_s

        # SAFEGUARD: THIS SHOULD NOT HAPPEN.
        s_eff = np.minimum(np.maximum(s_eff, 0.0), 1.0)

        # Get all the underground boundaries.
        (l0, l1, l2, l3) = self.porous.layers

        # Find the indexes of each underground layer.
        # NB: To ensure continuity in the interpolating values
        # we keep the equality symbol (>=) on all three layers.
        soil_layer_idx = np.where((z >= l0) & (z <= l1))
        sapr_layer_idx = np.where((z >= l1) & (z <= l2))
        wbed_layer_idx = np.where((z >= l2) & (z <= l3))

        # Create an 1_[d x m] array.
        ones_dm = np.ones(psi.shape)

        # Initialize the Unsaturated Hydraulic Conductivity.
        K = self.k_hc.sat_soil * (s_eff ** self.k_hc.lambda_exponent)

        # Initialize the Background Hydraulic Conductivity.
        k_bkg = self.k_hc.sat_soil * ones_dm

        # SOIL LAYER
        if soil_layer_idx:
            # Number of cells in the soil layer.
            n_soil = soil_layer_idx[0].size

            # Set the mean values of $Kbkg_soil$.
            mean_soil = self.k_hc.sat_soil * np.ones(n_soil)

            # Soil noise variables: !!! REDUCED !!!
            rnd_soil = 0.05 * n_rnd[soil_layer_idx]

            # Weight function:
            s_eff_soil = s_eff[soil_layer_idx]

            # Noise variance.
            sigma_soil = self.k_hc.sigma_noise * (1.0 - s_eff_soil)

            # Vectorized version.
            if dim_m:
                # Update the value of $Kbkg_soil$.
                k_bkg[soil_layer_idx, :] = logN_rnd(mean_soil.repeat(dim_m).reshape(n_soil, dim_m),
                                                    sigma_soil,
                                                    rnd_soil.repeat(dim_m).reshape(n_soil, dim_m))

                # Hydraulic conductivity for Soil layer.
                K[soil_layer_idx, :] = (s_eff_soil ** self.k_hc.lambda_exponent) * k_bkg[soil_layer_idx, :]
            else:
                # Update the value of $Kbkg_soil$.
                k_bkg[soil_layer_idx] = logN_rnd(mean_soil, sigma_soil, rnd_soil)

                # Hydraulic conductivity for Soil layer.
                K[soil_layer_idx] = (s_eff_soil ** self.k_hc.lambda_exponent) * k_bkg[soil_layer_idx]
            # _end_if_
        # _end_if_

        # SAPROLITE LAYER:
        if sapr_layer_idx:
            # Number of cells in the soil layer.
            n_sapr = sapr_layer_idx[0].size

            # Test domain for the interpolation function.
            z_test = np.arange(l1, l2+1)

            # Make a interpolation function.
            fun_sapr = interp1d(z_test, np.linspace(self.k_hc.sat_soil,
                                                    self.k_hc.sat_saprolite,
                                                    z_test.size))
            # Compute the mean values of $Kbkg_sap$.
            mean_sapr = fun_sapr(z[sapr_layer_idx])

            # Saprolite noise variables: !!! REDUCED !!!
            rnd_sapr = 0.10 * n_rnd[sapr_layer_idx]

            # Weight function:
            s_eff_sapr = s_eff[sapr_layer_idx]

            # Noise variance.
            sigma_sapr = self.k_hc.sigma_noise * (1.0 - s_eff_sapr)

            # Vectorized version.
            if dim_m:
                # Update the value of $Kbkg_sapr$.
                k_bkg[sapr_layer_idx, :] = logN_rnd(mean_sapr.repeat(dim_m).reshape(n_sapr, dim_m),
                                                    sigma_sapr,
                                                    rnd_sapr.repeat(dim_m).reshape(n_sapr, dim_m))

                # Hydraulic conductivity for Saprolite layer.
                K[sapr_layer_idx, :] = (s_eff_sapr ** self.k_hc.lambda_exponent) * k_bkg[sapr_layer_idx, :]
            else:
                # Update the value of $Kbkg_sapr$.
                k_bkg[sapr_layer_idx] = logN_rnd(mean_sapr, sigma_sapr, rnd_sapr)

                # Hydraulic conductivity for Saprolite layer.
                K[sapr_layer_idx] = (s_eff_sapr ** self.k_hc.lambda_exponent) * k_bkg[sapr_layer_idx]
            # _end_if_
        # _end_if_

        # WEATHERED BEDROCK LAYER:
        if wbed_layer_idx:
            # Number of cells in the soil layer.
            n_wbed = wbed_layer_idx[0].size

            # Compute the mean values of $Kbkg_wbed$.
            z_test = np.arange(l2, l3+1)

            # Compute the two parameters of the exponential function:
            p0 = self.k_hc.sat_saprolite
            p1 = np.log(p0 / self.k_hc.sat_fresh_bedrock) / l3

            # Construct the interpolation function.
            fun_wbed = interp1d(z_test, p0 * np.exp(-np.linspace(0, l3, z_test.size) * p1))

            # Get the values at the 'z' (weathered bedrock only).
            mean_wbed = fun_wbed(z[wbed_layer_idx])

            # Weathered bedrock noise variables:
            rnd_wbed = n_rnd[wbed_layer_idx]

            # Weight function:
            s_eff_wbed = s_eff[wbed_layer_idx]

            # Noise variance.
            sigma_wbed = self.k_hc.sigma_noise * (1.0 - s_eff_wbed)

            # Vectorized version.
            if dim_m:
                # Update the value of $Kbkg_wbed$.
                k_bkg[wbed_layer_idx, :] = logN_rnd(mean_wbed.repeat(dim_m).reshape(n_wbed, dim_m),
                                                    sigma_wbed,
                                                    rnd_wbed.repeat(dim_m).reshape(n_wbed, dim_m))

                # Hydraulic conductivity for Saprolite layer.
                K[wbed_layer_idx, :] = (s_eff_wbed ** self.k_hc.lambda_exponent) * k_bkg[wbed_layer_idx, :]
            else:
                # Update the value of $Kbkg_wbed$.
                k_bkg[wbed_layer_idx] = logN_rnd(mean_wbed, sigma_wbed, rnd_wbed)

                # Hydraulic conductivity for Saprolite layer.
                K[wbed_layer_idx] = (s_eff_wbed ** self.k_hc.lambda_exponent) * k_bkg[wbed_layer_idx]
            # _end_if_
        # _end_if_

        # SAFEGUARD:
        K[id_sat] = k_bkg[id_sat]
        K = np.minimum(K, k_bkg)

        # Compute the Specific Moisture Capacity [dTheta/dPsi].
        C = (self.m * self.n) * self.alpha * \
            delta_s * (s_eff ** (1.0 / self.m + 1.0)) * (self.alpha * np.abs(psi)) ** (self.n-1.0)

        # Set explicitly the saturated cells to the minimum value.
        C[id_sat] = self.epsilon

        # Bad conditions for the Specific moisture capacity.
        bad_condition = (C < self.epsilon) | ~np.isfinite(C) | (np.imag(C) != 0.0)

        # Replace with a small positive number to improve numerical stability.
        C[bad_condition] = self.epsilon

        # Space discretization [L: cm].
        dz = self.dz

        # After: (Collins and Bras, 2007).
        # Here we assume that the equation is solved for 30min
        # time intervals, hence: dt = 0.5 and dz/dz --> 2.0*dz
        q_inf_max = np.minimum(2.0*(porous_z[0] - q[0])*dz, k_bkg[0])

        # Tuple with all the related variables.
        return np.atleast_1d(q, K, C, k_bkg, q_inf_max)
    # _end_def_

# _end_class_
