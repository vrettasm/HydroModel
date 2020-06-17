import numpy as np
from scipy.interpolate import interp1d
from code.src.utilities import logN_rnd
from code.src.models.hydrological_model import HydrologicalModel

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

        :param psi: pressure head (suction) [dim_d x dim_m]

        :param z: depth values (increasing downwards) [dim_d x 1]

        :param args: in here we pass additional parameters for the noise model.

        :return: q (water content), K (unsaturated hydraulic conductivity),
        C (specific moisture capacity) and q_inf_max (max infiltration capacity).
        """

        # Make sure the input 'psi' has at least shape (d, 1).
        if len(psi.shape) == 1:
            psi = np.reshape(psi, (psi.size, 1))
        # _end_if_

        # Make sure the input 'z' has shape (d, 1).
        if len(z.shape) == 1:
            z = np.reshape(z, (z.size, 1))
        # _end_if_

        # Get the dimensions of the input array.
        dim_d, dim_m = psi.shape

        # Check the input dimensions (of the vertical domain).
        if dim_d != z.size:
            raise ValueError(" Input size dimensions do not match:"
                             " {0} not equal to {1}".format(dim_d, z.size))
        # _end_if_

        # Get the porosity at 'z'.
        porous_z = self.porous(z)

        # Initialize 'q' (volumetric water content) variable.
        q = np.zeros(dim_d)

        # Initialise at None.
        n_rnd = None

        # Extract additional parameters.
        for arg_i in args:
            # Check if the argument is a dictionary.
            if isinstance(arg_i, dict):
                # Extract the random noise vector.
                if "n_rnd" in arg_i:
                    n_rnd = arg_i["n_rnd"]
                # _end_if_
            # _end_if_
        # _end_for_

        # Repeat if necessary (for vectorization).
        if dim_m > 1:
            # These will be (d, m) arrays.
            q = np.repeat(q, dim_m, 1)
            n_rnd = np.repeat(n_rnd, dim_m, 1)
            porous_z = np.repeat(porous_z, dim_m, 1)
        # _end_if_

        # Pre-compute constant parameters.
        delta_s = porous_z - self.theta_res

        # Check if there are saturated cells.
        id_sat = psi >= self.psi_sat

        # Loop for each input vector.
        for i in range(dim_m):
            # Check if there are saturated cells.
            j = id_sat[:, i]

            # Volumetric water content in saturated cells.
            q[j, i] = porous_z[j, i]

            # Compute the volumetric moisture content in unsaturated cells.
            q[~j, i] = self.theta_res + \
                       delta_s[~j, i] * (1.0 + (self.alpha * np.abs(psi[~j, i])) ** self.n) ** (-self.m)
        # _end_for_

        # Compute the effective saturation (Se \in [0,1]).
        # (i.e. the "normalized" water content)
        s_eff = (q - self.theta_res) / delta_s

        # SAFEGUARD: THIS SHOULD NOT HAPPEN.
        s_eff = np.minimum(np.maximum(s_eff, 0.0), 1.0)

        # Get all the underground boundaries.
        (l0, l1, l2, l3) = porous_z.layers

        # Find the indexes of each underground layer.
        # NB: To ensure continuity in the interpolating values
        # we keep the equality symbol (>=) on all three layers.
        soil_layer_idx = np.array((z >= l0) & (z <= l1), dtype=bool)
        sapr_layer_idx = np.array((z >= l1) & (z <= l2), dtype=bool)
        wbed_layer_idx = np.array((z >= l2) & (z <= l3), dtype=bool)

        # Create an 1_[d x m] array.
        ones_dm = np.ones((dim_d, dim_m))

        # Initialize the Unsaturated Hydraulic Conductivity.
        K = self.k_hc.sat_soil * (s_eff ** (self.k_hc.lambda_exponent * ones_dm))

        # Initialize the Background Hydraulic Conductivity.
        Kbkg = self.k_hc.sat_soil * ones_dm

        # SOIL LAYER:
        if np.any(soil_layer_idx):
            # Number of cells in the soil layer.
            n_soil = soil_layer_idx.sum()

            # Set the mean values of $Kbkg_soil$.
            mean_soil = self.k_hc.sat_soil * np.ones((n_soil, dim_m))

            # Soil noise variables: !!! REDUCED !!!
            rnd_soil = 0.05 * n_rnd[soil_layer_idx]

            # Weight function:
            s_eff_soil = s_eff[soil_layer_idx]

            # Replicate 'mean' and 'sigma' parameters.
            sigma_soil = self.k_hc.sigma_noise * (1.0 - s_eff_soil)

            # Update the value of $Kbkg_soil$.
            Kbkg[soil_layer_idx] = logN_rnd(mean_soil, sigma_soil, rnd_soil)

            # Hydraulic conductivity for Soil layer.
            K[soil_layer_idx] = (s_eff_soil ** self.k_hc.lambda_exponent) * Kbkg[soil_layer_idx]
        # _end_if_

        # SAPROLITE LAYER:
        if np.any(sapr_layer_idx):
            # Test domain for the interpolation function.
            z_test = np.arange(l1, l2)

            # Make a interpolation function.
            fun_sapr = interp1d(z_test,
                                np.linspace(self.k_hc.sat_soil,
                                            self.k_hc.sat_saprolite,
                                            z_test.size))

            # Compute the mean values of $Kbkg_sap$.
            mean_sapr = fun_sapr(z[sapr_layer_idx])

            # Saprolite noise variables: !!! REDUCED !!!
            rnd_sapr = 0.10 * n_rnd[sapr_layer_idx]

            # Weight function:
            s_eff_sapr = s_eff[sapr_layer_idx]

            # Make sure mean_sapr has at least shape (d, 1).
            if len(mean_sapr.shape) == 1:
                mean_sapr = np.reshape(mean_sapr, (mean_sapr.size, 1))
            # _end_if_

            # Replicate 'mean' and 'sigma' parameters
            # for the vectorized version of the code.
            mean_sapr = np.repeat(mean_sapr, dim_m, 1)
            sigma_sapr = self.k_hc.sigma_noise * (1.0 - s_eff_sapr)

            # Update the value of $Kbkg_sapr$.
            Kbkg[sapr_layer_idx] = logN_rnd(mean_sapr, sigma_sapr, rnd_sapr)

            # Hydraulic conductivity for Saprolite layer.
            K[sapr_layer_idx] = (s_eff_sapr ** self.k_hc.lambda_exponent) * Kbkg[sapr_layer_idx]
        # _end_if_

        # WEATHERED BEDROCK LAYER:
        if np.any(wbed_layer_idx):
            # Compute the mean values of $Kbkg_wbed$.
            z_test = np.arange(l2, l3)

            # Compute the two parameters of the exponential function:
            p0 = self.k_hc.sat_saprolite
            p1 = np.log(p0 / self.k_hc.sat_fresh_bedrock) / l3

            # Construct the interpolation function.
            fun_wbed = interp1d(z_test,
                                p0 * np.exp(-np.linspace(0, l3, z_test.size) * p1))

            # Get the values at the 'z' (weathered bedrock only).
            mean_wbed = fun_wbed(z[wbed_layer_idx])

            # Weathered bedrock noise variables:
            rnd_wbed = n_rnd[wbed_layer_idx]

            # Weight function:
            s_eff_wbed = s_eff[wbed_layer_idx]

            # Make sure mean_wbed has at least shape (d, 1).
            if len(mean_wbed.shape) == 1:
                mean_wbed = np.reshape(mean_wbed, (mean_wbed.size, 1))
            # _end_if_

            # Replicate 'mean' and 'sigma' parameters
            # for the vectorized version of the code.
            mean_wbed = np.repeat(mean_wbed, dim_m, 1)
            sigma_wbed = self.k_hc.sigma_noise * (1.0 - s_eff_wbed)

            # Update the value of $Kbkg_wbed$.
            Kbkg[wbed_layer_idx] = logN_rnd(mean_wbed, sigma_wbed, rnd_wbed)

            # Hydraulic conductivity for Saprolite layer.
            K[wbed_layer_idx] = (s_eff_wbed ** self.k_hc.lambda_exponent) * Kbkg[wbed_layer_idx]
        # _end_if_

        # SAFEGUARD:
        K[id_sat] = Kbkg[id_sat]

        # SAFEGUARD:
        K = np.minimum(K, Kbkg)

        # Compute the Specific Moisture Capacity [dTheta/dPsi].
        C = (self.m * self.n) * self.alpha *\
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
        q_inf_max = np.minimum(2.0*(porous_z[0] - q[0])*dz, Kbkg[0])

        # Tuple with all the related variable.
        return q, K, C, q_inf_max
    # _end_def_

# _end_class_
