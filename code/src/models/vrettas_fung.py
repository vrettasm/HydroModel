import numpy as np
from code.src.models.hydrological_model import HydrologicalModel

class VrettasFung(HydrologicalModel):
    """
        This class represents a hydrologic model based on the Vrettas-Fung model.

        -- REFERENCE(S) --
        Vrettas, M. D., and I. Y. Fung (2015), "Toward a new parameterization of hydraulic conductivity in climate
        models: Simulation of rapid groundwater fluctuations in Northern California", Journal of Advances in Modeling
        Earth Systems, 07, doi:10.1002/2015MS000516.

        Vrettas, M. D., and I. Y. Fung (2017), "Sensitivity of transpiration to subsurface properties: Exploration
        with a 1-D model", Journal of Advances in Modelling Earth Systems, 9, doi:10.1002/2016MS000901.
        -- REFERENCE(S) --

    """

    def __init__(self, soil, porous, k_sat, theta_res, dz):
        # Call the constructor of the parent class.
        super().__init__(soil, porous, k_sat, theta_res, dz)
    # _end_def_

    def __call__(self, psi, z):
        """

        :param psi:

        :param z:

        :return:
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

        # Create a vector with the K_{sat} values.
        k_bkg = self.k_sat * np.ones(dim_d)

        # Get the porosity at 'z'.
        porous_z = self.porous(z)

        # Initialize 'q' (water) variable.
        q = np.zeros(dim_d)

        # Repeat if necessary (for vectorization).
        if dim_m > 1:
            q = np.repeat(q, dim_m, 1)
            k_bkg = np.repeat(k_bkg, dim_m, 1)
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

        # Pre-compute the $m-th$ root.
        mth_s_eff = s_eff ** (1.0 / self.m)

        # Get *ALL* the underground boundaries.
        (l0, l1, l2, l3) = porous_z.layers

        # Find the indexes of each underground layer.
        # NB: To ensure continuity in the interpolating values
        # we keep the equality symbol (>=) on all three layers.
        soil_layer_idx = np.array((z >= l0) & (z <= l1), dtype=bool)
        sapr_layer_idx = np.array((z >= l1) & (z <= l2), dtype=bool)
        wbed_layer_idx = np.array((z >= l2) & (z <= l3), dtype=bool)

        # SAFEGUARD:
        K[id_sat] = k_bkg[id_sat]

        # SAFEGUARD:
        K = np.minimum(K, k_bkg)

        # Compute the Specific Moisture Capacity [dTheta/dPsi].
        C = (self.m * self.n) * self.alpha *\
            delta_s * (s_eff ** (1.0 / self.m + 1.0)) * (self.alpha * np.abs(psi)) ** (self.n-1.0)

        # Set explicitly the saturated cells to the minimum value.
        C[id_sat] = self.epsilon

        # bad conditions for the Specific moisture capacity.
        bad_condition = (C < self.epsilon) | ~np.isfinite(C) | (np.imag(C) != 0.0)

        # Replace with a small positive number to improve numerical stability.
        C[bad_condition] = self.epsilon

        # Space discretization [L: cm].
        dz = self.dz

        # After: (Collins and Bras, 2007).
        # Here we assume that the equation is solved for 30min
        # time intervals, hence: dt = 0.5 and dz/dz --> 2.0*dz
        q_inf_max = np.minimum(2.0*(porous_z[0] - q[0])*dz, k_bkg[0])

        # Tuple with all the related variable.
        return q, K, C, q_inf_max
    # _end_def_

# _end_class_
