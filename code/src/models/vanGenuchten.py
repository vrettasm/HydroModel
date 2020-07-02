import numpy as np

from .hydrological_model import HydrologicalModel


class vanGenuchten(HydrologicalModel):
    """
    This class uses the 'van Genuchten' - Mualem model to create an object that computes
    the unsaturated hydraulic conductivity along with other parameters related to the
    unsaturated flow of water in the vadose zone, using Richards' PDE.

    -- REFERENCE --
    M. Th. van Genuchten (1980), A closed-form equation for predicting the hydraulic conductivity
    of unsaturated soils. Soil Science Society of America Journal 44, pp: 892-898.
    -- REFERENCE --

    """

    def __init__(self, soil, porous, k_sat, theta_res, dz):
        # Call the constructor of the parent class.
        super().__init__(soil, porous, k_sat, theta_res, dz)
    # _end_def_

    def __call__(self, psi, z, *args):
        """
        A direct call to an object of this class will return the water content, along other
        related quantities, at a specific depth 'z', given the input pressure head (suction).

        :param psi: pressure head (suction) [dim_d x dim_m].

        :param z: depth values (increasing downwards) [dim_d x dim_m].

        :param args: (for compatibility) to keep a uniform interface.

        :return: q (water content), K (unsaturated hydraulic conductivity),
        C (specific moisture capacity), Ksat (saturated hydraulic conductivity),
        and q_inf_max (max infiltration capacity [dim_d x dim_m].

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
        if dim_d != z.size:
            raise ValueError(" {0}: Input size dimensions do not match:"
                             " {1} not equal to {2}.".format(self.__class__.__name__, dim_d, z.size))
        # _end_if_

        # Create a vector with the K_{sat} values.
        k_sat = self.k_hc.sat_soil * np.ones(psi.shape)

        # Get the porosity field at 'z'.
        porous_z, *_ = self.porous(z)

        # Make sure the porosity is at least 1-D.
        porous_z = np.atleast_1d(porous_z)

        # Vectorized version.
        if dim_m:
            porous_z = porous_z.repeat(dim_m).reshape(dim_d, dim_m)
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

        # Pre-compute the $m-th$ root.
        mth_s_eff = s_eff ** (1.0 / self.m)

        # Compute the Unsaturated Hydraulic Conductivity.
        K = k_sat * np.sqrt(s_eff) * (1.0 - (1.0 - mth_s_eff) ** self.m) ** self.n

        # SAFEGUARD:
        K[id_sat] = k_sat[id_sat]
        K = np.minimum(K, k_sat)

        # Compute the Specific Moisture Capacity [dTheta/dPsi].
        C = (self.m * self.n) * self.alpha *\
            delta_s * (s_eff ** (1.0 / self.m + 1.0)) * (self.alpha * np.abs(psi)) ** (self.n - 1.0)

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
        q_inf_max = np.minimum(2.0*(porous_z[0] - q[0])*dz, k_sat[0])

        # Tuple with all the related variables.
        return np.atleast_1d(q, K, C, k_sat, q_inf_max)
    # _end_def_

# _end_class_
