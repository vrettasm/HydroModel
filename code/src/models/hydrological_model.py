import numpy as np
from copy import deepcopy


class HydrologicalModel(object):
    """
    This class represents a general (base) class of a hydrological model.
    All specific models must "inherit" from this class. In here we keep a
    copy of the necessary data (+ objects) that are required to construct
    a hydrological model.
    """

    def __init__(self, soil, porous, k_hc, theta_res, dz):
        """
        Default constructor.

        :param soil: soil properties object
        :param porous: porosity object
        :param k_hc: hydraulic conductivity object
        :param theta_res: water content (residual value)
        :param dz: spatial discretization [L: cm]
        """

        # Extract soil parameters.
        self.n = soil.n
        self.m = soil.m
        self.alpha = soil.alpha
        self.psi_sat = soil.psi_sat
        self.epsilon = np.maximum(soil.epsilon, 1.0e-8)

        # Copy the porosity and hydraulic conductivity objects.
        self.porous = deepcopy(porous)
        self.k_hc = deepcopy(k_hc)

        # Make a copy of the other input parameters.
        self.theta_res = theta_res
        self.dz = dz
    # _end_def_

    def pressure_head(self, theta, z):
        """
        Returns the pressure head at depth(s) 'z', given as input the volumetric
        water content theta.

        :param theta: volumetric water content.

        :param z: depth that we want to get the pressure head.

        :return: pressure head (suction) at depth(s) 'z'.
        """
        # Make sure the input 'psi' has at least shape (d, 1).
        if len(theta.shape) == 1:
            theta = np.reshape(theta, (theta.size, 1))
        # _end_if_

        # Make sure the input 'z' has shape (d, 1).
        if len(z.shape) == 1:
            z = np.reshape(z, (z.size, 1))
        # _end_if_

        # Get the dimensions of the input array.
        dim_d, dim_m = theta.shape

        # Check the input dimensions (of the vertical domain).
        if dim_d != z.size:
            raise ValueError(" Input size dimensions do not match:"
                             " {0} not equal to {1}".format(dim_d, z.size))
        # _end_if_

        # Get the porosity field at 'z'.
        porous_z, _, _ = self.porous(z)

        # Repeat if necessary (for vectorization).
        if dim_m > 1:
            porous_z = np.repeat(porous_z, dim_m, 1)
        # _end_if_

        # Pre-compute constant parameters.
        delta_s = porous_z - self.theta_res

        # Volumetric water content [-] can't drop below the minimum level
        # and also go above the porosity profile values at each depth 'z'.
        q = np.minimum(np.maximum(theta, self.theta_res), porous_z)

        # Compute the effective saturation (Se \in [0,1]).
        # (i.e. the "normalized" water content)
        s_eff = (q - self.theta_res) / delta_s

        # N.B: Here we do not allow 's_eff' to be equal to zero because
        # raising zero to the power of '-1/m' will result in Inf errors.
        s_eff = np.minimum(np.maximum(s_eff, self.epsilon), 1.0)

        # Check for saturated cells.
        # Here allow for some small errors to exist.
        id_sat = s_eff >= 0.99998

        # % Initialize return array.
        psi_z = np.zeros((dim_d, dim_m))

        # Not easily vectorized (because of the index "id_sat").
        for i in range(dim_m):
            # Extract the saturated locations of the i-th vector.
            j = id_sat[:, i]

            # Compute the pressure head (psi) on the unsaturated soil.
            psi_z[~j, i] = -((s_eff[~j, i] ** (-1.0 / self.m) - 1.0) ** (1.0 / self.n)) / self.alpha

            # Compute the pressure head (psi) on the saturated soil.
            psi_z[~j, i] = np.arange(0, np.sum(j)-1.0) * self.dz
        # _end_if_

        # SAFEGUARD:
        psi_z[~np.isfinite(psi_z)] = -1.0e+5

        # Pressure head (suction), Effective Saturation.
        return psi_z, s_eff
    # _end_def_

# _end_class_
