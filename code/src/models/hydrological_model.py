import numpy as np


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

        :param soil: soil properties object.
        :param porous: porosity object.
        :param k_hc: hydraulic conductivity object.
        :param theta_res: water content (residual value).
        :param dz: spatial discretization [L: cm].
        """

        # Extract soil parameters.
        self.n = soil.n
        self.m = soil.m
        self.alpha = soil.alpha
        self.psi_sat = soil.psi_sat
        self.epsilon = np.maximum(soil.epsilon, 1.0e-8)

        # Copy the porosity and hydraulic conductivity objects.
        self.k_hc = k_hc
        self.porous = porous

        # Make a copy of the other input parameters.
        self.dz = dz
        self.theta_res = theta_res
    # _end_def_

    def pressure_head(self, theta, z):
        """
        Returns the pressure head at depth(s) 'z', given as input
        the volumetric water content 'theta'.

        :param theta: volumetric water content [dim_d x dim_m].

        :param z: depth that we want to get the pressure head [dim_d].

        :return: pressure head (suction) at depth(s) 'z' [dim_d x dim_m].

        :raises ValueError: if there is a mismatch in the input dimensions.
        """

        # Ensure the input is 1-D.
        z, theta = np.atleast_1d(z, theta)

        # Get the vector size.
        dim_d, dim_m = theta.shape[0], None

        # Check if the input is 2D.
        if theta.ndim == 2:
            dim_m = theta.shape[1]
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
        psi_z = np.zeros(theta.shape)

        # Compute the pressure head (psi) on the unsaturated soil.
        psi_z[~id_sat] = -((s_eff[~id_sat] ** (-1.0 / self.m) - 1.0) ** (1.0 / self.n)) / self.alpha

        # Compute the pressure head (psi) on the saturated soil.
        psi_z[id_sat] = np.arange(0, np.sum(id_sat)) * self.dz

        # SAFEGUARD:
        psi_z[~np.isfinite(psi_z)] = -1.0e+5

        # Pressure head (suction), Effective Saturation.
        return psi_z, s_eff
    # _end_def_

# _end_class_
