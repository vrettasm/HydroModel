import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import gamma as sc_gamma

class TreeRoots(object):
    """
    This class models tree roots. It computes the roots density for a given
    depth and also returns roots efficiency, as described in refs. [1, 2].

    To have more consistent results, the root depth in here is mapped to a
    predefined range [1 to 100].

    -- REFERENCE(S) --
    [1] Chun-Ta Lai and Gabriel Katul (2000): The dynamic role of root-water
    uptake in coupling potential to actual transpiration. Advances in Water
    Resources, Vol. 23, pages: 427-439.

    [2] D. B. G. Collins and R. L. Bras (2007): Plant rooting strategies in
    water-limited ecosystems. Water Resources Research, Vol. 43, W06407,
    doi:10.1029/2006WR005541, pages: 1-10.
    -- REFERENCE(S) --

    """

    __slots__ = ("profile", "r_model", "max_depth_cm", "f_interp", "dz", "porous")

    def __init__(self, ln, dz, r_model, porous):
        """
        Default constructor of the root density object class.

        :param ln: number of discrete cells that form the total root zone.

        :param dz: space discretization [L: cm].

        :param r_model: root density model (type).
        """

        # Store the space discretization.
        self.dz = dz

        # Store the root density model.
        self.r_model = r_model

        # Store the porosity object.
        self.porous = porous

        # Store the maximum root depth [L: cm].
        self.max_depth_cm = (ln * dz)

        # Discrete cells with uniform 'dz' spacing.
        # This maps the depth in: [1-100].
        n_cells = np.linspace(1, 100, ln)

        # Create the profile according to the selected type.
        if str.upper(r_model) == "UNIFORM":

            root_pdf = np.ones(ln)/ln

        elif str.upper(r_model) == "NEGATIVE_EXP":

            # [WARNING]: The range should be between [1..20].
            # Exponential (mean) parameter:
            mu = 15

            # Point to evaluate the exponential function.
            kp = n_cells/mu

            # Set the exponential profile.
            root_pdf = np.exp(-kp)/mu

        elif str.upper(r_model) == "GAMMA_PDF":

            # [WARNING]: The range should be between [1..10].
            # Shape parameter 'a':
            shape = 2.5

            # Scale parameter 'b':
            scale = 5.0

            # The gamma 'pdf' is given by (in LaTex format):
            # $f(x|a,b) = \frac{1}{b^a Gamma(a)} x^{a-1}e^{-x/b}$
            root_pdf = sc_gamma.pdf(n_cells, a=shape, scale=scale)

        elif str.upper(r_model) == "MIXTURE":

            # -(A)-
            # [WARNING]: The range should be between [1..20].
            # Exponential (mean) parameter:
            mu = 15

            # Point to evaluate the exponential function.
            kp = n_cells/mu

            # Set the exponential profile.
            root_a = np.exp(-kp)/mu

            # -(B)-
            # [WARNING]: The range should be between [1..10].
            # Shape parameter 'a':
            shape = 2.5

            # Scale parameter 'b':
            scale = 5.0

            # The gamma 'pdf' is given by (in LaTex format):
            # $f(x|a,b) = \frac{1}{b^a Gamma(a)} x^{a-1}e^{-x/b}$
            root_b = sc_gamma.pdf(n_cells, a=shape, scale=scale)

            # Prior of the first pdf.
            w1 = 0.15

            # Prior of the second pdf.
            w2 = 0.85

            # Apply mixture model.
            root_pdf = (w1 * root_a) + (w2 * root_b)

        else:
            raise ValueError(" Wrong root density profile type: {0}".format(r_model))
        # _end_if_

        # Cut-off very small values before the normalization.
        root_pdf = np.maximum(root_pdf, 1.0e-10)

        # Compute the integral of the root profile term.
        total_area = np.sum(root_pdf)*dz

        # Assign the normalized root profile.
        self.profile = root_pdf/total_area

        # Interpolate the profile on the total root depth [L: cm].
        self.f_interp = interp1d(np.linspace(0.0, ln * dz, ln), self.profile)
    # _end_def_

    def __call__(self, z_new=None):
        """
        A call of the tree roots object will return the root pdf profile for a
        specific input 'z_new'. If no new depth values are given, the call will
        return the full profile.

        :param z_new: the depth(s) at which we want the new profile [L: cm].

        :return: the density profile at depth(s) 'z'.
        """

        # Check if we have given input.
        if z_new is None:
            # Return the profile (over the whole z-domain).
            return self.profile
        else:
            # Return the profile at the new input 'z_new'.
            return self.f_interp(z_new)
        # _end_if_
    # _end_def_

    def efficiency(self, theta_z, z_roots):
        """
        This method computes the root efficiency as the product of two terms:

        ( i) alpha01 -> maximum efficiency, when the water is not limiting,
        (ii) alpha02 -> root shutdown.

        The method is vectorized for faster execution, when the input matrix
        "theta" has multiple columns [dim_d x dim_m].

        :param theta_z: soil moisture values at depth(s) 'z'. [dim_d x dim_m]

        :param z_roots: depth values at which we want to return the root efficiency. [dim_d x 1]

        :return: (1) root efficiency (alpha01 * alpha02) [dim_d x dim_m]
                 (2) available water in the root zone    [1 x 1]
        """

        # Make sure the theta input is at least (dim_d x 1)
        if len(theta_z.shape) == 1:
            theta_z = theta_z.reshape(-1, 1)
        # _end_if_

        # Get the input dimensions.
        dim_d, dim_m = theta_z.shape

        # Make sure the dimensions match.
        if dim_d != z_roots.size:
            raise RuntimeError(" Input dimensions do not match:"
                               " {0} not equal to {1}".format(dim_d, z_roots.size))
        # _end_if_

        # Compute porosity, field capacity and wilting points at 'z'.
        porous_z, f_cap_z, wlt_z = self.porous(z_roots)

        # Make sure the porosity is at least (dim_d x 1)
        if len(porous_z.shape) == 1:
            porous_z = porous_z.reshape(-1, 1)
        # _end_if_

        # Make sure the field capacity is at least (dim_d x 1)
        if len(f_cap_z.shape) == 1:
            f_cap_z = f_cap_z.reshape(-1, 1)
        # _end_if_

        # Make sure the wilting point is at least (dim_d x 1)
        if len(wlt_z.shape) == 1:
            wlt_z = wlt_z.reshape(-1, 1)
        # _end_if_

        # Replicate for vectorization.
        if dim_m > 1:
            porous_z = porous_z.repeat(dim_m, 1)
            f_cap_z = f_cap_z.repeat(dim_m, 1)
            wlt_z = wlt_z.repeat(dim_m, 1)
        # _end_if_

        # Compute the available water (for extraction) in the root-zone.
        water_k = np.sum(theta_z - wlt_z, axis=0) * self.dz

        # Check if there is water that can be extracted.
        if np.all(water_k > 0.0):
            # Compute the cumulative sum of - $\Sum_{z=0}^{zb} \theta(z)$, from the
            # top of the surface (z = 0), until the depth of the root zone (z = zb).
            local_z = np.cumsum(theta_z, axis=0) * self.dz

            # Extract the last row (corresponds to the total integral).
            total_z = local_z[-1, :]

            # Safeguard: Avoid division by zero.
            total_z[total_z == 0.0] = 1.0

            # Replicate the vector: [dim_d x dim_m]
            total_z = np.ones((dim_d, dim_m)) * total_z

            # Compute the denominator of the first term in alpha_01.
            # [dim_d x 1]
            denom_01 = porous_z - wlt_z

            # Safeguard: Avoid division by zero.
            denom_01[denom_01 == 0.0] = 1.0

            # Maximum efficiency when water is not limiting.
            alpha_01 = np.maximum(theta_z/denom_01, local_z/total_z)

            # Preallocate root shutdown (default is zero).
            # Satisfies the cases where: $\theta <= \theta_{wlt}$
            alpha_02 = np.zeros((dim_d, dim_m))

            # Find the indexes (in the soil moisture),
            # where $\theta_{wlt} < \theta <= \theta_{flc}$
            cond_01 = (wlt_z < theta_z) & (theta_z <= f_cap_z)

            # Compute the denominator in alpha_02.
            denom_02 = f_cap_z[cond_01] - wlt_z[cond_01]

            # Safeguard: (this should never happen).
            denom_02[denom_02 == 0.0] = 1.0

            # Compute the root extraction for each layer.
            alpha_02[cond_01] = (theta_z[cond_01] - f_cap_z[cond_01])/denom_02

            # The final condition checks whether $theta(z) > theta_{flc}$,
            # and sets the values of 'alpha_02' to one.
            alpha_02[theta_z > f_cap_z] = 1.0

            # Safeguard: Make sure [0.0 <= alpha_02 <= 1.0].
            alpha_02 = np.minimum(np.maximum(alpha_02, 0.0), 1.0)

            # Here we make sure that if all the root zone is flooded with water at
            # full capacity  (theta(z) > theta_{flc}), then we artificially reduce
            # the ET demand to prevent the roots from drowning.
            alpha_02[:, np.all(alpha_02 == 1.0)] = 0.01

            # Root efficiency for each layer: $\rho(\theta(z))$
            rho_theta = np.abs(alpha_01 * alpha_02)

            # Compute the integral of the root efficiency term.
            tot_theta = np.sum(rho_theta, axis=0) * self.dz

            # Constraint No.1: normalize it so the integral is equal to one.
            rho_theta = rho_theta / tot_theta
        else:
            # Here we set everything to zero.
            water_k = 0.0
            rho_theta = np.zeros((dim_d, dim_m))
        # _end_if_

        return rho_theta, water_k
    # _end_def_

    def __str__(self):
        """
        Override to print a readable string presentation of the Tree Root object.
        At the moment this will include its id(), along with its profile type.

        :return: a string representation of a Root density object.
        """
        return " Tree-Roots Id({0}): Type={1},"\
               " Max-Depth={2} [cm]".format(id(self), self.r_model, self.max_depth_cm)
    # _end_def_

    @property
    def max_root_depth(self):
        """
        Accessor of the maximum root depth [L: cm].
        :return: the depth in centimeters.
        """
        return self.max_depth_cm
    # _end_def_

    @property
    def root_type(self):
        """
        Accessor of the root density model type.
        :return: the modelled pdf type (string).
        """
        return self.r_model
    # _end_def_

# _end_class_
