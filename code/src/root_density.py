import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import gamma as sc_gamma

class RootDensity(object):
    """
    This function computes the root density for a given depth of size 'L'.
    To have more consistent results, here it maps every depth to the range
    of a pre defined range of 1 to 100.
    """

    __slots__ = ("profile", "r_model", "max_depth_m", "f_interp")

    def __init__(self, ln, dz, r_model):
        """
        Default constructor of the root density object class.

        :param ln: number of discrete cells that form the total root zone.

        :param dz: space discretization [L: cm].

        :param r_model: root density model (type).

        """

        # Store the root density model.
        self.r_model = r_model

        # Store the maximum root depth [L: m].
        self.max_depth_m = (ln * dz)/100.0

        # Discrete cells with uniform 'dz' spacing.
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
            root_pdf = w1 * root_a + w2 * root_b

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
        A call of the root density object will return the root pdf profile for a
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

    def __str__(self):
        """
        Override to print a readable string presentation of the Root density object.
        At the moment this will include its id(), along with its profile type.

        :return: a string representation of a Root density object.
        """
        return " Root density Id({0}): Type={1},"\
               " Max-Depth={2} [m]".format(id(self), self.r_model, self.max_depth_m)
    # _end_def_

    @property
    def max_root_depth(self):
        """
        Accessor of the maximum root depth [L: m].
        :return: the depth in meters.
        """
        return self.max_depth_m
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
