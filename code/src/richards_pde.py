import numpy as np
from scipy.integrate import solve_ivp

class RichardsPDE(object):
    """

    """
    __slots__ = ("m_data", "x_mesh", "x_mid", "mid_i", "xzmp1", "zxmp1",
                 "nx", "h_model", "ic0")

    def __init__(self, m_data=None):
        """
        Default constructor of the class.

        :param m_data: dictionary with all the simulation parameters.
        """
        # Check if we have been given input.
        if not m_data:
            raise ValueError(" No input is given. Richards' model cannot initialize.")
        # _end_if_

        # Make a reference to the input data structure.
        self.m_data = m_data

        # Set the initial conditions array.
        self.ic0 = m_data["init_cond"]

        # Get the hydrological model.
        self.h_model = m_data["hydro_model"]

        # # PREPARE THE MESH-GRIDS AND STORE THEM IN THE OBJECT # #

        # Spatial discretized grid [L x 1].
        self.x_mesh = self.m_data["z_grid"]

        # Get the number of discrete points.
        self.nx = self.x_mesh.size

        # Get the distance between space points.
        dx = np.diff(self.x_mesh, axis=0)

        # Make sure the space domain is increasing (downwards).
        if np.any(dx <= 0.0):
            raise RuntimeError(" Space domain is not increasing.")
        # _end_if_

        # Initialize the $nx-1$ mid-points where functions will be evaluated.
        self.x_mid = self.x_mesh[0:-1] + 0.5 * dx

        # Mid-point indexes.
        self.mid_i = np.arange(1, self.nx)

        # % Interior grid (mid)-points.
        self.xzmp1 = np.zeros(self.nx)
        self.xzmp1[1:] = self.x_mesh[1:self.nx] - self.x_mid
        self.zxmp1 = self.x_mid - self.x_mesh[0:-1]
    # _end_def_

    def __call__(self, t, y):
        """
        This function implements the space discretization of the PDE.
        It is vectorized for best performance, so the ode solver can
        use this property to speed up the solution.

        :param t:

        :param y:

        :return: dydt
        """

        # Make sure the input 'y' has at least shape (d, 1).
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
        # _end_if_

        # Get the number of different vectors.
        _, dim_m = y.shape

        # Preallocate for efficiency.
        dydt = np.zeros((self.nx, dim_m))

        # Evaluate the PDE at the top.
        y0, dy0 = midpoints(self.x_mesh[0], y[0, :],
                            self.x_mesh[1], y[1, :])
    # _end_def_

    # _end_def_

    def pde_fun(self, z, t, y, dydz, **kwargs):

        # Initial (dimensional) values.
        dim_d, dim_m = y.size, 1

        # Number of different columns.
        if len(y.shape) > 1:
            dim_d, dim_m = y.shape
        # _end_if_

        # Compute the hydraulic conductivity 'K(.)' and the specific moisture
        # capacity 'C(.)'.  Additionally return the soil moisture at the same
        # depths and the porosity. These are used for the root efficiency and
        # the hydraulic lift processes.
        theta, K, C, _ = self.h_model(y, z, **kwargs)

        # Compute the flux term.
        f = K*(np.minimum(dydz, 1.01) - 1.0)

        # Sink term is initialized to zero.
        s = np.zeros((dim_d, dim_m))

        # Make sure the length exceeds one cell.
        if dim_d > 1:
            # Get the Month (Mn) and the Hour (Hr) of the current time.
            Mn, Hr = kwargs["time_k"].month, kwargs["time_k"].hour

            # In Normal Mode: SPINUP == FALSE.
            if self.m_data["sim_flags"]["SPINUP"]:

                # Get the discretization step [L: cm]
                dz = self.m_data["dz"]

                # Compute day-light hours. We assume that the daylight is
                # between [06:00:00] (morning) and [17:59:59] (afternoon).
                daylight = Hr in np.arange(6, 18)

                # If either of this processes is active than
                # we need the rooting depth profiles at 'z'.
                if self.m_data["sim_flags"]["HLIFT"] or self.m_data["sim_flags"]["ET"]:
                    # Get the root density object.
                    root_pdf = self.m_data["root_pdf"]

                    # Find the indexes of the root-zone.
                    r_cells = z <= root_pdf.max_root_depth

                    # Get the roots density at depths (z).
                    roots_z = root_pdf(z)
                # _end_if_

                # Hydraulic Redistribution.
                # This runs ONLY in during the night-time!
                if self.m_data["sim_flags"]["HLIFT"] and (not daylight):
                    # Parameter for HR:
                    a_star = 1800

                    # Inverse of $\psi_{50}$:
                    inv_psi_50 = self.m_data["iPsi_50"]

                    # Saturated hydraulic conductance:
                    c_sat = a_star * self.m_data["env_param"]["Leaf_Area_Index"]

                    # Hydraulic conductance parameter:
                    c_HR = c_sat * ((1.0 - inv_psi_50 * y[r_cells, :]) ** 2) * roots_z

                    # Pressure difference.
                    dy = dydz * dz

                    # Update the flux 'f' with the new HR term. NB: We need to
                    # scale by '0.5' because the q_{HR} has [cm/h] and we need
                    # per 0.5h.
                    f[r_cells, :] += 0.5 * c_HR * dy[r_cells, :]
                # _end_if_

                # Evapo-transpiration (Tree Roots Water Uptake)
                # This runs ONLY during day-time!
                if self.m_data["sim_flags"]["ET"] and daylight:
                    pass
                # _end_if_

                # Lateral (subsurface) Runoff:
                if self.m_data["sim_flags"]["LF"]:
                    pass
                # _end_if_

            # _end_if_
        # _end_if_

        return C, s, f

    # _end_def_

    def ic_fun(self, y0, **kwargs):
        """
        This function can model the initial conditions vector.
        In this case is a simple assignment, but it can handle
        mode complex initializations if required.

        :param y0: initial conditions vector y(z,t=0) [dim_d x 1].

        :param kwargs: to pass additional parameters if requested.

        :return: initial conditions vector.
        """
        # Update the object status
        self.ic0 = y0

        # Return the vector.
        return self.ic0
    # _end_def_

    def bc_fun(self, z_left, y_left, z_right, y_right, **kwargs):
        """

        :param z_left:

        :param y_left:

        :param z_right:

        :param y_right:

        :param kwargs:

        :return:
        """

        # Get the maximum infiltration capacity from the model.
        theta_left, _, _, q_inf_max = self.h_model(y_left, z_left, **kwargs)

        # Initial value assumes 1-D.
        dim_m = 1

        # Number of different columns.
        if len(y_left.shape) > 1:
            _, dim_m = y_left.shape
        # _end_if_

        # [TOP]: Neumann boundary condition.
        q_left = np.ones(dim_m)

        # Initialize to zero flux.
        p_left = np.zeros(dim_m)

        # Extract the value of the forcing (rainfall).
        # (NB. Make sure is non-negative .... >= 0.0)
        q_rain_flux = np.abs(kwargs["precip_k"])

        # Net rain flux that enters the domain is given by:
        # >> netInputFlux = Precipitation - Interception.
        net_input_flux = (1.0 - kwargs["intercept_k"]) * q_rain_flux

        # If the top cells are not saturated set an upper bound.
        p_left[y_left < self.m_data["soil"].psi_sat] = np.minimum(net_input_flux, q_inf_max)

        # [BOTTOM]: Neumann condition is set to zero flux.
        p_right = np.zeros(dim_m)
        q_right = np.ones(dim_m)

        # Return the boundary values.
        return p_left, q_left, p_right, q_right
    # _end_def_

    def solve(self, t, y0, rtol=1.0e-4, atol=1.0e-4):
        """
        This function numerically integrates a system of ordinary differential
        equations (self), within the time span "t", given an initial value y0.

        :param t: time window to integrate the ode (t0, tf).

        :param y0: initial conditions array [L x 1].

        :param z: spatial domain at which we want to solve the equation. [L x 1].

        :param rtol: absolute tolerance.

        :param atol: absolute tolerance.

        :return: dydt [L x n]
        """

        # Trials counter
        n_trials = 5

        # Try to solve the interval "n_trials" times.
        while n_trials > 0:
            # Current solution
            sol_t = solve_ivp(self, t_span=t, y0=y0,
                              method='LSODA', atol=atol, rtol=rtol)

            # Check if the solver terminated successfully.
            if sol_t.success:
                # Exit the loop here.
                break
            else:
                # Reduce the noise in the model.
                # Reduce the noise in the model.
                # Reduce the noise in the model.
                # Reduce the noise in the model.

                # Decrease counter by one.
                n_trials -= 1

                # If we reach here the solver failed to integrate successfully.
                if n_trials == 0:
                    print(" The solver failed with message: {0}".format(sol_t.message))
                # _end_if_
            # _end_if_

        # _end_while_

        # Return the (integrated) time derivative: dydt.
        return sol_t.y[:, -1]
    # _end_def_

# _end_class_

# Helper module method.
def midpoints(x_left, u_left, x_right, u_right):
    """
    Interpolation (helper) function that is used
    in the discretization of the Richards' PDE.

    :param x_left: [dim_d x 1]

    :param u_left: [dim_d x dim_m]

    :param x_right: [dim_d x 1]

    :param u_right: [dim_d x dim_m]

    :return: derivatives and the mid-points.
    """

    # Check for input mis-match.
    if u_left.shape != u_right.shape:
        raise RuntimeError(" midpoints: Input 'u' dimensions do not match."
                           " {0} not equal to {1}".format(u_left.shape, u_right.shape))
    # _end_if_

    # Check for input mis-match.
    if x_left.shape != x_right.shape:
        raise RuntimeError(" midpoints: Input 'x' dimensions do not match."
                           " {0} not equal to {1}".format(x_left.shape, x_right.shape))
    # _end_if_

    # Use a simple (arithmetic) average to approximate
    # the 'mid-points' between 'u_left' and 'u_right'.
    u_mid = 0.5 * (u_left + u_right)

    # Spacings between the two input points.
    # This will be either scalar, or vector.
    dx = x_right - x_left

    # Vectorization code.
    if len(u_mid.shape) > 1:

        # Get the number of columns.
        _, m = u_mid.shape

        # If the x-input is vector.
        # !! Otherwise is scalar !!
        if len(dx.shape) == 1:

            # Convert it to 1d array.
            dx = dx.reshape(-1, 1)

            # Repeat if necessary.
            if m > 1:
                # These will be (d, m) arrays.
                dx = np.repeat(dx, m, 1)
            # _end_if_

        # _end_if_

    # _end_if_

    # Central Difference Formula:
    du_mid = (u_right - u_left)/dx

    # Return the derivative and the mid-points.
    return u_mid, du_mid

# _end_def_
