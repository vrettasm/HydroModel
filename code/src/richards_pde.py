import numpy as np
from scipy.integrate import solve_ivp

class RichardsPDE(object):
    """

    """
    __slots__ = ("m_data", "x_mesh", "x_mid", "mid_i",
                 "xzmp1", "zxmp1", "nx")

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

    # This is the discretePDEs(t, y)
    def __call__(self, t, y):
        """

        :param t:

        :param y:

        :return:
        """

        # Get the number of different vectors.
        _, dim_m = y.shape

        # Preallocate for efficiency.
        dydt = np.zeros((self.nx, dim_m))

    # _end_def_

    @staticmethod
    def pde_mid(x_left, u_left, x_right, u_right):
        """

        :param x_left:

        :param u_left:

        :param x_right:

        :param u_right:

        :return:
        """

        # Use a simple (arithmetic) average to approximate
        # the mid-points between 'u_left' and 'u_right'.
        u_mid = 0.5*(u_left + u_right)

        # Spacing between the two input points.
        dx = x_right - x_left

        # Get the number of columns.
        _, m = u_left.shape

        # Repeat if necessary (for vectorization).
        if m > 1:
            # These will be (d, m) arrays.
            dx = np.repeat(dx, m, 1)
        # _end_if_

        # Central Difference Formula:
        # u(i+0.5*h) - u(i-0.5*h):
        du_mid = (u_right - u_left)/dx

        # Return the derivative and the mid-points.
        return u_mid, du_mid
    # _end_def_

    def pde_fun(self):
        pass
    # _end_def_

    def ic_fun(self):
        pass
    # _end_def_

    def bc_fun(self):
        pass
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
