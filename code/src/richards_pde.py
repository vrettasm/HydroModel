import numpy as np
from scipy.integrate import solve_ivp
from code.src.utilities import find_wtd

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

        # Evaluate the PDE at the top (1/2).
        y0, dy0 = midpoints(self.x_mesh[0], y[0, :],
                            self.x_mesh[1], y[1, :])

        # Evaluate the PDE at the top (2/2).
        cL, sL, fL, *_ = self.pde_fun(self.x_mid[0], y0, dy0)

        # Evaluate the boundary conditions.
        pL, qL, pR, qR = self.bc_fun(self.x_mesh[0], y[0, :], self.x_mesh[-1], y[-1, :])

        # TOP BOUNDARY:
        if qL == 0:
            dydt[0, :] = pL
        else:
            # Compute the contribution of C(.):
            denom = qL * self.zxmp1[0, :] * cL

            # Avoid division by zero.
            denom[denom == 0.0] = 1.0

            # Compute the derivative at $z = 0$:
            dydt[0, :] = (pL + qL * (fL + self.zxmp1[0, :] * sL)) / denom
        # _end_if_

        # INTERIOR POINTS (vectorized computation):
        y_mid, dy_mid = midpoints(self.x_mesh[self.x_mid+0], y[self.x_mid+0, :],
                                  self.x_mesh[self.x_mid+1], y[self.x_mid+1, :])
        # PDE evaluation.
        cR, sR, fR, *_ = self.pde_fun(self.x_mid(self.x_mid), y_mid, dy_mid)

        # WARNING: DO NOT EDIT THESE LINES
        cLi = np.append(cL, cR, axis=0)
        fLi = np.append(fL, fR, axis=0)
        sLi = np.append(sL, sR, axis=0)
        # WARNING: DO NOT EDIT THESE LINES

        # Compute the contribution of C(.):
        denom = self.zxmp1[self.x_mid] * cR + self.xzmp1[self.x_mid] * cLi[0:-2]

        # Avoid division by zero.
        denom[denom == 0.0] = 1.0

        # Compute the derivatives at $z = [1:-2]$:
        dydt[self.x_mid, :] = (fR - fLi[0:-2, :] + (self.zxmp1[self.x_mid] * sR +
                                                    self.xzmp1[self.x_mid] * sLi[0:-2, :])) / denom
        # BOTTOM BOUNDARY:
        if qR == 0.0:
            dydt[-1, :] = pR
        else:
            # Compute the contribution of C(.):
            denom = -qR * self.xzmp1[-1] * cLi[-1, :]

            # Avoid division by zero.
            denom[denom == 0.0] = 1.0

            # Compute the derivative at $z = end$:
            dydt[-1, :] = (pR + qR * (fLi[-1, :] - self.xzmp1[-1] * sLi[-1, :])) / denom
        # _end_if_

        # Return the derivative.
        return dydt
    # _end_def_

    def pde_fun(self, z, y, dydz, **kwargs):
        """

        Richards' Equation (PDE - 1d).

        :param z:

        :param y:

        :param dydz:

        :param kwargs:

        :return:
        """

        # Make sure the 'y' input is at least (dim_d x 1)
        if len(y.shape) == 1:
            y = y.reshape(-1, 1)
        # _end_if_

        # Get the input dimensions.
        dim_d, dim_m = y.shape

        # Compute the hydraulic conductivity 'K(.)' and the specific moisture
        # capacity 'C(.)'.  Additionally return the soil moisture at the same
        # depths and the porosity. These are used for the root efficiency and
        # the hydraulic lift processes.
        theta, K, C, _ = self.h_model(y, z, **kwargs)

        # Compute the flux term.
        flux = K*(np.minimum(dydz, 1.01) - 1.0)

        # Sink term is initialized to zero.
        sink = np.zeros((dim_d, dim_m))

        # Define some output variables.
        tran_k, latf_k = None, None

        # Make sure the length exceeds one cell.
        if dim_d > 1:
            # Get the Month (Mn) and the Hour (Hr) of the current time.
            Mn, Hr = kwargs["time_k"].month, kwargs["time_k"].hour

            # In Normal Mode: SPINUP == FALSE.
            if not self.m_data["sim_flags"]["SPINUP"]:

                # Get the discretization step [L: cm]
                dz = self.m_data["dz"]

                # Compute day-light hours. We assume that the daylight is
                # between [06:00:00] (morning) and [17:59:59] (afternoon).
                daylight = Hr in np.arange(6, 18)

                # If either of this processes is active than
                # we need the rooting depth profiles at 'z'.
                if self.m_data["sim_flags"]["HLIFT"] or self.m_data["sim_flags"]["ET"]:
                    # Get the root density object.
                    tree_roots = self.m_data["tree_roots"]

                    # Find the indexes of the root-zone.
                    r_cells = z <= tree_roots.max_root_depth

                    # Get the roots density at depths (z).
                    roots_z = tree_roots(z)
                # _end_if_

                # Hydraulic Redistribution.
                # This runs ONLY in during the night-time!
                if self.m_data["sim_flags"]["HLIFT"] and (not daylight):
                    # Parameter for HR:
                    a_star = 1800

                    # Inverse of $\psi_{50}$:
                    inv_psi_50 = self.m_data["iPsi_50"]

                    # Saturated hydraulic conductance:
                    c_sat = a_star * self.m_data["Trees"]["Leaf_Area_Index"]

                    # Hydraulic conductance parameter:
                    c_HR = c_sat * ((1.0 - inv_psi_50 * y[r_cells, :]) ** 2) * roots_z

                    # Pressure difference.
                    dy = dydz * dz

                    # Update the flux 'f' with the new HR term. NB: We need to
                    # scale by '0.5' because the q_{HR} has [cm/h] and we need
                    # per 0.5h.
                    flux[r_cells, :] += 0.5 * c_HR * dy[r_cells, :]
                # _end_if_

                # Evapo-transpiration (Tree Roots Water Uptake)
                # This runs ONLY during day-time!
                if self.m_data["sim_flags"]["ET"] and daylight:
                    # Compute the root efficiency.
                    rho_theta, Wk = tree_roots.root_efficiency(theta[r_cells, :], z[r_cells])

                    # Get the product of the two terms.
                    x_out = rho_theta * roots_z

                    # Compute the integral of:
                    # $\int_{z} (root_efficiency x root_fraction) dz$
                    tot_x = np.sum(x_out, axis=0) * dz

                    #  Constraint No.2:
                    if np.all(tot_x > 1.0):
                        # Normalize it so it's equal to one.
                        x_out = x_out / tot_x

                        # Recompute the integral.
                        tot_x = np.sum(x_out, axis=0) * dz
                    # _end_if_

                    # Compute  the transpiration parameter only if there
                    # is not root shutdown! Otherwise don't remove water.
                    if np.all(tot_x > 0.0):
                        # Total transpiration demand.
                        tot_tr = np.minimum(kwargs["atm"], Wk)

                        # Compute the rate of potential transpiration
                        # by dividing the total atmospheric demand,at
                        # the current time, with the integrated value
                        # of the (root_efficiency x root_fraction).
                        tr_pot = tot_tr / tot_x

                        # Compute the uptake distribution.
                        uptake = tr_pot * x_out

                        # Remove the root uptake from the sink term.
                        sink[r_cells, :] = -uptake

                        # Store the transpiration as function of depth.
                        tran_k = uptake
                    # _end_if_
                # _end_transpiration_if_

                # Lateral (subsurface) Runoff:
                if self.m_data["sim_flags"]["LF"]:
                    # Find saturated cell indexes.
                    id_sat = (y >= self.m_data["soil"].psi_sat)

                    # Switch according to the running mode.
                    if self.m_data["sim_flags"]["PREDICT"]:
                        # Predictive (running) mode.
                        # Set the sink coefficient accordingly:
                        if Mn in [10, 11, 12, 1, 2, 3]:
                            # Wet season coefficient.
                            alpha_low = -2.5e-3
                        else:
                            # Dry season coefficient.
                            alpha_low = -1.5e-3
                        # _end_if_

                        # Number of cells, from the bottom of the well, that stay
                        # always saturated. That number can vary from well to well
                        # and from year to year. This is to prevent the spatial
                        # domain from draining completely during extended droughts.
                        sat_cells = kwargs["sat_cells"] - 1

                        # Auxiliary variable.
                        low_lim = dim_d - sat_cells

                        # Exponent range.
                        nu = np.linspace(1.5, 0.0, low_lim)

                        # Repeat for all multiple entries of 'y'.
                        for i in np.range(dim_m):
                            # Find "wtd_est" (index).
                            wtd_est = find_wtd(id_sat[:, i])

                            # If the "wtd_est" is above a pre-defined depth value.
                            if wtd_est < low_lim:
                                # Set the cell indexes that we remove water from.
                                j = np.arange(wtd_est, wtd_est + 1)

                                # Update scale with current estimate of the wtd.
                                # The higher the exponent the faster it drains!!
                                # alpha_lat = alpha_low * (low_lim - wtd_est) / low_lim
                                alpha_lat = alpha_low * (1.0 - (j / low_lim) ** nu[j])

                                # Compute the sink term proportional to y(z,t).
                                sink[j, i] = np.minimum(alpha_lat * y[j, i], sink[j, i])

                                # Store the lateral flow (runoff), along with
                                # the locations (indexes) in the space domain.
                                latf_k = (j, np.abs(sink[j, 1]))
                            # _end_if_
                        # _end_for_
                    else:
                        # Monitoring (running) mode.
                        # Current water table observation.
                        wtd_obs = np.minimum(kwargs["WTDi"], dim_d-1)

                        # Sink (scale) coefficient:
                        alpha_lat = -2.5e-4

                        # Repeat for all multiple entries of 'y'.
                        for i in range(dim_m):
                            # Find "wtd_est" (index).
                            wtd_est = find_wtd(id_sat[:, i])

                            # If the "wtd_est" is inside the space domain.
                            if (wtd_est < dim_d) & (wtd_est < wtd_obs):
                                # Find  the locations between the estimated
                                # $wtdEst$ and the actual value of $wtdObs$.
                                j = np.arange(wtd_est, wtd_obs)

                                # Compute the sink term proportional to y(z,t).
                                sink[j, i] = np.minimum(alpha_lat * y[j,i], sink[j, i])

                                # Store the lateral flow (runoff), along with
                                # the locations (indexes) in the space domain.
                                latf_k = (j, np.abs(sink[j, 1]))
                            # _end_if_
                        # _end_for_
                    # _end_if_
                # _end_if_
            # _end_if_
        # _end_if_

        # Exit:
        return C, sink, flux, tran_k, latf_k
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

    def solve(self, t, y0, r_tol=1.0e-4, a_tol=1.0e-4):
        """
        This function numerically integrates a system of ordinary differential
        equations (self), within the time span "t", given an initial value y0.

        :param t: time window to integrate the ode (t0, tf).

        :param y0: initial conditions array [L x 1].

        :param z: spatial domain at which we want to solve the equation. [L x 1].

        :param r_tol: absolute tolerance.

        :param a_tol: absolute tolerance.

        :return: dydt [L x n]
        """

        # Trials counter
        n_trials = 5

        # Try to solve the interval "n_trials" times.
        while n_trials > 0:
            # Current solution
            sol_t = solve_ivp(self, t_span=t, y0=y0,
                              method='LSODA', atol=a_tol, rtol=r_tol)

            # Check if the solver terminated successfully.
            if sol_t.success:
                # Exit the loop here.
                break
            else:
                # Reduce the noise in the model.
                self.m_data["n_rnd"] *= 0.8

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
