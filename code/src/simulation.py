import json
import time
import h5py
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.random import default_rng, SeedSequence

from .hydraulic_conductivity import HydraulicConductivity
from .models.vanGenuchten import vanGenuchten
from .models.vrettas_fung import VrettasFung
from .porosity import Porosity
from .richards_pde import RichardsPDE
from .soil_properties import SoilProperties
from .tree_roots import TreeRoots
from .utilities import find_wtd
from .water_content import WaterContent


class Simulation(object):
    """
    Main simulation class. The normal workflow is as follow:

    1) Create a simulation object. You can also give a name
    that will be used when saving the data.

        >> sim_01 = Simulation("Sim_01")

    2) Setup its parameters and water data (precipitation, etc.).

        >> sim_01.setupModel(params, water_data)

    3) Run the simulation (this step might take a while).

        >> sim_01.run()

    4) Finally save the results in a hdf5 (compressed) file.

        >> sim_01.saveResults()
    """

    __slots__ = ("name", "mData", "rng", "pde_model", "output")

    def __init__(self, name=None, seed=None):
        """
        Default constructor of the Simulation class.

        :param name: (string) is optional but it will be used for constructing
        a meaningful filename to save the results at the end of the simulation.

        :param seed: (int) is used to initialize the random generator. If None,
        the generator will be initialized at random from the OS.
        """
        # Check if a simulation name has been given.
        if name is not None:
            self.name = name
        else:
            self.name = "ID_None"
        # _end_if_

        # This dictionary will hold all the simulation data.
        self.mData = {}

        # Create a random number generator.
        if seed is not None:
            self.rng = default_rng(SeedSequence(seed))
        else:
            self.rng = default_rng()
        # _end_if_

        # Place holder for the pde model.
        self.pde_model = None

        # Place holder for the output storage.
        self.output = {}
    # _end_def_

    def setupModel(self, params, data):
        """
        This method is called BEFORE the run() and sets up all the variables for
        the simulation.  It is also responsible for checking the validity of the
        input parameters before use.

        :param params: (dict) contains all the given parameters. If None, then it
        will use the default parameters to initialize object. However this is not
        recommended since the user will have no control on the parameters.

        :param data: (pandas.DataFrame) that will contain the water data for the
        simulation. These include, time-series of the dates, water table depths
        precipitation values, etc.

        :return: None.

        :raises ValueError: if some data are missing.

        :raises RuntimeError: TBD.
        """

        # The spacing here defines a Uniform-Grid [L: cm].
        dz = 5.0

        # Copy the well number that is being simulated.
        self.mData["Well_No"] = params["Well_No"]

        # Open the file in "Read Only" mode.
        with open(Path(params["Site_Information"]), 'r') as site_file:
            # Load the site information.
            site_info = json.load(site_file)
        # _end_with_

        # Check if the Well number exists in the site info file.
        if not str(self.mData["Well_No"]) in site_info["Well"]:
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The selected well does not exist in the site information file.")
        # _end_if_

        # Extract the Well information.
        well = site_info["Well"][str(self.mData["Well_No"])]

        # Make sure the well is not defined as fully saturated.
        if well["sat_depth"] >= well["max_depth"]:
            raise RuntimeError(f" {self.__class__.__name__}:"
                               f" The well seems fully saturated.")
        # _end_if_

        # Compute the number of continuously saturated cell (from the bottom).
        self.mData["sat_cells"] = np.ceil(well["sat_depth"] / dz)

        # Spatial domain parameters for the underground layers [L:cm]:
        layers = (well["soil"], well["saprolite"], well["weathered"], well["max_depth"])

        # Add the underground layers to the structure.
        self.mData["layers"] = layers

        # Add the spatial discretization to the structure.
        self.mData["dz"] = dz

        # Create a vertical grid (increasing downwards)
        z_grid = np.arange(well["soil"], well["max_depth"] + dz, dz)

        # Add the spatial domain to the structure.
        self.mData["z_grid"] = z_grid

        # Create a soil properties object.
        try:
            soil = SoilProperties(params["Soil_Properties"]["n"],
                                  params["Soil_Properties"]["a0"],
                                  params["Soil_Properties"]["psi_sat"],
                                  params["Soil_Properties"]["epsilon"])
        except Exception as e0:
            # Display the error message.
            print(f" SoilProperties failed to initialize: {e0}."
                  f" It will use default initialization parameters.")

            # Default initialization.
            soil = SoilProperties()
        # _end_try_

        # Add soil properties to the dictionary.
        self.mData["soil"] = soil

        # Create a water content object.
        try:
            theta = WaterContent(params["Water_Content"]["Theta_Min"],
                                 params["Water_Content"]["Theta_Max"],
                                 params["Water_Content"]["Theta_Residual"],
                                 params["Water_Content"]["Wilting_Point_cm"],
                                 params["Water_Content"]["Field_Capacity_cm"])
        except Exception as e0:
            # Display the error message.
            print(f" WaterContent failed to initialize: {e0}."
                  f" It will use default initialization parameters.")

            # Default initialization.
            theta = WaterContent()
        # _end_try_

        # Add water content to the dictionary.
        self.mData["theta"] = theta

        # Create a hydraulic conductivity object.
        try:
            K = HydraulicConductivity(params["Hydraulic_Conductivity"]["Sat_Soil"],
                                      params["Hydraulic_Conductivity"]["Sat_Saprolite"],
                                      params["Hydraulic_Conductivity"]["Sat_Fresh_Bedrock"],
                                      params["Hydraulic_Conductivity"]["Sigma_Noise"],
                                      params["Hydraulic_Conductivity"]["Lambda_Exponent"])
        except Exception as e0:
            # Display the error message.
            print(f" HydraulicConductivity failed to initialize: {e0}."
                  f" It will use default initialization parameters.")

            # Default initialization.
            K = HydraulicConductivity()
        # _end_try_

        # Add hydraulic conductivity to the dictionary.
        self.mData["K"] = K

        # Add the environmental parameters.
        self.mData["env_param"] = params["Environmental"]

        # Add the simulation (execution) flags.
        self.mData["sim_flags"] = params["Simulation_Flags"]

        # Create a porosity object.
        porous = Porosity(z_grid, layers, theta, soil,
                          params["Hydrological_Model"]["Porosity_Profile"])

        # Add it to the dictionary.
        self.mData["porosity"] = porous

        # Create and add the root density object.
        self.mData["tree_roots"] = TreeRoots(np.ceil(params["Trees"]["Max_Root_Depth_cm"]/dz),
                                             dz, params["Trees"]["Root_Pdf_Profile"], porous)

        # Add the selected hydrological model.
        if str.upper(params["Hydrological_Model"]["Name"]) == "VRETTAS_FUNG":
            # Create the hydro-object.
            self.mData["hydro_model"] = VrettasFung(soil, porous, K, theta.res, dz)

            # Print a message.
            print(" Selected model: Vrettas-Fung")
        else:
            # Create the hydro-object.
            self.mData["hydro_model"] = vanGenuchten(soil, porous, K, theta.res, dz)

            # Print a message.
            print(" Selected model: vanGenuchten")
        # _end_if_

        # Extract the observational data from the pandas.Dataframe:
        r_datenum = data.loc[:, "Datenum"]

        # We need to convert the dates (from MATLAB to Python).
        # NOTE: The value 719529 is MATLAB's datenum value of the "Unix epoch"
        # start (1970-01-01), which is the default origin for pd.to_datetime().
        timestamps = pd.to_datetime(r_datenum - 719529, unit='D')

        # Store timestamps in the dictionary.
        self.mData["time"] = [t.round(freq="s") for t in timestamps]

        # Convert the meters to [L:cm] before storing them.
        z_wtd_cm = np.array(np.abs(np.round(100.0 * data.loc[:, "WTD_m"])))

        # Check if there are NaN values.
        if np.any(np.isnan(z_wtd_cm)):
            raise RuntimeError(f" {self.__class__.__name__}:"
                               f" Water table depth observations contain NaN values.")
        # _end_if_

        # Since the observational data are not "gridded" we put them
        # on the z-grid [L: cm], before store them in the dictionary.
        self.mData["zWtd_cm"] = np.array([z_grid[z_grid >= k][0] for k in z_wtd_cm])

        # The precipitation data are already in [L: cm].
        precip_cm = np.array(data.loc[:, "Precipitation_cm"])

        # Check if there are NaN values.
        if np.any(np.isnan(precip_cm)):
            raise ValueError(f" {self.__class__.__name__}:"
                             f" Precipitation observations contain NaN values.")
        # _end_if_

        # Store to dictionary.
        self.mData["precipitation_cm"] = precip_cm

        # TRANSPIRATION AND HYDRAULIC LIFT PARAMETERS:

        # Compute the total atmospheric demand (i.e. root water uptake)
        # as percentage of the total precipitation over the whole year.
        atm_demand = params["Environmental"]["Atmospheric_Demand"]
        interception = params["Environmental"]["Interception_pct"]

        # Store to dictionary.
        self.mData["interception"] = interception

        # total_atm_demand_cm = np.sum(atm_demand*(1.0 - interception)*precip_cm)
        # This is a test case where we use the precipitation demand from the year
        # 2008-2009. The '13.4253' [cm] correspond to the 10% demand of that year.
        total_atm_demand_cm = 13.4253 * 10 * atm_demand

        # Make sure evapo-transpiration percentage "et_pct" is within bounds.
        # NOTE: Wet_Season_pct declares how much of the losses are due to ET.
        et_pct = np.minimum(np.maximum(params["Environmental"]["Wet_Season_pct"], 0.0), 1.0)

        # Get the amount of ET per season [L: cm].
        wet_et_cm = et_pct * total_atm_demand_cm
        dry_et_cm = (1.0 - et_pct) * total_atm_demand_cm

        # Dry counter.
        n_dry = 0

        for tk in timestamps:
            # Check if the month belongs in the "dry-months" list.
            if tk.month in [4, 5, 6, 7, 8, 9]:
                n_dry += 1
            # _end_if_
        # _end_for_

        # Get the number of time-points.
        dim_t = timestamps.size

        # Add it to the dictionary.
        self.mData["dim_t"] = dim_t

        # Wet counter (complementary to the Dry counter)
        n_wet = dim_t - n_dry

        # Preallocate the (plant) transpiration array.
        plant_tp = np.zeros(dim_t)

        # Assign each time-point the value for the ET.
        for i, tk in enumerate(timestamps):
            # Dry season ET:
            # Distribute the dry season amount evenly
            # in the time points that fall in that season.
            if tk.month in [4, 5, 6, 7, 8, 9]:
                plant_tp[i] = 2.0*dry_et_cm/n_dry
            # _end_if_

            # Wet season ET:
            if tk.month in [10, 11, 12, 1, 2, 3]:
                plant_tp[i] = 2.0*wet_et_cm/n_wet
            # _end_if_

            # N.B.: Here we multiply by 2 each value because half of the
            # points fall in the day and the other half during the night.
        # _end_for_

        # SAFEGUARD: Replace NaN with zero.
        plant_tp[np.isnan(plant_tp)] = 0.0

        # Hydraulic lift related parameters.
        theta_50 = np.maximum(0.5 ** (1.0 / K.lambda_exponent), 0.05)

        # Inverse of $\psi_{50}$: assumes m=0.5, n=2.
        inv_psi_50 = soil.alpha / np.sqrt(theta_50 ** (-2.0) - 1.0)

        # Store tree-root related parameters in this structure.
        self.mData["atm"] = plant_tp
        self.mData["LAI"] = params["Trees"]["Leaf_Area_Index"]
        self.mData["iPsi_50"] = inv_psi_50

        # EVAPORATION:
        # We assume a flat value for the surface evaporation throughout the whole
        # year. Since the discrete time-points are half during the night and half
        # during the day we multiply the values by 2.

        # Set the evaporation uniformly.
        self.mData["surface_evap"] = 2.0 * np.sum(params["Environmental"]["Evaporation_pct"] * precip_cm) / dim_t

        # Create a Richards' equation object.
        self.pde_model = RichardsPDE(self.mData)

        # Get the filename (or None).
        ic_data_file = params["IC_Filename"]

        # Check if there is an initial conditions file.
        if ic_data_file is None:
            # Create an initial conditions vector.
            self.mData["initial_cond"] = self.initial_conditions()
        else:
            # Open the initial conditions file.
            with open(Path(ic_data_file), 'r') as input_file:
                # The file should have four columns.
                init_cond = pd.read_csv(input_file, names=["IC"])

                # Make sure its numpy.array.
                init_cond = np.array(init_cond.loc[:, "IC"])

                # Check if the dimensions match.
                if z_grid.shape != init_cond.shape:
                    raise RuntimeError(f" {self.__class__.__name__}:"
                                       f" IC vector's dimensions do not match the spatial grid.")
                # _end_if_

                # Print a message.
                print(f" IC vector was loaded successfully from: {ic_data_file}.")
            # _end_with_

            # Create an initial conditions vector.
            self.mData["initial_cond"] = init_cond
        # _end_if_

    # _end_def_

    def initial_conditions(self):
        """
        Creates an initial conditions vector, by running the PDE
        model forward in time for 'burn-in' number of iterations.

        :return: y0 state vector [dim_d].
        """

        # [WARNING] Set the flag to TRUE!
        self.mData["sim_flags"]["SPINUP"] = True

        # Extract the well number.
        well_no = self.mData["Well_No"]

        # Display info.
        print(f" [Initial Conditions for Well no. {well_no}]"
              f" Burn in period started ...")

        # Spatial domain.
        z = self.mData["z_grid"]

        # Update the wtd parameter in the structure.
        wtd_i = np.where(z == self.mData["zWtd_cm"][0])

        # Get the porosity profile.
        q_0, *_ = self.mData["porosity"]()

        # Compute the pressure head values.
        y0, *_ = self.mData["hydro_model"].pressure_head(q_0, z)

        # Set a maximum number of iterations.
        burn_in = 1500

        # Early stop flag.
        early_stop = False

        # Initial random vector.
        n_rnd = self.rng.standard_normal(z.size)

        # Create a local dictionary with parameters
        # for the specific iteration.
        args_0 = {"wtd": wtd_i[0][0],
                  "n_rnd": n_rnd,
                  "atm": self.mData["atm"][0],
                  "time": self.mData["time"][0],
                  "interception": self.mData["interception"],
                  "precipitation": self.mData["precipitation_cm"][0]}

        # Local copy (reference).
        dz = self.mData["dz"]

        # Local copy (reference).
        psi_sat = self.mData["soil"].psi_sat

        # Burn-in loop.
        for j in range(burn_in):

            # During the burn-in period we do not update the
            # random field. Moreover, the actual times don't
            # matter. Time span can be fixed.
            t_span = (0, 1)

            # Solve the PDE.
            y_j = self.pde_model.solve(t_span, y0, args_0)

            # Find estimated wtd at the j-th iteration.
            wtd_est = find_wtd(y_j >= psi_sat)

            # Compute the absolute error: |wtd_obs - wtd_est|
            abs_error = np.abs(self.mData["zWtd_cm"][0] - z[wtd_est])

            # Find the MSE between the two state vector solutions.
            mse_0 = np.mean((y_j - y0) ** 2)

            # Set the vector for the next iteration.
            y0 = y_j.copy()

            # Repeat the burn-in integration as long as the distance
            # between the two solutions is above a threshold value.
            if abs_error <= (2.0 * dz) and (mse_0 <= 0.01):
                # Change the flag.
                early_stop = True

                # Display final message.
                print(f" [Initial Conditions for Well no. {well_no}]"
                      f" finished at [itr: {j}] with [abs(error): {abs_error}]"
                      f" and [MSE: {mse_0}]")

                # Exit the loop.
                break
            # _end_if_
        # _end_for_

        # At this point the algorithm has reached maximum number of iterations.
        if not early_stop:
            print(f" [Initial Conditions for Well no. {well_no}]"
                  f" finished at maximum number of iterations.")
        # _end_of_

        # [WARNING] Set the flag to TRUE!
        self.mData["sim_flags"]["SPINUP"] = False

        # Return the state vector.
        return y0
    # _end_if_

    def run(self):
        """
        Runs the PDE model forward in time. All the output information
        is store in self.output dictionary, where we can later save it
        to the disk for further analysis.

        :return: None
        """
        # Check if the model parameters have been initialized.
        if not self.mData:
            raise RuntimeError(f" {self.__class__.__name__}:"
                               f" Simulation data structure 'mData' is empty.")
        # _end_if_

        # Transpiration and lateral flow values.
        # Reset to empty lists before each run.
        transpiration, lateral_flow = [], []

        # Get the initial conditions vector.
        y0 = self.mData["initial_cond"].copy()

        # Extract pressure head at saturation.
        psi_sat = self.mData["soil"].psi_sat

        # Extract hydrological model.
        h_model = self.mData["hydro_model"]

        # Spatial domain.
        z = self.mData["z_grid"]

        # Space discretization size.
        dim_d = z.size

        # Time domain size.
        dim_t = self.mData["dim_t"]

        # Volumetric water content.
        theta_vol = np.zeros((dim_t, dim_d))

        # Pressure head (suction).
        psi = np.zeros((dim_t, dim_d))

        # Background hydraulic conductivity.
        # Note: This is the equivalent of the
        # saturated in the vanGenuchten model.
        k_bkg = np.zeros((dim_t, dim_d))

        # Unsaturated hydraulic conductivity.
        k_hrc = np.zeros((dim_t, dim_d))

        # Water table depths (estimated).
        wtd_est = np.zeros(dim_t, dtype=int)

        # Error (absolute) of water table depth.
        abs_error = np.zeros(dim_t)

        # Initial wtd value.
        wtd_est[0] = find_wtd(y0 >= psi_sat)

        # Update the error (at t=t0).
        abs_error[0] = np.abs(self.mData["zWtd_cm"][0] - z[wtd_est[0]])

        # Save the initial conditions at time $t = 0$.
        psi[0] = y0.copy()

        # Create a random vector.
        n_rnd = self.rng.standard_normal(dim_d)

        # Run the hydrological model to get the initial water content.
        theta_vol[0], k_hrc[0], _, k_bkg[0], *_ = h_model(y0, z, {"n_rnd": n_rnd})

        # Virtual time array.
        tk = np.arange(0, dim_t)

        # Create a Richards' equation object.
        pde_model = self.pde_model

        # Start the timer.
        time_t0 = time.time()

        # Integrate the 'PDE' in tk. Each iteration corresponds to 0.5hr interval.
        for i in range(1, tk.size):

            # Find the index (location) of the target water table depth.
            wtd_i = np.where(z == self.mData["zWtd_cm"][i])

            # Check if the index has been found.
            if wtd_i[0].size == 0:
                # Print a warning message.
                print(" Warning: Observation {0} cannot be found in the spatial domain."
                      " Skipping iteration {1} ...".format(self.mData["zWtd_cm"][i], i))
                # SKip.
                continue
            # _end_if_

            # Create a local dictionary with parameters for the specific iteration.
            args_i = {"wtd": wtd_i[0][0],
                      "n_rnd": n_rnd,
                      "atm": self.mData["atm"][i],
                      "time": self.mData["time"][i],
                      "interception": self.mData["interception"],
                      "precipitation": self.mData["precipitation_cm"][i]}

            # Update the random field (at least) daily (~24hr):
            if (args_i["precipitation"] > 0.5) or (np.mod(i, 48) == 0):
                # Standard normal random variables ~ N(0,1):
                args_i["n_rnd"] = self.rng.standard_normal(dim_d)
            # _end_if_

            # This time-window corresponds to '30' minutes time intervals
            # between observations.
            t_span = (tk[i-1], tk[i])

            # Solve the PDE and return the solution at the final time point.
            y_i = pde_model.solve(t_span, y0, args_i)

            # Find $wtd$ at the k-th iteration.
            wtd_est[i] = find_wtd(y_i >= psi_sat)

            # Compute the absolute error: |wtd_{obs} - wtd_{est}|
            abs_error[i] = np.abs(self.mData["zWtd_cm"][i] - z[wtd_est[i]])

            # Store the pressure head to the array.
            psi[i] = y_i.copy()

            # Recover the values of the volumetric water content $\theta$
            # and the hydraulic conductivities $K(\theta)$ and 'Kbkg' for
            # the values of the final solution of the PDE.
            theta_vol[i], k_hrc[i], _, k_bkg[i], *_ = h_model(y_i, z, args_i)

            # Update the initial condition vector for the next iteration.
            y0 = y_i.copy()

            # Store cumulative output from the PDE.
            lateral_flow.append(pde_model.arg_out["lateral_flow"])
            transpiration.append(pde_model.arg_out["transpiration"])

            # Display summary statistics (every ~100 iterations).
            if np.mod(i, 100) == 0:
                # Below the water table.
                w_side = '+'

                # Above the water table.
                if self.mData["zWtd_cm"][i] < z[wtd_est[i]]:
                    w_side = '-'
                # _end_if_

                # Compute the current Mean Absolute Error.
                mae = np.mean(abs_error[0:i])

                # Display message.
                print(" [Well No. {0}] {1}: MAE = {2:.2f} cm,"
                      " [{3}]".format(self.mData["Well_No"], i, mae, w_side))
            # _end_if_
        # _end_for_

        # Stop the timer.
        time_tf = time.time()

        # Print duration.
        print(" Elapsed time: {0:.2f} seconds.\n".format(time_tf - time_t0))

        # Compute the effective saturation (normalized water content).
        s_eff = theta_vol / self.mData["porosity"]()[0]

        # Prepare the output.
        self.output["K_hrc"] = k_hrc
        self.output["K_bkg"] = k_bkg
        self.output["S_eff"] = s_eff
        self.output["psi_press"] = psi
        self.output["theta_vol"] = theta_vol
        self.output["abs_error"] = abs_error
        self.output["wtd_est_cm"] = z[wtd_est]

        # Convert lists to arrays, before adding them.
        self.output["lateral_flow"] = np.array(lateral_flow)
        self.output["transpiration"] = np.array(transpiration)
    # _end_def_

    def saveResults(self):
        """
        Saves the simulation results to a file. All the data should be stored
        inside the self.output dictionary and be of type "numpy.ndarray". For
        the moment the file is saved in the same directory as the main program.

        :return: None.
        """

        # Check if the output dictionary is empty.
        if not self.output:
            # Print a message and do not save anything.
            print(f" {self.__class__.__name__}:"
                  f" Simulation data structure 'output' is empty.")
        else:
            # Initial message.
            print(f" Saving the results to: {self.name}.h5")

            # Create the output filename. Remove spaces (if any).
            file_out = Path(self.name.strip().replace(" ", "_") + ".h5")

            # Save the data to an 'HDF5' file format.
            # NOTE:  Create file; truncate if exists.
            with h5py.File(file_out, 'w') as out_file:
                # Local reference.
                data = self.output

                # Extract all the data.
                for key in data:
                    # Default compressions level is '4'.
                    out_file.create_dataset(key, data=data[key],
                                            shape=data[key].shape,
                                            compression='gzip')
                # _end_for_
            # _end_with_

        # _end_if_
    # _end_def_

# _end_class_

# Auxiliary function.
def loadResults(filename=None):
    """
    Loads the simulation data that have been stored in the hdf5 file
    using the object's method saveResults().

    :param filename: (string) is the '.h5' file that contains the data.

    :return: a dictionary with all the data.
    """

    # Check if we have given input file.
    if filename is None:
        raise RuntimeError(" load_data: No input file is given.")
    # _end_if_

    # Simulation data.
    sim_data = {}

    # Open the file for read only.
    with h5py.File(Path(filename), 'r') as input_file:

        # Extract all the data.
        for key in input_file:
            sim_data[key] = np.array(input_file[key])
        # _end_for_

    # _end_with_

    # Return the dictionary.
    return sim_data
# _end_def_
