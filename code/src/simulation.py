import json
import numpy as np
import pandas as pd
from pathlib import Path

from code.src.porosity import Porosity
from code.src.water_content import WaterContent
from code.src.soil_properties import SoilProperties
from code.src.hydraulic_conductivity import HydraulicConductivity

class Simulation(object):
    """
    TBD
    """

    def __init__(self, name=None):
        """
        Default constructor of the Simulation class.

        :param name: (string) is optional but it will be used for constructing
        a meaningful filename to save the results at the end of the simulation.
        """
        # Check if a simulation name has been given.
        if not name:
            self.name = "ID_None"
        else:
            self.name = name
        # _end_if_

        # This dictionary will carry all the simulation data.
        self.mData = {}
    # _end_def_

    def setupModel(self, params=None, data=None):
        """
        This method is called BEFORE the run() and sets up all the variables for the simulation.
        It is also responsible for checking the validity of the input parameters before use.

        :param params: (dict) contains all the given parameters. If None, then it
        will use the default parameters to initialize object. However this is not
        recommended since the user will have no control on the parameters.

        :param data: (pandas.DataFrame) that will contain the water data for the
        simulation. These include, time-series of the dates, water table depths
        precipitation values, etc.

        :return: None.
        """

        # make sure we have input parameters.
        if not params:
            raise RuntimeError(" Can't setup the model parameters.")
        # _end_if_

        # Make sure we have water data.
        if not data:
            raise RuntimeError(" Can't setup the model data.")
        # _end_if_

        # Copy the well number that is being simulated.
        self.mData["Well_No"] = params["Well_No"]

        # Open the file in "Read Only" mode.
        with open(Path(params["Site_Information"]), 'r') as site_file:
            # Load the site information.
            site_info = json.load(site_file)
        # _end_with_

        # Check if the Well number exists in the site info file.
        if not str(self.mData["Well_No"]) in site_info["Well"]:
            raise ValueError(" The selected well does not exist in the site information file.")
        # _end_if_

        # Extract the Well information.
        well = site_info["Well"][str(self.mData["Well_No"])]

        # Make sure the well is not defined as fully saturated.
        if well["sat_depth"] >= well["max_depth"]:
            raise ValueError(" The well seems fully saturated.")
        # _end_if_

        # Spatial domain parameters for the underground layers [L:cm]:
        layers = (well["soil"], well["saprolite"], well["weathered"], well["max_depth"])

        # Add the underground layers to the structure.
        self.mData["layers"] = layers

        # The spacing here defines a Uniform-Grid [L:cm].
        dz = 5.0

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
            print(" SoilProperties failed to initialize: {}.\n"
                  " It will use default initialization parameters.".format(e0))

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
            print(" WaterContent failed to initialize: {}.\n"
                  " It will use default initialization parameters.".format(e0))

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
            print(" HydraulicConductivity failed to initialize: {}.\n"
                  " It will use default initialization parameters.".format(e0))

            # Default initialization.
            K = HydraulicConductivity()
        # _end_try_

        # Add hydraulic conductivity to the dictionary.
        self.mData["K"] = K

        # Add the environmental parameters.
        self.mData["env_param"] = params["Environmental"]

        # Add the selected hydrological model.
        self.mData["hydro_model"] = params["Hydrological_Model"]

        # Add the simulation (execution) flags.
        self.mData["sim_flags"] = params["Simulation_Flags"]

        # Create and add the porosity object.
        self.mData["porosity"] = Porosity(z_grid, layers, theta, soil,
                                          params["Hydrological_Model"]["Porosity_Profile"])

        # Extract the observational data from the pandas.Dataframe:
        r_datenum = data.loc[:, "Datenum"]

        # We need to convert the dates (from MATLAB to Python).
        # NOTE: The value 719529 is MATLAB's datenum value of the "Unix epoch"
        # start (1970-01-01), which is the default origin for pd.to_datetime().
        timestamps = pd.to_datetime(r_datenum-719529, unit='D')

        # Get the number of time-points.
        dim_t = timestamps.size

        # Store timestamps in the dictionary.
        self.mData["time"] = [t.round(freq="s") for t in timestamps]

        # Convert the meters to [L:cm] before storing them.
        z_wtd_cm = np.array(np.abs(np.round(100.0*data.loc[:, "WTD_m"])))

        # Check if there are NaN values.
        if np.any(np.isnan(z_wtd_cm)):
            raise ValueError(" Water table depth observations contain NaN.")
        # _end_if_

        # Since the observational data are not "gridded" we put them
        # on the z-grid [L: cm], before store them in the dictionary.
        self.mData["zWtd_cm"] = np.array([z_grid[z_grid >= k][0] for k in z_wtd_cm])

        # The precipitation data are already in [L: cm].
        precip_cm = np.array(data.loc[:, "Precipitation_cm"])

        # Check if there are NaN values.
        if np.any(np.isnan(precip_cm)):
            raise ValueError(" Precipitation observations contain NaN.")
        # _end_if_

        # Store to dictionary.
        self.mData["precip_cm"] = precip_cm

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
        wet_et_cm = et_pct*total_atm_demand_cm
        dry_et_cm = (1.0 - et_pct)*total_atm_demand_cm

        # Dry counter.
        n_dry = 0

        for tk in timestamps:
            # Check if the month belongs in the "dry-months" list.
            if tk.month in [4, 5, 6, 7, 8, 9]:
                n_dry += 1
            # _end_if_
        # _end_for_

        # Wet counter (complementary to the Dry counter)
        n_wet = dim_t-n_dry

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
        self.mData["iPsi_50"] = inv_psi_50

        # EVAPORATION:
        # We assume a flat value for the surface evaporation throughout the whole
        # year. Since the discrete time-points are half during the night and half
        # during the day we multiply the values by 2.

        # Set the evaporation uniformly.
        self.mData["surface_evap"] = 2.0 * np.sum(params["Environmental"]["Evaporation_pct"] * precip_cm) / dim_t

    # _end_def_

    def run(self):
        # Check if the model parameters have been initialized.
        if not self.mData:
            print(" Simulation data structure ('mData') is empty.")
        # _end_if_
    # _end_def_

    def saveResults(self):
        pass
    # _end_def_

# _end_class_
