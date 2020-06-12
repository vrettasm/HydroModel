# Import locally the json package.
import json
import numpy as np
from pathlib import Path
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
            self.name = "ID:None"
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

        # Create a vertical grid (increasing downwards)
        z_grid = np.arange(well["soil"], well["max_depth"]+dz, dz)

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
        self.mData["Env"] = params["Environmental"]

        # Add the selected hydrological model.
        self.mData["hModel"] = params["Hydrological_Model"]

        # Add the simulation (execution) flags.
        self.mData["sim_flags"] = params["Simulation_Flags"]

        # Add the water (well/precipitation/etc) data.
        self.mData["wy_data"] = data
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
