# Import locally the json package.
import json
import numpy as np
from pathlib import Path
from code.src.soil_properties import SoilProperties

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

        # Length of the vertical grid.
        dim_z = z_grid.size

        # Create a soil properties object.
        soil = SoilProperties(params["Soil_Properties"]["n"],
                              params["Soil_Properties"]["a0"],
                              params["Soil_Properties"]["psi_sat"],
                              params["Soil_Properties"]["epsilon"])
        return dim_z
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
