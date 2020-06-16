#!/usr/bin/env python3

"""
Berkeley Hydrological Simulation

TBD
"""

# Load the required packages.
import sys
import pandas as pd
from pathlib import Path
from code.src.simulation import Simulation

# INFO:
__author__ = "Michail Vrettas, PhD"
__email__ = "michail.vrettas@gmail.com"


def validateInputParametersFile(filename):
    """
    Validates an input (json) file to check if it contains
    the required keys. It does not validate the values of
    the keys.

    :param filename: Is a "Path" object that contains the input
    model parameters for the simulation.

    :return: A dictionary loaded from the input file.

    :raises ValueError if a keyword is missing from the file.
    """

    # Import locally the json package.
    import json

    # Open the file in "Read Only" mode.
    with open(filename, 'r') as input_file:

        # Load the model parameters.
        model_params = json.load(input_file)

        # Required keys in the json file.
        required_keys = ["IC_Filename", "Data_Filename", "Well_No",
                         "Water_Content", "Hydraulic_Conductivity",
                         "Hydrological_Model",  "Simulation_Flags",
                         "Environmental", "Soil_Properties", "Site_Information"]

        # Check the keywords for membership in the file.
        for k in required_keys:

            # The order in here doesn't matter.
            if k not in model_params:
                raise ValueError(" Key: {0}, is not given.".format(k))
            # _end_if_

        # _end_for_

        # Show message.
        print(" Model parameters are given correctly.")
    # _end_with_

    # The dictionary will contain all the input parameters.
    return model_params
# _end_def_


# Main function.
def main(params_file=None, data_file=None):
    """
    As the name suggests, this is the main function that is called to initiate the simulation run.

    :param params_file: (string) that points to the input file for the parameters.

    :param data_file: (string) that points to the input file for the water data.

    :return: None
    """

    # Check if we got model parameters.
    if params_file:
        try:
            # Make sure params_file is a Path object.
            params_file = Path(params_file)

            # Check if everything is ok.
            params = validateInputParametersFile(params_file)
        except ValueError as e0:
            # Show the error message.
            print(e0)

            # Exit the program.
            sys.exit(1)
        # _end_try_
    else:
        print(" The simulation can't run without input parameters.")

        # Exit the program.
        sys.exit(1)
    # _end_if_

    # Check if we got simulation water data. Make sure its a Path object.
    if data_file:
        data_file = Path(data_file)
    else:
        data_file = Path(params["Data_Filename"])
    # _end_if_

    # Display where we got the water data from.
    print(" Simulation water data file: {0}".format(data_file))

    try:
        # Open the water data in "Read Only" mode.
        with open(data_file, 'r') as input_file:
            # The file should have four columns.
            water_data = pd.read_csv(input_file,
                                     names=["ID", "Datenum", "Precipitation_cm", "WTD_m"])
        # _end_with_

        # Create a simulation object.
        sim_01 = Simulation("Sim_01")

        # Setup its parameters and water data.
        sim_01.setupModel(params, water_data)

        # Run the simulation.
        sim_01.run()

        # Save the results.
        sim_01.saveResults()
    except Exception as e1:
        print(e1)

        # Exit the program.
        sys.exit(1)
    # _end_try_

# _end_maim_


# Run the script.
if __name__ == "__main__":

    # Check if we have given input parameters.
    if len(sys.argv) > 1:
        # Local import.
        import argparse

        # Create a parser object
        parser = argparse.ArgumentParser(description=" Berkeley Hydrological Simulation ")

        # Input file with simulation parameters.
        parser.add_argument("--params",
                            help=" Input file (.json) with simulation parameters.")

        # Input file with simulation data.
        parser.add_argument("--data",
                            help=" Input file (.csv) with simulation data (e.g.: precipitation, wtd).")

        # Parse the arguments.
        args = parser.parse_args()

        # Call the main function.
        main(args.params, args.data)

        # Display final info.
        print(' Simulation completed.')
    else:
        sys.exit('Error: Not enough input parameters.')
    # _end_if_

# _end_program_
