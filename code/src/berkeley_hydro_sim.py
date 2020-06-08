#!/usr/bin/env python3

"""
Berkeley Hydrological Simulation

TBD
"""

# Load the required packages.
import os
import sys
import pandas as pd
from pathlib import Path

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

    :return: Nothing.

    :raises ValueError if a keyword is missing from the file.
    """

    # Import locally the json package.
    import json

    # Open the file in "Read Only" mode.
    with open(filename, 'r') as input_file:

        # Load tha model parameters.
        model_params = json.load(input_file)

        # Required keys in the json file.
        required_keys = ["IC_filename", "Data_filename", "Well_no",
                         "Water_Content", "Hydraulic_Conductivity",
                         "Hydrological_Model", "Simulation_Flags",
                         "Environmental"]

        # Check the keywords for membership in the file.
        for k in required_keys:

            # The order in here doesn't matter.
            if k not in model_params:
                raise ValueError(" Key: {0}, is not given.".format(k))
            # _end_if_

        # _end_for_
    # _end_with_

# _end_def_


# Main function.
def main(model_parameters=None, simulation_data=None):
    # Check if we got model parameters.
    if model_parameters:

        # Make sure model_parameters is a Path object.
        model_parameters = Path(model_parameters)

        try:
            # Check if everything is ok.
            validateInputParametersFile(model_parameters)

            # Show message.
            print(" Model parameters are given correctly.")
        except ValueError as e0:
            # Show the error message.
            print(e0)

            # Exit the program.
            sys.exit(1)
        # _end_try_

    else:
        if simulation_data:
            print(" The simulation will run with default parameters.")
        else:
            print(" The simulation can't run without datafile.")

            # Exit the program.
            sys.exit(1)
    # _end_if_

    # Check if we got simulation data.
    if simulation_data:
        print(" Simulation data from: {0}".format(simulation_data))
    else:
        print(" No data.")
    # _end_if_

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
