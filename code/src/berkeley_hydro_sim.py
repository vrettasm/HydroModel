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


# Main function.
def main(model_parameters=None, simulation_data=None):
    # Check if we got model parameters.
    if model_parameters:
        print(" Model parameters are given.")
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
