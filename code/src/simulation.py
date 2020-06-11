# Import packages here.

class Simulation(object):
    """
    TBD
    """

    def __init__(self, name=None):
        """
        Default constructor of the Simulation class.

        :param name: (string) is optional but it will be used for constructing
        meaningful a filename to save te results at the end of the simulation.
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
        pass
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
