# Import packages here.

class Simulation(object):
    """

    """

    def __init__(self, name=None):
        # Check if a simulation name has been given.
        if not name:
            self.name = "ID:None"
        else:
            self.name = name
        # _end_if_

        # This dictionary will carry all the simulation data.
        self.mData = {}
    # _end_def_

    def setupModelParameters(self, params=None, data=None):
        pass
    # _end_def_

    def run(self):
        # Check if the model parameters have been initialized.
        if not self.mData:
            print(" Simulation mData is empty.")
        # _end_if_

    # _end_def_

    def saveResults(self):
        pass
    # _end_def_

# _end_class_
