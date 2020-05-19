from numpy import minimum as np_min
from numpy import maximum as np_max

class WaterContent(object):
    """

    """

    # TBD
    __slots__ = ["minimum", "maximum", "residual", "wilting", "field_cap"]

    # Default constructor.
    def __init__(self, minimum: float, maximum: float, residual: float,
                 wilting: float, field_cap: float):
        """
        Constructs an object that will hold all the water content variables.

        :param minimum: (float) minimum water content.

        :param maximum: (float) maximum water content.

        :param residual: (float) residual water content, defined as the water content for
        which the gradient $dTheta/dh$ becomes zero, and $theta_{sat}$ is the saturated
        water content, which is equivalent to porosity, $phi$.

        :param wilting: (float) permanent wilting point or wilting point is defined as the
        minimum amount of water in the soil that the plant requires not to wilt.

        :param field_cap: (float) field capacity is the amount of soil moisture or water
        content held in the soil after excess water has drained away and the rate
        of downward movement has decreased.
        """

        # Do not allow these values to be bellow zero.
        self.minimum = self.inRange(minimum)

        self.maximum = self.inRange(maximum)

        self.residual = self.inRange(residual)

        # These are in length units [L:cm] and can have negative values.
        self.wilting = wilting

        self.field_cap = field_cap
    # _end_def_

    @staticmethod
    def inRange(new_value: float):
        """

        :param new_value:

        :return:
        """
        return np_max(np_min(1.0, new_value), 0.0)
    # _end_def_

    @property
    def min(self):
        return self.minimum
    # _end_def_

    @min.setter
    def min(self, new_value):
        self.minimum = self.inRange(new_value)
    # _end_def_

    @property
    def max(self):
        return self.maximum
    # _end_def_

    @max.setter
    def max(self, new_value):
        self.maximum = self.inRange(new_value)
    # _end_def_

    @property
    def mid(self):
        return 0.5*(self.maximum + self.minimum)
    # _end_def_

    @property
    def res(self):
        return self.residual
    # _end_def_

    @property
    def wlt(self):
        return self.wilting
    # _end_def_

    @property
    def flc(self):
        return self.field_cap
    # _end_def_
