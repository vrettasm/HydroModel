from numpy import maximum as np_max
from numpy import minimum as np_min


class WaterContent(object):
    """
    This class represents Water content, or soil moisture content.
    """

    __slots__ = ("minimum", "maximum", "residual", "wilting", "field_cap")

    # Default constructor.
    def __init__(self, minimum: float = 0.08, maximum: float = 0.30,
                 residual: float = 0.05, wilting: float = -1500.0,
                 field_cap: float = 340.0):
        """
        Constructs an object that will hold all the water content variables.

        :param minimum: (float) minimum water content (dimensionless [-]).

        :param maximum: (float) maximum water content (dimensionless [-]).

        :param residual: (float) residual water content, defined as the water
        content for which the gradient $dTheta/dh$ becomes zero, and $theta_{sat}$
        is the saturated water content, which is equivalent to porosity, $phi$
        (dimensionless [-]).

        :param wilting: (float) permanent wilting point or wilting point is
        defined as the minimum amount of water in the soil that the plant
        requires not to wilt.

        :param field_cap: (float) field capacity is the amount of soil moisture
        or water content held in the soil after excess water has drained away
        and the rate of downward movement has decreased.

        :raises ValueError: if there checkValues method fails.
        """

        # Do not allow these values to be bellow zero, or above one.
        self.minimum = WaterContent.inRange(minimum)
        self.maximum = WaterContent.inRange(maximum)
        self.residual = WaterContent.inRange(residual)

        # Check if the input values are in the right order.
        if not self._checkValues():
            raise ValueError(f" {self.__class__.__name__}:"
                             f" The volumetric water content input values are incorrect.")
        # _end_if_

        # These are in length units [L:cm] and can have negative values.
        self.wilting = wilting
        self.field_cap = field_cap
    # _end_def_

    @staticmethod
    def inRange(new_value: float):
        """
        Sets the new input value within the [0-1] range. This ensures
        that no negative values will enter the object.

        :param new_value: (float) represents a water content value
        (dimensionless [-]).

        :return: a new value in [0-1].
        """
        return np_max(np_min(1.0, new_value), 0.0)
    # _end_def_

    def _checkValues(self):
        """
        Check whether the water content values satisfy specific conditions.

        :return: True if the condition is satisfied.
        """
        return 0.0 <= self.residual < self.minimum < self.maximum <= 1.0
    # _end_def_

    @property
    def min(self):
        """
        Accessor method.

        :return: the minimum water content.
        """
        return self.minimum
    # _end_def_

    @min.setter
    def min(self, new_value):
        """
        Accessor method.

        :param new_value: for the minimum water content.

        :return: None.
        """
        # Keep the old minimum value.
        old_value = self.minimum

        # Make the change.
        self.minimum = new_value

        # Check if the change has broken the correct order of the values.
        if not self._checkValues():
            # Set back the old value.
            self.minimum = old_value

            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}: The new minimum: "
                             f" {new_value}, is not consistent with the rest"
                             f" of the values.")
        # _end_if_
    # _end_def_

    @property
    def max(self):
        """
        Accessor method.

        :return: the maximum water content.
        """
        return self.maximum
    # _end_def_

    @max.setter
    def max(self, new_value):
        """
        Accessor method.

        :param new_value: for the maximum water content.

        :return: None.
        """
        # Keep the old minimum value.
        old_value = self.maximum

        # Make the change.
        self.maximum = new_value

        # Check if the change has broken the correct order of the values.
        if not self._checkValues():
            # Set back the old value.
            self.maximum = old_value

            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}: The new maximum: {new_value}"
                             f" is not consistent with the rest of the values.")
        # _end_if_
    # _end_def_

    @property
    def mid(self):
        """
        Auxiliary water content function.

        :return: the mid-point value between the min/max values.
        """
        return 0.5*(self.maximum + self.minimum)
    # _end_def_

    @property
    def res(self):
        """
        Accessor method.

        :return: the residual water content.
        """
        return self.residual
    # _end_def_

    @res.setter
    def res(self, new_value):
        """
        Accessor method.

        :param new_value: for the residual water content.

        :return: None.
        """
        # Keep the old minimum value.
        old_value = self.residual

        # Make the change.
        self.residual = new_value

        # Check if the change has broken the correct order of the values.
        if not self._checkValues():
            # Set back the old value.
            self.residual = old_value

            # Raise an error with a message.
            raise ValueError(f" {self.__class__.__name__}: The new residual: {new_value}"
                             f" is not consistent with the rest of the values.")
        # _end_if_
    # _end_def_

    @property
    def wlt(self):
        """
        Accessor method.

        :return: the wilting point.
        """
        return self.wilting
    # _end_def_

    @wlt.setter
    def wlt(self, new_value):
        """
        Accessor method.

        :param new_value: for the wilting point.

        :return: None.
        """
        self.wilting = new_value
    # _end_def_

    @property
    def flc(self):
        """
        Accessor method.

        :return: the field capacity.
        """
        return self.field_cap
    # _end_def_

    @flc.setter
    def flc(self, new_value):
        """
        Accessor method.

        :param new_value: for the field capacity.

        :return: None.
        """
        self.field_cap = new_value
    # _end_def_

    # Auxiliary.
    def __str__(self):
        """
        Override to print a readable string presentation of the object.
        This will include its id(), along with its fields values.

        :return: a string representation of a WaterContent object.
        """
        # New line character.
        new_line = '\n'

        return f" WaterContent Id({id(self)}): {new_line}"\
               f" Minimum={self.minimum}, Maximum={self.maximum} {new_line}"\
               f" Residual={self.residual}, Wilting={self.wilting} {new_line}"\
               f" Field capacity={self.field_cap}"
    # _end_def_

# _end_class_
