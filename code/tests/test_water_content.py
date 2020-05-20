import unittest
from code.src.water_content import WaterContent

class TestWaterContent(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        print(" >> TestWaterContent - START -")
    # _end_def_

    @classmethod
    def tearDownClass(cls) -> None:
        print(" >> TestWaterContent - FINISH -")
    # _end_def_

    def setUp(self) -> None:
        """
        Creates a test object with the default input parameters.
        :return: None.
        """
        self.test_obj = WaterContent()
    # _end_def_

    def tearDown(self) -> None:
        """
        Check the final status of the object before destruction.
        :return: None
        """
        self.assertTrue(self.test_obj._checkValues())
    # _end_def_

    def test_max(self):
        """
        Test the 'max' accessor of the object. A value more
        than one (1.0) should raise a ValueError.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.max += 1.0
        # _end_with_
    # _end_def_

    def test_min(self):
        """
        Test the 'min' accessor of the object. A value less
        than zero (0.0) should raise a ValueError.
        :return: None
        """
        with self.assertRaises(ValueError):
            self.test_obj.min -= 1.0
        # _end_with_
    # _end_def_

    def test_res(self):
        """
        Test the 'res' accessor of the object. Residual water content
        must be less then the minimum value (but larger then zero).
        :return: None
        """
        with self.assertRaises(ValueError):
            # Setting the residual value '+10%' higher then
            # the minimum should raise an error (ValueError).
            self.test_obj.res = 1.1*self.test_obj.min
        # _end_with_
    # _end_def_

    def test_mid(self):
        """
        Test the 'mid' accessor of the object.
        :return: None
        """
        # Get the mid value manually.
        tmp_mid = (self.test_obj.min + self.test_obj.max)/2.0

        # Mid operator should return the same value.
        self.assertEqual(tmp_mid, self.test_obj.mid)
    # _end_def_

    def test_wrong_init_params(self):
        """
        Test an object initialization with wrong input parameters.
        :return: None
        """
        with self.assertRaises(ValueError):
            _ = WaterContent(minimum=0.2, maximum=0.1)
        # _end_with_
    # _end_def_


if __name__ == '__main__':
    unittest.main()
