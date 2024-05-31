import unittest

from brdr.utils import get_breakpoints_zerostreak


class TestUtils(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_series = series = [0,1,2,3,4,5,6,7,8,9,10]

    def test_get_breakpoints_zerostreak(self):
        sample_diff = [0, 5, 5, 3, 2, 2, 2, 6, 7, 7, 6]
        breakpoints, zerostreaks = get_breakpoints_zerostreak(self.sample_series,sample_diff)
        assert len(breakpoints) != 0
        assert len(zerostreaks) != 0

    def test_get_breakpoints_zerostreak_no_zerostreaks(self):
        breakpoints, zerostreaks = get_breakpoints_zerostreak(self.sample_series,self.sample_series)
        assert len(breakpoints) != 0
        assert len(zerostreaks) == 0
