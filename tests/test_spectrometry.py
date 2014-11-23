"""Spectrometry tests """
import nose 
from nose.tools import assert_equal, assert_true
import warnings

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import gammaspec

gspec1 = gammaspec.read_spe_file("test.spe")
eff_coeff = [-2.818615042612040000, -0.727352820018942000, -0.039579888648190400,
             -0.059230525466409600, 0.023772637347443000, 0.032530647507267100]


def test_read_spe():
    assert_equal(gspec1.real_time, 300.0)
    assert_equal(gspec1.live_time, 274.0)
    assert_equal(gspec1.start_time, "12:18:05  ")
    assert_equal(gspec1.start_date, "15-Mar-2012")
    assert_equal(gspec1.dead_time, 26.0)
    assert_equal(gspec1.det_id, "1")
    assert_equal(gspec1.det_descp, "DSPEC PRO 1")
    assert_equal(gspec1.start_chan_num, 0)
    assert_equal(gspec1.num_channels, 16384)
    
def test_calib():
    assert_equal(gammaspec.calc_e_eff(1, eff_coeff, 1), 0.059688551591347033)

def test_str():
    s = str(gspec1)
    assert_true(len(s) > 0)

if __name__ == "__main__":
    nose.runmodule()
