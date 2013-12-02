""" spectrometry tests """
import nose 
from nose.tools import assert_equal

import gammaSpec
import spectAnalysis

gspec1=gammaSpec.readspefile("test.spe")

# gspec1.print_parameters()
# spectAnalysis.plot_spectra(gspec1)
# spectAnalysis.plot_peak(gspec1, 843)


def test_read_times():
    assert_equal(gspec1.real_time, 300.0)
    assert_equal(gspec1.live_time, 274.0)
    assert_equal(gspec1.start_time, "12:18:05  ")
    assert_equal(gspec1.start_date, "15-Mar-2012")
    assert_equal(gspec1.dead_time, 26.0)
    
    
def test_read_det():
    assert_equal(gspec1.det_ID, "1")
    assert_equal(gspec1.det_descp, "DSPEC PRO 1")

if __name__ == "__main__":
    nose.runmodule()
