"""Tests srim reader"""

import nose 
from nose.tools import assert_equal, assert_true, assert_almost_equal, assert_raises
from pyne import srim_read as sr


sr_read1 = sr.read_srim_output('Helium')

def test_read_srim_output():
    assert_equal(sr_read1.bragg_corr,-0.83)
    assert_equal(sr_read1.sr_read1.filename,'SRIM Outputs\\Helium in Mylar (ICRU-222)')
    assert_equal(sr_read1.ion,'Helium')
    assert_equal(sr_read1.sr_read1.ion_mass,4.003)
    assert_equal(sr_read1.date,'November 15, 2017')
    assert_equal(sr_read1.stopping_units,' MeV / (mg/cm2) ')
    assert_equal(sr_read1.target_density_grams,'1.3970E+00 g/cm3')
    assert_equal(sr_read1.target_density_atoms,'9.6311E+22 atoms/cm3')

