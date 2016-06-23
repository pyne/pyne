"""fispact tests """
import nose 
from nose.tools import assert_equal, assert_true, assert_almost_equal, assert_raises

import warnings
from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne import fispact

fispactii_path="fispii.out"

fo=fispact.read_fis_out(fispactii_path)

def test_read_fis_out():
    """test simple single values have been read """
    assert_equal(fo.file_name, fispactii_path)
    assert_equal(fo.version, "FISPACT-II")
    assert_true(fo.isFisII)
    assert_equal(fo.tot_irrad_time, 8.640000E+04)
    assert_equal(fo.tot_fluence, 8.640000E+10)
    assert_equal(fo.ave_flux, 1.000000E+06)
    assert_equal(fo.cpu_time, 14.612)
    assert_equal(fo.num_irrad_step, 1)

def test_read_time_step():
    """test reading time steps for fispact-II"""
    ts1 = fo.timestep_data[1]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]
    assert_equal(len(fo.timestep_data), 137)
    assert_equal(ts1.step_length, 8.6400E+04)
    assert_equal(ts2.step_length, 3.6000E+03)
    assert_equal(ts3.step_length, 2.2918E+07)
    assert_equal(ts1.alpha_act, 3.643021E+03)
    assert_equal(ts1.beta_act, 8.465725E+06)
    assert_equal(ts1.gamma_act, 1.290820E+06)
    assert_equal(ts1.total_act, 9.75930E+06)
    assert_equal(ts1.total_act_no_trit, 9.75928E+06)
    assert_equal(ts1.alpha_heat, 1.65855E-11)
    assert_equal(ts1.beta_heat, 7.88824E-10)
    assert_equal(ts1.gamma_heat, 1.76488E-09)
    assert_equal(ts1.total_heat, 2.57029E-09)
    assert_equal(ts1.total_heat_no_trit, 2.57029E-09)
    assert_equal(ts1.num_nuclides, 163)
    assert_equal(ts1.num_fissions, 2.61456E+11)
    assert_equal(ts1.neutron_flux, 1.00000E+06)
    assert_equal(ts1.initial_mass, 1.0)
    assert_equal(ts1.total_mass, 1.0)
    assert_equal(ts1.density, 19.3)
    assert_equal(ts1.actinide_burn, 8.55148E-12)

    assert_equal(ts3.appm_h1, 5.0059E-08)
    assert_equal(ts3.appm_h2, 9.0075E-09)
    assert_equal(ts3.appm_h3, 3.5742E-09)
    assert_equal(ts3.appm_he3, 7.5389E-09)
    assert_equal(ts3.appm_he4, 1.0894E-08)


def test_read_spectra():
    """test read of spectra data each time step for fispact-II """
    ts1 = fo.timestep_data[1]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]
    assert_equal(ts3.gspec[0], 1.13205E+01)
    assert_equal(len(ts3.gspec), 24)




def test_read_dominant():
    """test read of dominant data each time step for fispact-II"""
    ts1 = fo.timestep_data[1]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]

    assert_equal(len(ts1.dom_data[0]), 1)
    assert_equal(len(ts2.dom_data[0]), 1)
    assert_equal(len(ts3.dom_data[0]), 1)

    assert_equal(len(ts1.dom_data[0]), len(ts1.dom_data[1]))
    assert_equal(len(ts1.dom_data[1]), len(ts1.dom_data[2]))

    assert_equal(ts1.dom_data[0][0], " ")
    assert_equal(ts1.dom_data[1][0], 1)
    assert_equal(ts1.dom_data[2][0], 1)
    assert_equal(ts1.dom_data[3][0], " ")
    assert_equal(ts1.dom_data[4][0], 1)
    assert_equal(ts1.dom_data[5][0], 1)
    assert_equal(ts1.dom_data[6][0], " ")
    assert_equal(ts1.dom_data[7][0], 1)
    assert_equal(ts1.dom_data[8][0], 1)
    assert_equal(ts1.dom_data[9][0], " ")
    assert_equal(ts1.dom_data[10][0], 1)
    assert_equal(ts1.dom_data[11][0], 1)
    assert_equal(ts1.dom_data[12][0], " ")
    assert_equal(ts1.dom_data[13][0], 1)
    assert_equal(ts1.dom_data[14][0], 1)

    assert_equal(ts1.dom_data[0][-1], "rest")
    assert_equal(ts2.dom_data[0][-1], "rest")
    assert_equal(ts3.dom_data[0][-1], "rest")

def test_read_composition():
    """test read of composition data each time step for fispact-II"""
    ts1 = fo.timestep_data[1]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]

    assert_equal(len(ts1.composition[0]), 1)
    assert_equal(len(ts2.composition[0]), 1)
    assert_equal(len(ts3.composition[0]), 1)

def test_read_inv():
    """test read of inventory data for each time step for fispact-II """
    ts1 = fo.timestep_data[1]
    ts2 = fo.timestep_data[5]
    ts3 = fo.timestep_data[-1]
