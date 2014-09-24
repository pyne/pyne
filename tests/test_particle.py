"""Tally tests"""
import os
import warnings

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne.particle import is_valid,name,is_heavy_ion,pdc_number,mcnp,mcnp6, \
    fluka


def test_is_valid():
    assert_equal(is_valid("Proton"),True)
    assert_equal(is_valid("Protium"),True)
    assert_equal(is_valid("Hydrogen"),True)
    assert_equal(is_valid(2212),True)
    assert_equal(is_valid(-2212),True)
    assert_equal(is_valid("Neutron"),True)
    assert_equal(is_valid("AntiProton"),True)
    assert_equal(is_valid("AntiNeutron"),True)

def test_name():
    assert_equal(name("Proton"),"Proton")
    assert_equal(name("Neutron"),"Neutron")
    assert_equal(name(2112),"Neutron")
    assert_equal(name("Hydrogen"),"Proton")
    assert_equal(name("Protium"),"Proton")
    assert_equal(name(10010000),"Proton")
    assert_equal(name("10010000"),"Proton")

def test_is_heavy_ion():
    assert_equal(is_heavy_ion("Proton"),False)
    assert_equal(is_heavy_ion("Hydrogen"),False)
    assert_equal(is_heavy_ion("1H"),False)
    assert_equal(is_heavy_ion(10010000),False)
    assert_equal(is_heavy_ion("10010000"),False)
    assert_equal(is_heavy_ion("2H"),True)
    assert_equal(is_heavy_ion("3He"),True)
    assert_equal(is_heavy_ion("22Na"),True)

def test_pdc_number():
    assert_equal(pdc_number("Proton"),2212)
    assert_not_equal(pdc_number("AntiProton"),2212)
    assert_equal(pdc_number("AntiProton"),-2212)
    assert_not_equal(pdc_number("Proton"),-2212)
    assert_equal(pdc_number("Hydrogen"),2212)
    assert_equal(pdc_number("22Na"),0)

def test_mcnp_id():
    assert_equal(mcnp("Neutron"),"N")
    assert_equal(mcnp(2112),"N")
    assert_equal(mcnp("Photon"),"P")
    assert_equal(mcnp("Gamma"),"P")
    assert_equal(mcnp("Electron"),"E")

def test_mcnp6_id():
    assert_equal(mcnp6("Neutron"),"N")
    assert_equal(mcnp6(2112),"N")
    assert_equal(mcnp6("Photon"),"P")
    assert_equal(mcnp6("Gamma"),"P")
    assert_equal(mcnp6("Electron"),"E")
    assert_equal(mcnp6("Proton"),"H")
    assert_equal(mcnp6("Hydrogen"),"H")

def test_fluka_id():
    assert_equal(fluka("Neutron"),"NEUTRON")
    assert_equal(fluka(2112),"NEUTRON")
    assert_equal(fluka("Photon"),"PHOTON")
    assert_equal(fluka("Gamma"),"PHOTON")
    assert_equal(fluka("Electron"),"ELECTRON")
    assert_equal(fluka("Beta-"),"ELECTRON")
    assert_equal(fluka("Proton"),"PROTON")
    assert_equal(fluka("Hydrogen"),"PROTON")


# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
