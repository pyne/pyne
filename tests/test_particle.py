"""Tally tests"""
import os
import warnings

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne.particle import is_valid,name,is_heavy_ion,pdc_number


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


# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
