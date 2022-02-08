"""Tally tests"""
import os
import warnings

from unittest import TestCase
import nose

from nose.tools import (
    assert_equal,
    assert_not_equal,
    assert_raises,
    raises,
    assert_almost_equal,
    assert_true,
    assert_false,
    assert_in,
)

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne.particle import is_valid, name, is_heavy_ion, id, mcnp, mcnp6, fluka, geant4


def test_is_valid():
    assert_equal(is_valid("Proton"), True)
    assert_equal(is_valid("Protium"), True)
    assert_equal(is_valid("Hydrogen"), True)
    assert_equal(is_valid(2212), True)
    assert_equal(is_valid(-2212), True)
    assert_equal(is_valid("Neutron"), True)
    assert_equal(is_valid("AntiProton"), True)
    assert_equal(is_valid("AntiNeutron"), True)


def test_name():
    assert_equal(name("Proton"), "Proton")
    assert_equal(name("Neutron"), "Neutron")
    assert_equal(name(2112), "Neutron")
    assert_equal(name("Hydrogen"), "Proton")
    assert_equal(name("Protium"), "Proton")
    assert_equal(name(10010000), "Proton")
    assert_equal(name("10010000"), "Proton")


def test_is_heavy_ion():
    assert_equal(is_heavy_ion("Proton"), False)
    assert_equal(is_heavy_ion("Hydrogen"), False)
    assert_equal(is_heavy_ion("1H"), False)
    assert_equal(is_heavy_ion(10010000), False)
    assert_equal(is_heavy_ion("10010000"), False)
    assert_equal(is_heavy_ion("2H"), True)
    assert_equal(is_heavy_ion("3He"), True)
    assert_equal(is_heavy_ion("22Na"), True)


def test_id_number():
    assert_equal(id("Proton"), 2212)
    assert_not_equal(id("AntiProton"), 2212)
    assert_equal(id("AntiProton"), -2212)
    assert_not_equal(id("Proton"), -2212)
    assert_equal(id("Hydrogen"), 2212)
    assert_equal(id("22Na"), 0)


def test_mcnp_id():
    assert_equal(mcnp("Neutron"), "n")
    assert_equal(mcnp(2112), "n")
    assert_equal(mcnp("Photon"), "p")
    assert_equal(mcnp("Gamma"), "p")
    assert_equal(mcnp("Electron"), "e")


def test_mcnp6_id():
    assert_equal(mcnp6("Neutron"), "n")
    assert_equal(mcnp6(2112), "n")
    assert_equal(mcnp6("Photon"), "p")
    assert_equal(mcnp6("Gamma"), "p")
    assert_equal(mcnp6("Electron"), "e")
    assert_equal(mcnp6("Proton"), "h")
    assert_equal(mcnp6("Hydrogen"), "h")


def test_fluka_id():
    assert_equal(fluka("Neutron"), "NEUTRON")
    assert_equal(fluka(2112), "NEUTRON")
    assert_equal(fluka("Photon"), "PHOTON")
    assert_equal(fluka("Gamma"), "PHOTON")
    assert_equal(fluka("Electron"), "ELECTRON")
    assert_equal(fluka("Beta-"), "ELECTRON")
    assert_equal(fluka("Proton"), "PROTON")
    assert_equal(fluka("Hydrogen"), "PROTON")


def test_geant4_id():
    assert_equal(geant4("Neutron"), "neutron")
    assert_equal(geant4(2112), "neutron")
    assert_equal(geant4("Photon"), "gamma")
    assert_equal(geant4("Gamma"), "gamma")
    assert_equal(geant4("Electron"), "e-")
    assert_equal(geant4("Beta-"), "e-")
    assert_equal(geant4("Proton"), "proton")
    assert_equal(geant4("Hydrogen"), "proton")


# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
