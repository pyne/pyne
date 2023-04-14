"""Tally tests"""
import os
import warnings

from unittest import TestCase


from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne.particle import is_valid, name, is_heavy_ion, id, mcnp, mcnp6, fluka, geant4


def test_is_valid():
    assert is_valid("Proton") == True
    assert is_valid("Protium") == True
    assert is_valid("Hydrogen") == True
    assert is_valid(2212) == True
    assert is_valid(-2212) == True
    assert is_valid("Neutron") == True
    assert is_valid("AntiProton") == True
    assert is_valid("AntiNeutron") == True


def test_name():
    assert name("Proton") == "Proton"
    assert name("Neutron") == "Neutron"
    assert name(2112) == "Neutron"
    assert name("Hydrogen") == "Proton"
    assert name("Protium") == "Proton"
    assert name(10010000) == "Proton"
    assert name("10010000") == "Proton"


def test_is_heavy_ion():
    assert is_heavy_ion("Proton") == False
    assert is_heavy_ion("Hydrogen") == False
    assert is_heavy_ion("1H") == False
    assert is_heavy_ion(10010000) == False
    assert is_heavy_ion("10010000") == False
    assert is_heavy_ion("2H") == True
    assert is_heavy_ion("3He") == True
    assert is_heavy_ion("22Na") == True


def test_id_number():
    assert id("Proton") == 2212
    assert id("AntiProton") != 2212
    assert id("AntiProton") == -2212
    assert id("Proton") != -2212
    assert id("Hydrogen") == 2212
    assert id("22Na") == 0


def test_mcnp_id():
    assert mcnp("Neutron") == "n"
    assert mcnp(2112) == "n"
    assert mcnp("Photon") == "p"
    assert mcnp("Gamma") == "p"
    assert mcnp("Electron") == "e"


def test_mcnp6_id():
    assert mcnp6("Neutron") == "n"
    assert mcnp6(2112) == "n"
    assert mcnp6("Photon") == "p"
    assert mcnp6("Gamma") == "p"
    assert mcnp6("Electron") == "e"
    assert mcnp6("Proton") == "h"
    assert mcnp6("Hydrogen") == "h"


def test_fluka_id():
    assert fluka("Neutron") == "NEUTRON"
    assert fluka(2112) == "NEUTRON"
    assert fluka("Photon") == "PHOTON"
    assert fluka("Gamma") == "PHOTON"
    assert fluka("Electron") == "ELECTRON"
    assert fluka("Beta-") == "ELECTRON"
    assert fluka("Proton") == "PROTON"
    assert fluka("Hydrogen") == "PROTON"


def test_geant4_id():
    assert geant4("Neutron") == "neutron"
    assert geant4(2112) == "neutron"
    assert geant4("Photon") == "gamma"
    assert geant4("Gamma") == "gamma"
    assert geant4("Electron") == "e-"
    assert geant4("Beta-") == "e-"
    assert geant4("Proton") == "proton"
    assert geant4("Hydrogen") == "proton"

