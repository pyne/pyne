import os

from nose.tools import assert_true, assert_false, assert_in, assert_equal
from pyne.dbgen.materials_library import *
from pyne.material import Material
from pyne.material_library import MaterialLibrary
from pyne.pyne_config import pyne_conf
import pyne.nucname as nucname
import tables as tb
import numpy as np


def assert_close(obs, exp, margin=1e-6):
    # default margin is because we only have 1e-6 of precision
    assert(abs(obs-exp) < margin)


def test_is_comp_matname_or_density():
    comp = ["H", 1001, 1000, 0.101327, 0.583640, 0.068228, None]
    matname = ["1.  ", "A-150 Tissue-Equivalent Plastic (A150TEP)", None, None, None, None, None]
    density = ["Density (g/cm3) =", None, 1.127000, None, "Total atoms/b-cm =", None, 1.169E-01]
    badlines = [
        ["Formula =", "-", None, None, "Molecular weight (g/mole) =", None, "-"],
        ["The above density is estimated to be accurate to 4 significant digits.  Dealing with uncertainties is left to the user.", None, None, None, None, None, None],
        ["The following data was calculated from the input weight fractions.", None, None, None, None, None, None],
        [None, None, None, None, None, None, None],
        [None, None, None, "Weight", "Atom", "Atom", None],
        ["Element", "Neutron ZA", "Photon ZA", "Fraction", "Fraction", "Density",  None],
        ]
    assert_true(is_comp_matname_or_density(comp))
    assert_true(is_comp_matname_or_density(matname))
    assert_true(is_comp_matname_or_density(density))
    for line in badlines:
        assert_false(is_comp_matname_or_density(line))


def test_grab_materials_compendium():
    mats = grab_materials_compendium('../pyne/dbgen/materials_compendium.csv')
    assert(len(mats) == 372)

    # this tests a material where we don't do any element expansion
    a150tep_comp = mats["A-150 Tissue-Equivalent Plastic (A150TEP)"].comp
    expected_mat = Material({"H": 0.101327, "C": 0.775501, "N": 0.035057,
                             "O": 0.052316, "F": 0.017422, "Ca": 0.018378})
    expected_mat.normalize()
    expected_mat = expected_mat.expand_elements()
    for key, value in expected_mat.comp.items():
        assert_close(a150tep_comp[key], value)

    # this tests a material where we do do element expansion
    pubr = mats["Plutonium Bromide"]
    bromium = sum((frac for nuc, frac in pubr.comp.items() if nucname.zzzaaa(nuc) // 1000 == 35))
    assert_close(bromium, 0.500617)
    assert_close(pubr[942380000], 0.000250)


def test_against_nuc_data():
    nuc_data = pyne_conf.NUC_DATA_PATH
    if not os.path.isfile(nuc_data):
        raise RuntimeError("Tests require nuc_data.h5. Please run nuc_data_make.")
    obs_matslib = MaterialLibrary(nuc_data,
                                  datapath="/material_library/materials",
                                  nucpath="/material_library/nucid")
    gasoline = Material({
        "H": 0.157000,
        "C": 0.843000,
        },
        density=0.721,
        metadata={"name": "Gasoline"}).expand_elements()

    pubr3 = Material({
        "Br": 0.500617,
        "Pu-238": 0.000250,
        "Pu-239": 0.466923,
        "Pu-240": 0.029963,
        "Pu-241": 0.001998,
        "Pu-242": 0.000250
        },
        density=6.75,
        metadata={"name": "Plutonium Bromide"}).expand_elements()

    obs_gasoline = obs_matslib["Gasoline"]
   
    # remove empty elements
    assert_equal( { (x,val) for x,val in set(obs_gasoline.comp.items()) if val > 0 }, 
            set(gasoline.comp.items()))
    assert_equal(obs_gasoline.density, gasoline.density)
    assert_equal(obs_gasoline.metadata["name"], gasoline.metadata["name"])

    obs_pubr3 = obs_matslib["Plutonium Bromide"]
    # remove empty elements
    assert_equal( { (x,val) for x,val in set(obs_pubr3.comp.items()) if val > 0 }, 
            set(pubr3.comp.items()))
    assert_equal(obs_pubr3.density, pubr3.density)
    assert_equal(obs_pubr3.metadata["name"], pubr3.metadata["name"])
