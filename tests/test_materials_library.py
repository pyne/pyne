from pyne.dbgen.materials_library import *
from pyne.material import Material
from nose.tools import assert_true, assert_false
import pyne.nucname as nucname


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
