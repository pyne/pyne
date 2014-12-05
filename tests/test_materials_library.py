from pyne.dbgen.materials_library import *
from pyne.material import Material
from nose.tools import assert_true, assert_false, assert_equal
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


def test_elem_line_to_mat():
    mass = 0.101327
    nuclide = ["H", 1001, 1000, 0.101327, 0.583640, 0.068228, None]
    assert_equal(elem_line_to_mat(nuclide), Material({1001: 1.0}, mass))

    no_expand = ["H", 1000, 1000, 0.101327, 0.583640, 0.068228, None]
    assert_equal(elem_line_to_mat(no_expand), Material({1000: 1.0}, mass))

    expand = ["H", "-", 1000, 0.101327, 0.583640, 0.068228, None]
    assert_equal(elem_line_to_mat(expand).comp[10010000], 0.999885)
    assert_equal(elem_line_to_mat(expand).comp[10020000], 0.000115)


def test_grab_materials_compendium():
    mats = grab_materials_compendium('../pyne/dbgen/materials_compendium.csv')
    assert(len(mats) == 372)

    # this tests a material where we don't do any element expansion
    a150tep = mats[0]
    # gotta normalize our expectation!
    assert_close(a150tep.comp[nucname.id(1001)], (0.101327 / 1.000001))
    assert_close(a150tep.comp[nucname.id(6000)], (0.775501 / 1.000001))
    assert_close(a150tep.comp[nucname.id(7014)], (0.035057 / 1.000001))
    assert_close(a150tep.comp[nucname.id(8016)], (0.052316 / 1.000001))
    assert_close(a150tep.comp[nucname.id(9019)], (0.017422 / 1.000001))
    assert_close(a150tep.comp[nucname.id(20000)], (0.018378 / 1.000001))

    # this tests a material where we do do element expansion
    pubr = [mat for mat in mats if mat.metadata["name"] == "Plutonium Bromide"][0]
    bromium = sum((frac for nuc, frac in pubr.comp.items() if nucname.zzzaaa(nuc) // 1000 == 35))
    assert_close(bromium, 0.500617)
