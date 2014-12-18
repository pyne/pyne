import os

from nose.tools import assert_true, assert_false, assert_in
from pyne.dbgen.materials_library import *
from pyne.material import Material, MaterialLibrary
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


def get_comps_from_nuc_data(nuc_data):
    with tb.openFile(nuc_data, 'r') as f:
        comps = np.array([m['comp'] for m in f.root.material_library.materials])
        nucids = np.array(f.root.material_library.nucid)
    comps_zipped = (zip(nucids, comp) for comp in comps)
    comps_cleaned = ((p for p in comp if p[1] != 0) for comp in comps_zipped)
    obs_comps = filter(lambda x: x != {}, [dict(comp) for comp in comps_cleaned])
    return obs_comps


def test_output_h5():
    nuc_data = "test_nd.h5"
    if os.path.isfile(nuc_data):
        os.remove(nuc_data)
    matslib = make_matslib(os.path.join(os.path.split(__file__)[0],
                                        '../pyne/dbgen/materials_compendium.csv'))
    make_materials_compendium(nuc_data, matslib)

    obs_comps = get_comps_from_nuc_data(nuc_data)
    os.remove(nuc_data)
    exp_comps = [mat.comp for mat in list(matslib.values())]
    for obs_comp in obs_comps:
        assert(obs_comp in exp_comps)


def test_against_nuc_data():
    nuc_data = pyne_conf.NUC_DATA_PATH
    if not os.path.isfile(nuc_data):
        raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")
    matslib = make_matslib(os.path.join(os.path.split(__file__)[0],
                                        '../pyne/dbgen/materials_compendium.csv'))

    obs_comps = get_comps_from_nuc_data(nuc_data)
    exp_comps = [mat.comp for mat in list(matslib.values())]
    for obs_comp in obs_comps:
        assert_in(obs_comp, exp_comps)
