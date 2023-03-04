import os
import numpy as np
import tables as tb
import pytest

import pyne.nucname as nucname
from pyne.pyne_config import pyne_conf
from pyne.material_library import MaterialLibrary
from pyne.material import Material
from pyne.dbgen.materials_library import *
from pyne import nuc_data

nucvec = {
    10010000: 1.0,
    80160000: 1.0,
    691690000: 1.0,
    922350000: 1.0,
    922380000: 1.0,
    942390000: 1.0,
    942410000: 1.0,
    952420000: 1.0,
    962440000: 1.0,
}
leu = {922380000: 0.96, 922350000: 0.04}


def assert_mat_almost_equal(first, second, places=1E-7):
    assert first.mass == pytest.approx(second.mass, rel=places)
    assert first.density == pytest.approx(second.density, rel=places)
    assert first.atoms_per_molecule == pytest.approx(second.atoms_per_molecule, rel=places)
    assert first.metadata["name"] == second.metadata["name"]
    nucs = set(second.comp)
    assert set(first.comp) == nucs
    for nuc in nucs:
        assert first.comp[nuc] == pytest.approx(second.comp[nuc], rel=places)


def assert_close(obs, exp, margin=1e-6):
    # default margin is because we only have 1e-6 of precision
    assert abs(obs - exp) < margin


def test_is_comp_matname_or_density():
    comp = ["H", 1001, 1000, 0.101327, 0.583640, 0.068228, None]
    matname = [
        "1.  ",
        "A-150 Tissue-Equivalent Plastic (A150TEP)",
        None,
        None,
        None,
        None,
        None,
    ]
    density = [
        "Density (g/cm3) =",
        None,
        1.127000,
        None,
        "Total atoms/b-cm =",
        None,
        1.169e-01,
    ]
    badlines = [
        ["Formula =", "-", None, None, "Molecular weight (g/mole) =", None, "-"],
        [
            "The above density is estimated to be accurate to 4 significant digits.  Dealing with uncertainties is left to the user.",
            None,
            None,
            None,
            None,
            None,
            None,
        ],
        [
            "The following data was calculated from the input weight fractions.",
            None,
            None,
            None,
            None,
            None,
            None,
        ],
        [None, None, None, None, None, None, None],
        [None, None, None, "Weight", "Atom", "Atom", None],
        ["Element", "Neutron ZA", "Photon ZA", "Fraction", "Fraction", "Density", None],
    ]
    assert is_comp_matname_or_density(comp)
    assert is_comp_matname_or_density(matname)
    assert is_comp_matname_or_density(density)
    for line in badlines:
        assert not is_comp_matname_or_density(line)


def test_grab_materials_compendium():
    mats = grab_materials_compendium("../pyne/dbgen/materials_compendium.csv")
    assert len(mats) == 372

    # this tests a material where we don't do any element expansion
    a150tep_comp = mats["A-150 Tissue-Equivalent Plastic (A150TEP)"].comp
    expected_mat = Material(
        {
            "H": 0.101327,
            "C": 0.775501,
            "N": 0.035057,
            "O": 0.052316,
            "F": 0.017422,
            "Ca": 0.018378,
        }
    )
    expected_mat.normalize()
    expected_mat = expected_mat.expand_elements()
    for key, value in expected_mat.comp.items():
        assert_close(a150tep_comp[key], value)

    # this tests a material where we do do element expansion
    pubr = mats["Plutonium Bromide"]
    bromium = sum(
        (frac for nuc, frac in pubr.comp.items() if nucname.zzzaaa(nuc) // 1000 == 35)
    )
    assert_close(bromium, 0.500617)
    assert_close(pubr[942380000], 0.000250)


def test_against_nuc_data():
    nuc_data = pyne_conf.NUC_DATA_PATH
    if not os.path.isfile(nuc_data):
        raise RuntimeError("Tests require nuc_data.h5. Please run nuc_data_make.")
    obs_matslib = MaterialLibrary(nuc_data)
    gasoline = Material(
        {
            "H": 0.157000,
            "C": 0.843000,
        },
        density=0.721,
        metadata={"name": "Gasoline"},
    ).expand_elements()

    pubr3 = Material(
        {
            "Br": 0.500617,
            "Pu-238": 0.000250,
            "Pu-239": 0.466923,
            "Pu-240": 0.029963,
            "Pu-241": 0.001998,
            "Pu-242": 0.000250,
        },
        density=6.75,
        metadata={"name": "Plutonium Bromide"},
    ).expand_elements()

    obs_gasoline = obs_matslib["Gasoline"]
    # remove empty elements
    assert (
        {(x, val) for x, val in set(obs_gasoline.comp.items()) if val > 0} ==
        set(gasoline.comp.items()))
    assert obs_gasoline.density == gasoline.density
    assert obs_gasoline.metadata["name"] == gasoline.metadata["name"]

    obs_pubr3 = obs_matslib["Plutonium Bromide"]
    # remove empty elements
    assert (
        {(x, val) for x, val in set(obs_pubr3.comp.items()) if val > 0} ==
        set(pubr3.comp.items()))
    assert obs_pubr3.density == pubr3.density
    assert obs_pubr3.metadata["name"] == pubr3.metadata["name"]


def test_matlib_json():
    filename = "matlib.json"
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."
    lib = {"leu": Material(leu), "nucvec": nucvec, "aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_json(filename)

    rmatlib = MaterialLibrary()
    rmatlib.from_json(filename)
    assert set(wmatlib) == set(rmatlib)
    for key in rmatlib:
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)


def test_matlib_hdf5_nuc_data():
    matlib = MaterialLibrary()
    matlib.from_hdf5(nuc_data, "/materials")
    matlib.write_hdf5("matlib_test.h5")
    mat_lib_load_test = MaterialLibrary("matlib_test.h5")
    os.remove("matlib_test.h5")
    for key in matlib:
        assert_mat_almost_equal(matlib[key], mat_lib_load_test[key])


def test_matlib_hdf5():
    filename = "matlib.h5"
    if filename in os.listdir("."):
        os.remove(filename)
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."
    lib = {"leu": Material(leu), "nucvec": nucvec, "aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_hdf5(filename, "/mats1")
    rmatlib = MaterialLibrary()
    rmatlib.from_hdf5(filename, "/mats1")
    os.remove(filename)
    # Round trip!
    rmatlib.write_hdf5(filename, "/mats1")
    wmatlib = MaterialLibrary(filename, "/mats1")
    assert set(wmatlib) == set(rmatlib)
    for key in rmatlib:
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)


def test_hdf5_overwrite():
    filename = "matlib.h5"
    if filename in os.listdir("."):
        os.remove(filename)
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."
    lib = {"aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_hdf5(filename, "/mats1")

    lib = {"leu": Material(leu)}
    umatlib = MaterialLibrary(lib)
    # test error raise if path already exists
    pytest.raises(RuntimeError, umatlib.write_hdf5, filename, "/mats1", False)

    # test overwriting
    umatlib.write_hdf5(filename, "/mats1", h5_overwrite=True)
    rmatlib = MaterialLibrary()
    rmatlib.from_hdf5(filename, "/mats1")

    # Round trip!
    assert set(umatlib) == set(rmatlib)
    for key in rmatlib:
        assert_mat_almost_equal(umatlib[key], rmatlib[key])
    os.remove(filename)


def test_matlib_query():
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."
    mat_nucvec = Material(nucvec)
    mat_nucvec.metadata["name"] = "nucvec"
    lib = {"nucvec": nucvec, "aqua": water}
    matlib = MaterialLibrary(lib)

    matlib_aqua = matlib["aqua"]
    assert_mat_almost_equal(water, matlib_aqua)

    matlib_nucvec = matlib["nucvec"]
    assert_mat_almost_equal(mat_nucvec, matlib_nucvec)

    mat_leu = Material(leu)
    mat_leu.metadata["name"] = "leu"
    matlib["leu"] = mat_leu
    matlib_leu = matlib["leu"]
    assert_mat_almost_equal(mat_leu, matlib_leu)


def test_matlib_delete():
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."

    pubr3 = Material(
        {
            "Br": 0.500617,
            "Pu-238": 0.000250,
            "Pu-239": 0.466923,
            "Pu-240": 0.029963,
            "Pu-241": 0.001998,
            "Pu-242": 0.000250,
        },
        density=6.75,
        metadata={"name": "Plutonium Bromide"},
    ).expand_elements()

    lib = {"pubr3": pubr3, "aqua": water}
    matlib = MaterialLibrary(lib)

    # test delete methods
    matlib.remove_material("pubr3")

    assert len(matlib) == 1

    del matlib["aqua"]

    assert len(matlib) == 0
