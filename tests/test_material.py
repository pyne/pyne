"""Material tests"""
import os
from copy import deepcopy
import warnings

from unittest import TestCase
import pytest

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne import nuc_data
from pyne.material import (
    Material,
    from_atom_frac,
    from_activity,
    from_hdf5,
    from_text,
    MapStrMaterial,
    MultiMaterial,
)
from pyne import jsoncpp
from pyne import data
from pyne import nucname
from pyne import utils
from pyne import cram
import numpy as np
from numpy.testing import assert_array_equal
import tables as tb

if utils.use_warnings():
    utils.toggle_warnings()

nclides = 9
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


def assert_mat_almost_equal(first, second, places=1e-7):
    assert first.mass == pytest.approx(second.mass, rel=places)
    assert first.density == pytest.approx(second.density, rel=places)
    assert first.atoms_per_molecule == pytest.approx(
        second.atoms_per_molecule, rel=places
    )
    assert first.metadata == second.metadata
    nucs = set(second.comp)
    assert set(first.comp) == nucs
    for nuc in nucs:
        assert first.comp[nuc] == pytest.approx(second.comp[nuc], rel=places)


def make_mat_txt():
    """Helper for mat.txt"""
    with open("mat.txt", "w") as f:
        f.write("U235  0.05\nU238  0.95")


def make_mat_h5():
    """Helper for mat.h5"""
    f = tb.open_file("mat.h5", "w")
    f.create_group("/", "mat", "Mass Material Test")
    f.create_array("/mat", "Mass", np.array([1.0, 0.5, 0.0]), "Mass Test")
    f.create_array("/mat", "U235", np.array([1.0, 0.75, 0.0]), "U235 Test")
    f.create_array("/mat", "PU239", np.array([0.0, 0.25, 0.0]), "PU239 Test")
    f.close()


###############################################################################


def test_mat1():
    mat = Material("mat.txt")
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.mass == 1.0


def test_mat2():
    mat = Material("mat.txt", 42)
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.mass == 42.0


def test_mat3():
    mat = Material("mat.txt", -42)
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.mass == 1.0


def test_mat4():
    mat = Material({922350000: 0.05, 922380000: 0.95}, 15, metadata={"units": "kg"})
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.mass == 15.0
    assert mat.metadata["units"] == "kg"


def test_from_text():
    mat = Material(metadata={"units": "kg"})
    mat.from_text("mat.txt")
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.metadata["units"] == "kg"


def test_from_textdup():
    mat = Material(metadata={"units": "kg"})
    mat.from_text("matdup.txt")
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}
    assert mat.metadata["units"] == "kg"


def test_from_textelem():
    mat = Material(metadata={"units": "kg"})
    mat.from_text("matelem.txt")
    assert mat.comp == {10000000: 0.1, 80000000: 0.9}
    assert mat.metadata["units"] == "kg"


def test_write_text():
    if "leu.txt" in os.listdir("."):
        os.remove("leu.txt")

    leu = Material({"U235": 0.04, "U238": 0.96}, 42.0, 1.0, 1.0)
    leu.write_text("leu.txt")

    with open("leu.txt") as f:
        written = f.read()
    expected = (
        "Mass    42\n" "Density 1\n" "APerM   1\n" "U235    0.04\n" "U238    0.96\n"
    )
    assert written == expected

    read_leu = from_text("leu.txt")
    assert leu.mass == read_leu.mass
    assert leu.atoms_per_molecule == read_leu.atoms_per_molecule
    assert leu.comp == read_leu.comp

    os.remove("leu.txt")


def test_from_hdf5_protocol_0():
    mat = Material()
    mat.from_hdf5("mat.h5", "/mat", protocol=0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat.from_hdf5("mat.h5", "/mat", 0, 0)
    assert mat.mass == 1.0
    assert mat.comp == {922350000: 1.0, 942390000: 0.0}

    mat.from_hdf5("mat.h5", "/mat", 1, 0)
    assert mat.mass == 0.5
    assert mat.comp == {922350000: 0.75, 942390000: 0.25}

    mat.from_hdf5("mat.h5", "/mat", 2, 0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat.from_hdf5("mat.h5", "/mat", -1, 0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat.from_hdf5("mat.h5", "/mat", -2, 0)
    assert mat.mass == 0.5
    assert mat.comp == {922350000: 0.75, 942390000: 0.25}

    mat.from_hdf5("mat.h5", "/mat", -3, 0)
    assert mat.mass == 1.0
    assert mat.comp == {922350000: 1.0, 942390000: 0.0}


def test_hdf5_protocol_1_old():
    if "proto1.h5" in os.listdir("."):
        os.remove("proto1.h5")

    # Test material writing
    leu = Material({"U235": 0.04, "U238": 0.96}, 4.2, 2.72, 1.0)
    leu.metadata["comment"] = "first light"
    leu.write_hdf5("proto1.h5", nucpath="/nucid", chunksize=10)

    for i in range(2, 11):
        leu = Material({"U235": 0.04, "U238": 0.96}, i * 4.2, 2.72, 1.0 * i)
        leu.metadata["comment"] = "fire in the disco - {0}".format(i)
        leu.write_hdf5("proto1.h5")

    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5("proto1.h5", "/mat_name", -3, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 8.0
    assert m.mass == 33.6
    assert m.comp == {922350000: 0.04, 922380000: 0.96}
    assert m.metadata["comment"] == "fire in the disco - 8"

    m = from_hdf5("proto1.h5", "/mat_name", 3, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 4.0
    assert m.mass == 16.8
    assert m.comp == {922350000: 0.04, 922380000: 0.96}
    assert m.metadata["comment"] == "fire in the disco - 4"

    # Test mateiral old format detection
    leu = Material({"U235": 0.02, "U238": 0.98}, 6.34, 2.72, 1.0)
    leu.metadata["comment"] = "hitting the dancefloor"
    leu.write_hdf5("proto1.h5", datapath="/new_mat", chunksize=10)

    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5("proto1.h5", "/new_mat", -0, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 1.0
    assert m.mass == 6.34
    assert m.comp == {922350000: 0.02, 922380000: 0.98}
    assert m.metadata["comment"] == "hitting the dancefloor"

    os.remove("proto1.h5")

    leu = Material({"U235": 0.02, "U238": 0.98}, 6.34, 2.72, 1.0)
    leu.metadata["comment"] = "hitting the dancefloor"
    leu.write_hdf5("proto1.h5", datapath="/new_mat", nucpath="/nucid", chunksize=10)
    # Test material writing
    leu = Material({"U235": 0.04, "U238": 0.96}, 4.2, 2.72, 1.0)
    leu.metadata["comment"] = "first light"
    leu.write_hdf5("proto1.h5", datapath="/new_mat", chunksize=10)
    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5("proto1.h5", "/new_mat", -0, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 1.0
    assert m.mass == 6.34
    assert m.comp == {922350000: 0.02, 922380000: 0.98}
    assert m.metadata["comment"] == "hitting the dancefloor"
    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5("proto1.h5", "/new_mat", -1, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 1.0
    assert m.mass == 4.2
    assert m.comp == {922350000: 0.04, 922380000: 0.96}
    assert m.metadata["comment"] == "first light"

    os.remove("proto1.h5")


def test_hdf5_protocol_1():
    if "proto1.h5" in os.listdir("."):
        os.remove("proto1.h5")

    # Test material writing
    leu = Material({"U235": 0.04, "U238": 0.96}, 4.2, 2.72, 1.0)
    leu.metadata["comment"] = "first light"
    leu.write_hdf5("proto1.h5")

    for i in range(2, 11):
        leu = Material({"U235": 0.04, "U238": 0.96}, i * 4.2, 2.72, 1.0 * i)
        leu.metadata["comment"] = "fire in the disco - {0}".format(i)
        leu.write_hdf5("proto1.h5")

    # Loads with protocol 1 now.
    for i in range(2, 11):
        m = Material()
        m.from_hdf5("proto1.h5", "/mat_name", i - 1, 1)
        assert m.density == 2.72
        assert m.atoms_per_molecule == 1.0 * i
        assert m.mass == i * 4.2
        assert m.comp == {922350000: 0.04, 922380000: 0.96}
        assert m.metadata["comment"] == "fire in the disco - {0}".format(i)

    m = from_hdf5("proto1.h5", "/mat_name", -1, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 10.0
    assert m.mass == 42.0
    assert m.comp == {922350000: 0.04, 922380000: 0.96}
    assert m.metadata["comment"] == "fire in the disco - 10"

    m = from_hdf5("proto1.h5", "/mat_name", 0, 1)
    assert m.density == 2.72
    assert m.atoms_per_molecule == 1.0
    assert m.mass == 4.2
    assert m.comp == {922350000: 0.04, 922380000: 0.96}
    assert m.metadata["comment"] == "first light"
    os.remove("proto1.h5")


class TestMaterialMethods(TestCase):
    "Tests that the Material member functions work."

    def test_normalize(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        mat.normalize()
        assert mat.mass == 1.0

    def test_mult_by_mass(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        nucvec = mat.mult_by_mass()
        assert nucvec == {922350000: 0.75, 922380000: 14.25}

    def test_activity(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        obs = mat.activity()
        exp = {922350000: 59953.15151320379, 922380000: 177216.65219208886}
        assert set(obs) == set(exp)
        assert set(obs.values()) == set(exp.values())

    def test_from_activity_meth(self):
        mat = Material()
        mat.from_activity({922350000: 59953.15151320379, 922380000: 177216.65219208886})
        assert mat.atoms_per_molecule == -1.0
        assert mat.comp[922350000] == 0.05
        assert mat.comp[922380000] == 0.95
        assert mat.mass == 15.0
        assert mat.molecular_mass() == 237.8986180886295

        mat = Material()
        mat.from_activity({"H": 0.0, "H3": 1.0, "O16": 0.0})
        assert mat.comp == {10030000: 1.0}
        assert mat.mass == 2.809161396488766e-15

        mat = Material()
        hto = {"H": 1.0, "H3": 1.0, "O16": 1.0}
        pytest.raises(ValueError, mat.from_activity, hto)

    def test_decay_heat_stable(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        obs = mat.decay_heat()
        exp = {922350000: 4.48963565256e-14, 922380000: 1.2123912039e-13}
        assert set(obs) == set(exp)
        for key in exp:
            assert obs[key] == pytest.approx(exp[key])

    def test_decay_heat_metastable(self):
        mat = Material({471080001: 0.5, 551340001: 0.5}, 2)
        obs = mat.decay_heat()
        # decay heat values from external calculation using half-life: BNL
        # q-values: ORNL, atomic mass: AMDC
        exp = {471080001: 7.33578506845e-8, 551340001: 6.241109861949e-3}
        assert set(obs) == set(exp)
        for key in exp:
            assert obs[key] == pytest.approx(exp[key])

    def test_dose_per_g(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        # testing for default source
        obs1 = mat.dose_per_g("ext_air")
        exp1 = {922350000: 1.11264406283e-14, 922380000: 5.01315571163e-15}
        assert set(obs1) == set(exp1)
        for key in exp1:
            assert obs1[key] == pytest.approx(exp1[key])
        # testing for non-default source
        obs2 = mat.dose_per_g("ingest", 1)
        exp2 = {922350000: 27.11394777435299, 922380000: 77.59215574705213}
        assert set(obs2) == set(exp2)
        for key in exp2:
            assert obs2[key] == pytest.approx(exp2[key])

    def test_molecular_mass(self):
        mat_empty = Material({})
        assert mat_empty.molecular_mass() == 0.0

        mat_u238 = Material({922380000: 1.0})
        mw_u238 = mat_u238.molecular_mass()
        try:
            assert mw_u238 == pytest.approx(238.050786996)
        except AssertionError:
            assert mw_u238 == pytest.approx(238.0)

        mat_mixed = Material({922350000: 0.5, 922380000: 0.5})
        mw_mixed = mat_mixed.molecular_mass()
        try:
            assert mw_mixed / 236.547360417 == pytest.approx(1.0, rel=1e-4)
        except AssertionError:
            assert mw_mixed / 236.5 == pytest.approx(1.0, rel=1e-4)


def test_expand_elements1():
    natmat = Material(
        {"C": 1.0, 902320000: 0.5, "PU": 4.0, "U": 3.0}, metadata={"y": 1.0}
    )
    expmat = natmat.expand_elements()
    assert 60120000 in expmat.comp
    assert not 60000000 in expmat.comp
    assert natmat.metadata == expmat.metadata
    assert not natmat.metadata is expmat.metadata


def test_expand_elements2():
    """Inspired by #86"""
    natmat = Material({"C": 1.0})
    expmat = natmat.expand_elements()
    afrac = expmat.to_atom_frac()
    assert data.natural_abund(60120000) == pytest.approx(afrac[60120000])
    assert data.natural_abund(60130000) == pytest.approx(afrac[60130000])


def test_expand_elements3():
    natmat = Material({"C": 1.0})
    exception_ids = {nucname.id("C")}
    expmat = natmat.expand_elements(exception_ids)
    afrac = expmat.to_atom_frac()
    assert natmat[60000000] == pytest.approx(afrac[60000000])


def test_collapse_elements1():
    """Very simple test to combine nucids"""
    nucvec = {
        10010000: 1.0,
        80160000: 1.0,
        80160001: 1.0,
        691690000: 1.0,
        922350000: 1.0,
        922380000: 1.0,
        942390000: 1.0,
        952420000: 1.0,
        962440000: 1.0,
    }

    exception_ids = {
        nucname.id(1001),
        nucname.id("U-235"),
        nucname.id("U-238"),
        nucname.id("Pu-239"),
        nucname.id("Pu-241"),
    }

    mat = Material(nucvec)

    cmat = mat.collapse_elements(exception_ids)

    assert cmat.comp[80000000] == mat.comp[80160000] + mat.comp[80160001]
    assert cmat.comp[922350000] == mat.comp[922350000]
    assert cmat.comp[942390000] == mat.comp[942390000]
    assert cmat.comp[950000000] == mat.comp[952420000]
    assert cmat.comp[960000000] == mat.comp[952420000]


def test_mass_density():
    ethanol = from_atom_frac({"C": 2, "H": 6, "O": 1})
    atom_density_ethanol = 9.282542841e22  # atom density not molecule density
    mass_density = ethanol.mass_density(atom_density_ethanol)
    expected_mass_density = 0.78900
    assert mass_density == pytest.approx(expected_mass_density, rel=1e-4)


def test_number_density():
    ethanol = from_atom_frac({"C": 2, "H": 6, "O": 1}, density=0.78900)
    obs = ethanol.number_density()
    exp = 9.2825e22
    assert obs / exp == pytest.approx(1.0, rel=1e-4)


def test_set_mat_int_1():
    mat = Material(nucvec, -1)
    mat1 = mat.set_mat([922350000, 922380000, 80160000], 2)
    comp = 2 / (nclides + 3.0)
    assert mat1.comp[80160000] == pytest.approx(comp)
    assert mat1.comp[922350000] == pytest.approx(comp)
    assert mat1.comp[922380000] == pytest.approx(comp)
    assert mat1.mass == nclides + 3


def test_set_mat_int_2():
    mat = Material(nucvec)
    mat1 = mat.set_mat([922350000, 922380000, 80160000], 2)
    comp = 2.0 / (nclides + 3.0)
    assert mat1.comp[80160000] == pytest.approx(comp)
    assert mat1.comp[922350000] == pytest.approx(comp)
    assert mat1.comp[922380000] == pytest.approx(comp)
    assert mat1.mass == nclides + 3


class TestMassSubMaterialMethods(TestCase):
    "Tests that the Material sub-Material ter member functions work."

    def test_sub_mat_int_1(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_mat([922350000, 922380000, 80160000])
        assert mat1.comp[80160000] == pytest.approx(0.3333333333333)
        assert mat1.comp[922350000] == pytest.approx(0.3333333333333)
        assert mat1.comp[922380000] == pytest.approx(0.3333333333333)
        assert mat1.mass == 3.0

    def test_sub_mat_int_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_mat([922350000, 922380000, 80160000])
        assert mat1.comp[80160000] == pytest.approx(0.3333333333333)
        assert mat1.comp[922350000] == pytest.approx(0.3333333333333)
        assert mat1.comp[922380000] == pytest.approx(0.3333333333333)
        assert mat1.mass == 3.0

    def test_sub_mat_attr_1(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_mat(["U235", "U238", "80160", "H1"])
        assert mat1.comp[10010000] == pytest.approx(0.25)
        assert mat1.comp[80160000] == pytest.approx(0.25)
        assert mat1.comp[922350000] == pytest.approx(0.25)
        assert mat1.comp[922380000] == pytest.approx(0.25)
        assert mat1.mass == 4.0

    def test_sub_mat_attr_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_mat(["U235", "U238", "80160", "H1"])
        assert mat1.comp[10010000] == pytest.approx(0.25)
        assert mat1.comp[80160000] == pytest.approx(0.25)
        assert mat1.comp[922350000] == pytest.approx(0.25)
        assert mat1.comp[922380000] == pytest.approx(0.25)
        assert mat1.mass == 4.0

    def test_sub_elem_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_elem(8)
        assert mat1.comp == {80160000: 1.0}
        assert mat1.mass == 1.0

    def test_sub_elem_2(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_elem(8)
        assert mat1.comp == {80160000: 1.0}
        assert mat1.mass == 1.0

    def test_sub_u_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_u()
        assert mat1.comp == {922350000: 0.5, 922380000: 0.5}
        assert mat1.mass == 2.0

    def test_sub_u_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_u()
        assert mat1.comp == {922350000: 0.5, 922380000: 0.5}
        assert mat1.mass == 2.0

    def test_pu_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_pu()
        assert mat1.comp == {942390000: 0.5, 942410000: 0.5}
        assert mat1.mass == 2.0

    def test_pu_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_pu()
        assert mat1.comp == {942390000: 0.5, 942410000: 0.5}
        assert mat1.mass == 2.0

    def test_lan_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_lan()
        assert mat1.comp == {691690000: 1.0}
        assert mat1.mass == 1.0

    def test_lan_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_lan()
        assert mat1.comp == {691690000: 1.0}
        assert mat1.mass == 1.0

    def test_act_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_act()
        assert mat1.comp[922350000] == 1.0 / 6.0
        assert mat1.comp[922380000] == 1.0 / 6.0
        assert mat1.comp[942390000] == 1.0 / 6.0
        assert mat1.comp[942410000] == 1.0 / 6.0
        assert mat1.comp[952420000] == 1.0 / 6.0
        assert mat1.comp[962440000] == 1.0 / 6.0
        assert mat1.mass == 6.0

    def test_act_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_act()
        assert mat1.comp[922350000] == 1.0 / 6.0
        assert mat1.comp[922380000] == 1.0 / 6.0
        assert mat1.comp[942390000] == 1.0 / 6.0
        assert mat1.comp[942410000] == 1.0 / 6.0
        assert mat1.comp[952420000] == 1.0 / 6.0
        assert mat1.comp[962440000] == 1.0 / 6.0
        assert mat1.mass == 6.0

    def test_tru_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_tru()
        assert mat1.comp[942390000] == 1.0 / 4.0
        assert mat1.comp[942410000] == 1.0 / 4.0
        assert mat1.comp[952420000] == 1.0 / 4.0
        assert mat1.comp[962440000] == 1.0 / 4.0
        assert mat1.mass == 4.0

    def test_tru_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_tru()
        assert mat1.comp[942390000] == 1.0 / 4.0
        assert mat1.comp[942410000] == 1.0 / 4.0
        assert mat1.comp[952420000] == 1.0 / 4.0
        assert mat1.comp[962440000] == 1.0 / 4.0
        assert mat1.mass == 4.0

    def test_ma_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_ma()
        assert mat1.comp[952420000] == 1.0 / 2.0
        assert mat1.comp[962440000] == 1.0 / 2.0
        assert mat1.mass == 2.0

    def test_ma_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_ma()
        assert mat1.comp[952420000] == 1.0 / 2.0
        assert mat1.comp[962440000] == 1.0 / 2.0
        assert mat1.mass == 2.0

    def test_fp_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_fp()
        assert mat1.comp[10010000] == 1.0 / 3.0
        assert mat1.comp[80160000] == 1.0 / 3.0
        assert mat1.comp[691690000] == 1.0 / 3.0
        assert mat1.mass == 3.0

    def test_fp_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_fp()
        assert mat1.comp[10010000] == 1.0 / 3.0
        assert mat1.comp[80160000] == 1.0 / 3.0
        assert mat1.comp[691690000] == 1.0 / 3.0
        assert mat1.mass == 3.0

    def test_sub_range(self):
        mat = Material(nucvec)
        mat1 = mat.sub_range(920000000, 930000000)
        assert mat1.mass == 2.0
        for nuc in mat1:
            assert 920000000 <= nuc < 930000000

        mat1 = mat.sub_range(upper="U238")
        assert mat1.mass == 4.0
        for nuc in mat1:
            assert nuc < 922380000


class TestMaterialOperatorOverloading(TestCase):
    "Tests that the Material operator overloads work."
    u235 = Material({922350000: 1.0})
    u238 = Material({922380000: 1.0})

    def test_add_num(self):
        mat = self.u235 + 30.0
        assert mat.mass == 31.0

    def test_radd_num(self):
        mat = 90 + self.u235
        assert mat.mass == 91.0

    def test_add_mat(self):
        mat = self.u235 + self.u238
        assert mat.comp == {922350000: 0.5, 922380000: 0.5}
        assert mat.mass == 2.0

    def test_mul_num(self):
        mat = self.u235 * 2.0
        assert mat.mass == 2.0

    def test_rmul_num(self):
        mat = 150 * self.u235
        assert mat.mass == 150.0

    def test_div_num(self):
        mat = self.u235 / 10
        assert mat.mass == 0.1


#
# Test atom fraction functions
#
def test_to_atom_frac():
    h2o = {10010000: 0.11191487328888054, 80160000: 0.8880851267111195}
    mat = Material(h2o, atoms_per_molecule=3.0)
    af = mat.to_atom_frac()
    assert mat.atoms_per_molecule == 3.0
    assert af[10010000] == 2.0
    assert af[80160000] == 1.0
    assert mat.molecular_mass() == 18.01056468408


def test_from_atom_frac_meth():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[10010000] == 0.11191487328888054
    assert mat.comp[80160000] == 0.8880851267111195
    assert mat.mass == 18.01056468408
    assert mat.molecular_mass() == 18.01056468408

    h2 = Material({10010000: 1.0}, atoms_per_molecule=2.0)
    h2o = {"O16": 1.0, h2: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[10010000] == 0.11191487328888054
    assert mat.comp[80160000] == 0.8880851267111195
    assert mat.molecular_mass() == 18.01056468408

    mt1 = from_atom_frac({1001: 0.1, 6000: 0.8, 8016: 0.1})
    assert mt1.comp[10010000] == 0.008911815984758811
    assert mt1.comp[60000000] == 0.8496511972150773
    assert mt1.comp[80160000] == 0.14143698680016392
    assert mt1.molecular_mass() == 11.308862682573398

    ihm = Material()
    ihm.from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, "O16": 2.0}
    mat = Material()
    mat.from_atom_frac(uox)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[80160000] == pytest.approx(0.11912625441775175, rel=1e-16)
    assert mat.comp[922350000] == pytest.approx(0.4376375781743612, rel=1e-15)
    assert mat.comp[922380000] == pytest.approx(0.4432361674078869, rel=1e-15)
    assert mat.molecular_mass() / 268.53718683220006 == pytest.approx(1.0, rel=1e-15)


def test_to_atom_dens():
    h2o = {10010000: 0.11191487328888054, 80160000: 0.8880851267111195}
    mat = Material(h2o, density=1.0)
    ad = mat.to_atom_dens()
    assert ad[10010000] / (10.0**22) == pytest.approx(6.68734335169385)
    assert ad[80160000] / (10.0**22) == pytest.approx(3.34367167584692)


#
# Test mapping functions
#


def test_len():
    mat = Material(nucvec)
    assert len(mat) == 9


def test_contains():
    mat = Material(nucvec)
    assert 10010000 in mat
    assert 922350000 in mat
    assert not 92000000 in mat
    pytest.raises(TypeError, lambda: "word" in mat)


def test_getitem_int():
    mat = Material(nucvec)
    assert mat[922350000] == 1.0
    pytest.raises(KeyError, lambda: mat[42])

    mat = Material(leu)
    assert mat[922350000] == 0.04
    assert mat[922380000] == 0.96
    pytest.raises(KeyError, lambda: mat[922340000])


def test_getitem_str():
    mat = Material(nucvec)
    assert mat["U235"] == 1.0
    pytest.raises(RuntimeError, lambda: mat["word"])

    mat = Material(leu)
    assert mat["U235"] == 0.04
    assert mat["U238"] == 0.96
    pytest.raises(KeyError, lambda: mat["U234"])


def test_getitem_slice_int():
    mat = Material(nucvec)
    mat1 = mat[920000000:930000000]
    assert mat1.mass == 2.0
    for nuc in mat1:
        assert 920000000 <= nuc < 930000000

    mat1 = mat[:922380000]
    assert mat1.mass == 4.0
    for nuc in mat1:
        assert nuc < 922380000

    mat1 = mat[922350000:]
    assert mat1.mass == 6.0
    for nuc in mat1:
        assert 922350000 <= nuc


def test_getitem_slice_str():
    mat = Material(nucvec)
    mat1 = mat["U":"NP"]
    assert mat1.mass == 2.0
    for nuc in mat1:
        assert 920000000 <= nuc < 930000000

    mat1 = mat[:"U238"]
    assert mat1.mass == 4.0
    for nuc in mat1:
        assert nuc < 922380000

    mat1 = mat["U235":]
    assert mat1.mass == 6.0
    for nuc in mat1:
        assert 922350000 <= nuc


def test_getitem_sequence():
    mat = Material(nucvec)
    mat1 = mat[922380000, 922350000]
    assert mat1.mass == 2.0
    assert set(mat1) == set([922380000, 922350000])

    mat1 = mat[922380, "H2", "h1"]
    assert mat1.mass == 2.0
    assert set(mat1) == set([922380000, 10010000])


def test_setitem_int():
    mat = Material(nucvec)
    assert mat.mass == 9.0
    assert mat[922350000] == 1.0
    mat[922350000] = 2.0
    assert mat.mass == 10.0
    assert mat[922350000] == 2.0

    mat = Material(leu)
    assert mat.mass == 1.0
    assert mat[922350000] == 0.04
    assert mat[922380000] == 0.96
    pytest.raises(KeyError, lambda: mat[922340000])
    mat[922340000] = 17.0
    assert mat.mass == 18.0
    assert mat[922340000] == 17.0
    assert mat[922350000] == 0.04
    assert mat[922380000] == 0.96


def test_setitem_str():
    mat = Material(nucvec)
    assert mat.mass == 9.0
    assert mat[922350000] == 1.0
    mat["U235"] = 2.0
    assert mat.mass == 10.0
    assert mat[922350000] == 2.0

    mat = Material(leu)
    assert mat.mass == 1.0
    assert mat[922350000] == 0.04
    assert mat[922380000] == 0.96
    pytest.raises(KeyError, lambda: mat[922340000])
    mat["U234"] = 17.0
    assert mat.mass == 18.0
    assert mat[922340000] == 17.0
    assert mat[922350000] == 0.04
    assert mat[922380000] == 0.96


def test_setitem_slice_int():
    mat = Material(nucvec)
    mat_id = id(mat)
    mat[920000000:930000000] = 42
    assert mat_id == id(mat)
    assert mat.mass == 91.0
    assert mat[10010000] == 1.0
    assert mat[922350000] == 42.0
    assert mat[922380000] == 42.0

    mat = Material(nucvec)
    mat[:922380000] = 0.0
    assert mat.mass == 5.0
    for nuc in mat:
        if nuc < 922380000:
            assert mat[nuc] == 0.0
        else:
            assert mat[nuc] == 1.0

    mat = Material(nucvec)
    mat[922350000:] = 2
    assert mat.mass == 15.0
    for nuc in mat:
        if 922350000 <= nuc:
            assert mat[nuc] == 2.0
        else:
            assert mat[nuc] == 1.0


def test_setitem_slice_str():
    mat = Material(nucvec)
    mat["U":"Np"] = 42
    assert mat.mass == 91.0
    assert mat[10010000] == 1.0
    assert mat[922350000] == 42.0
    assert mat[922380000] == 42.0

    mat = Material(nucvec)
    mat[:"U238"] = 0.0
    assert mat.mass == 5.0
    for nuc in mat:
        if nuc < 922380000:
            assert mat[nuc] == 0.0
        else:
            assert mat[nuc] == 1.0

    mat = Material(nucvec)
    mat["U235":] = 2
    assert mat.mass == 15.0
    for nuc in mat:
        if 922350000 <= nuc:
            assert mat[nuc] == 2.0
        else:
            assert mat[nuc] == 1.0


def test_setitem_sequence():
    mat = Material(nucvec)
    mat_id = id(mat)
    mat[922380000, 922350000] = 42
    assert mat_id == id(mat)
    assert mat.mass == 91.0
    assert mat[10010000] == 1.0
    assert mat[922350000] == 42.0
    assert mat[922380000] == 42.0

    mat = Material(nucvec)
    mat[922380000, "H2", "h1"] = 0.0
    assert mat.mass == 7.0
    assert len(mat) == 10
    for nuc in mat:
        if nuc in [10010000, 10020000, 922380000]:
            assert mat[nuc] == 0.0
        else:
            assert mat[nuc] == 1.0

    mat = Material(nucvec)
    mat["U235", "H2", "h1"] = 2
    assert mat.mass == 13.0
    assert len(mat) == 10
    for nuc in mat:
        if nuc in [10010000, 10020000, 922350000]:
            assert mat[nuc] == 2.0
        else:
            assert mat[nuc] == 1.0


def test_delitem_int():
    mat = Material(nucvec)
    assert mat[922350000] == 1.0
    del mat[922350000]
    pytest.raises(KeyError, lambda: mat[922350000])

    mat = Material(leu)
    assert mat[922350000] == 0.04
    del mat[922350000]
    assert mat.mass == 0.96
    assert mat.comp[922380000] == 1.0
    pytest.raises(KeyError, lambda: mat[922350000])


def test_delitem_str():
    mat = Material(nucvec)
    assert mat[922350000] == 1.0
    del mat["U235"]
    pytest.raises(KeyError, lambda: mat[922350000])

    mat = Material(leu)
    assert mat[922350000] == 0.04
    del mat["U235"]
    assert mat.mass == 0.96
    assert mat.comp[922380000] == 1.0
    pytest.raises(KeyError, lambda: mat[922350000])


def test_delitem_slice_int():
    mat = Material(nucvec)
    del mat[920000000:930000000]
    assert mat.mass == 7.0
    assert mat[10010000] == 1.0
    pytest.raises(KeyError, lambda: mat[922350000])
    pytest.raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[:922380000]
    assert mat.mass == 5.0
    for nuc in mat:
        if nuc < 922380000:
            pytest.raises(KeyError, lambda: mat[nuc])
        else:
            assert mat[nuc] == 1.0

    mat = Material(nucvec)
    del mat[922350000:]
    assert mat.mass == 3.0
    for nuc in mat:
        if 922350000 <= nuc:
            pytest.raises(KeyError, lambda: mat[nuc])
        else:
            assert mat[nuc] == 1.0


def test_delitem_slice_str():
    mat = Material(nucvec)
    del mat["U":"Np"]
    assert mat.mass == 7.0
    assert mat[10010000] == 1.0
    pytest.raises(KeyError, lambda: mat[922350000])
    pytest.raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[:"U238"]
    assert mat.mass == 5.0
    for nuc in mat:
        if nuc < 922380000:
            pytest.raises(KeyError, lambda: mat[nuc])
        else:
            assert mat[nuc] == 1.0

    mat = Material(nucvec)
    del mat["U235":]
    assert mat.mass == 3.0
    for nuc in mat:
        if 922350000 <= nuc:
            pytest.raises(KeyError, lambda: mat[nuc])
        else:
            assert mat[nuc] == 1.0


def test_delitem_sequence():
    mat = Material(nucvec)
    del mat[922380000, 922350000]
    assert mat.mass == 7.0
    assert mat[10010000] == 1.0
    pytest.raises(KeyError, lambda: mat[922350000])
    pytest.raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[922380000, "H2", "h1"]
    assert mat.mass == 7.0
    pytest.raises(KeyError, lambda: mat[10010000])
    pytest.raises(KeyError, lambda: mat[10020000])
    pytest.raises(KeyError, lambda: mat[922380000])


def test_iter():
    mat = Material(nucvec)
    for nuc in mat:
        assert mat[nuc] == 1.0

    mat = Material(leu)
    keys = set([922350000, 922380000])
    values = set([0.04, 0.96])
    items = set([(922350000, 0.04), (922380000, 0.96)])

    assert set(mat) == keys
    assert set(mat.values()) == values
    assert set(mat.items()) == items


#
# test material generation functions
#


def test_from_atom_frac_func():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = from_atom_frac(h2o)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[10010000] == 0.11191487328888054
    assert mat.comp[80160000] == 0.8880851267111195
    assert mat.mass == 18.01056468408
    assert mat.molecular_mass() == 18.01056468408

    h2 = Material({10010000: 1.0}, atoms_per_molecule=2.0)
    h2o = {"O16": 1.0, h2: 1.0}
    mat = from_atom_frac(h2o)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[10010000] == 0.11191487328888054
    assert mat.comp[80160000] == 0.8880851267111195
    assert mat.molecular_mass() == 18.01056468408

    ihm = from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, "O16": 2.0}
    mat = from_atom_frac(uox)
    assert mat.atoms_per_molecule == 3.0
    assert mat.comp[80160000] == pytest.approx(0.11912625441775175, rel=1e-16)
    assert mat.comp[922350000] == pytest.approx(0.4376375781743612, rel=1e-15)
    assert mat.comp[922380000] == pytest.approx(0.4432361674078869, rel=1e-15)
    assert mat.molecular_mass() / 268.53718683220006 == pytest.approx(1.0, rel=1e-15)


def test_from_activity_func():
    leu = {922350000: 59953.15151320379, 922380000: 177216.65219208886}
    mat = from_activity(leu)
    assert mat.atoms_per_molecule == -1.0
    assert mat.comp[922350000] == 0.05
    assert mat.comp[922380000] == 0.95
    assert mat.mass == 15.0
    assert mat.molecular_mass() == 237.8986180886295

    hto = {"H": 0.0, "H3": 1.0, "O16": 0.0}
    mat = from_activity(hto)
    assert mat.comp == {10030000: 1.0}
    assert mat.mass == 2.809161396488766e-15

    hto = {"H": 1.0, "H3": 1.0, "O16": 1.0}
    pytest.raises(ValueError, from_activity, hto)


def test_from_hdf5_func_protocol_0():
    mat = from_hdf5("mat.h5", "/mat", protocol=0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat = from_hdf5("mat.h5", "/mat", 0, 0)
    assert mat.mass == 1.0
    assert mat.comp == {922350000: 1.0, 942390000: 0.0}

    mat = from_hdf5("mat.h5", "/mat", 1, 0)
    assert mat.mass == 0.5
    assert mat.comp == {922350000: 0.75, 942390000: 0.25}

    mat = from_hdf5("mat.h5", "/mat", 2, 0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat = from_hdf5("mat.h5", "/mat", -1, 0)
    assert mat.mass == 0.0
    assert mat.comp == {922350000: 0.0, 942390000: 0.0}

    mat = from_hdf5("mat.h5", "/mat", -2, 0)
    assert mat.mass == 0.5
    assert mat.comp == {922350000: 0.75, 942390000: 0.25}

    mat = from_hdf5("mat.h5", "/mat", -3, 0)
    assert mat.mass == 1.0
    assert mat.comp == {922350000: 1.0, 942390000: 0.0}


def test_from_text_func():
    mat = from_text("mat.txt")
    assert mat.comp == {922350000: 0.05, 922380000: 0.95}


def test_map_str_material():
    m = MapStrMaterial()
    m["leu"] = Material(leu)
    m["heu"] = Material({"U238": 0.01, "U235": 0.99}, 42.0)
    assert len(m) == 2
    assert m["leu"].mass == 1.0
    assert m["leu"]["U235"] == 0.04
    assert m["heu"].mass == 42.0
    assert m["heu"]["U238"] == 0.42

    m = MapStrMaterial(
        {"leu": Material(leu), "heu": Material({"U238": 0.01, "U235": 0.99}, 42.0)}
    )
    assert len(m) == 2
    assert m["leu"].mass == 1.0
    assert m["leu"]["U235"] == 0.04
    assert m["heu"].mass == 42.0
    assert m["heu"]["U238"] == 0.42

    n = MapStrMaterial(m, False)
    assert len(n) == 2
    assert n["leu"].mass == 1.0
    assert n["leu"]["U235"] == 0.04
    assert n["heu"].mass == 42.0
    assert n["heu"]["U238"] == 0.42

    # points to the same underlying map
    n["other"] = Material({"PU239": 15.0})
    assert m["other"].mass == 15.0
    assert m["other"]["PU239"] == 15.0

    assert n["leu"].mass == 1.0
    assert n["leu"]["U235"] == 0.04
    assert n["heu"].mass == 42.0
    assert n["heu"]["U238"] == 0.42


def test_metadata():
    mat = Material(leu)
    assert len(mat.metadata) == 0
    mat.metadata["units"] = "kg"
    assert len(mat.metadata) == 1
    assert mat.metadata["units"] == "kg"

    mat.metadata = {"comment": "rawr", "amount": 42.0}
    assert mat.metadata.keys() == ["amount", "comment"]
    assert isinstance(mat.metadata, jsoncpp.Value)

    aview = mat.metadata
    aview["omnomnom"] = [1, 2, 5, 3]
    assert len(mat.metadata) == 3
    assert list(mat.metadata["omnomnom"]) == [1, 2, 5, 3]


#
# Test MultiMaterial
#
def test_multimaterial_mix_composition():
    mat1 = Material(
        nucvec={120240000: 0.3, 300000000: 0.2, 10010000: 0.1}, density=2.71
    )
    mat2 = Material(nucvec={60120000: 0.2, 280640000: 0.5, 10010000: 0.12}, density=8.0)
    mix = MultiMaterial({mat1: 0.5, mat2: 0.21})
    mat3 = mix.mix_by_mass()
    mat4 = mix.mix_by_volume()

    assert mat3.comp[10010000] == 0.16065498683155846
    assert mat3.comp[60120000] == 0.0721401580212985
    assert mat3.comp[120240000] == 0.352112676056338
    assert mat3.comp[280640000] == 0.18035039505324627
    assert mat3.comp[300000000] == 0.2347417840375587

    assert mat4.comp[10010000] == 0.15541581280722197
    assert mat4.comp[60120000] == 0.13501024631333625
    assert mat4.comp[120240000] == 0.2232289950576606
    assert mat4.comp[280640000] == 0.33752561578334067
    assert mat4.comp[300000000] == 0.14881933003844042


def test_multimaterial_mix_density():
    mat1 = Material(nucvec={120240000: 0.3, 300000000: 0.2, 10010000: 0.1}, density=1.0)
    mat2 = Material(nucvec={60120000: 0.2, 280640000: 0.5, 10010000: 0.12}, density=2.0)
    # mixing by hand to get density 1.5 when 50:50 by volume of density 1 & 2
    mix = MultiMaterial({mat1: 0.5, mat2: 0.5})
    mat3 = mix.mix_by_volume()
    # calculated mass fracs by hand, same problem as above stated in mass fraction terms
    # rather than volume fraction terms.
    mix = MultiMaterial({mat1: 1 / 3.0, mat2: 2 / 3.0})
    mat4 = mix.mix_by_mass()

    assert mat3.density == 1.5
    assert mat4.density == 1.5

    assert mat3.density == mat4.density


def test_deepcopy():
    x = Material(
        {"H1": 1.0},
        mass=2.0,
        density=3.0,
        atoms_per_molecule=4.0,
        metadata={"name": "loki"},
    )
    y = deepcopy(x)
    assert x == y
    y.comp[10010000] = 42.0
    y[80160000] = 21.0
    y.density = 65.0
    y.atoms_per_molecule = 28.0
    y.metadata["wakka"] = "jawaka"
    y.mass = 48.0
    y.metadata["name"] = "odin"
    assert x != y
    assert x == Material(
        {"H1": 1.0},
        mass=2.0,
        density=3.0,
        atoms_per_molecule=4.0,
        metadata={"name": "loki"},
    )


def test_openmc():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "LEU",
        },
        density=19.1,
    )

    mass = leu.openmc(indent_lvl=0)
    mass_exp = (
        '<material id="2" name="LEU" >\n'
        '  <density value="19.1" units="g/cc" />\n'
        '  <nuclide name="U235" wo="4.0000e-02" />\n'
        '  <nuclide name="U238" wo="9.6000e-01" />\n'
        "</material>\n"
    )
    assert mass == mass_exp

    atom = leu.openmc(frac_type="atom", indent_lvl=0)
    atom_exp = (
        '<material id="2" name="LEU" >\n'
        '  <density value="19.1" units="g/cc" />\n'
        '  <nuclide name="U235" ao="4.0491e-02" />\n'
        '  <nuclide name="U238" ao="9.5951e-01" />\n'
        "</material>\n"
    )
    assert atom == atom_exp

    # check write/read consistency
    leu.write_hdf5("leu.h5")

    leu_read = Material()
    leu_read.from_hdf5("leu.h5", "/mat_name")

    mass = leu.openmc(indent_lvl=0)
    assert mass == mass_exp

    atom = leu.openmc(frac_type="atom", indent_lvl=0)
    assert atom == atom_exp


def test_openmc_mat0():
    leu = Material(
        nucvec={"U235": 0.04, "U236": 0.0, "U238M": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92236": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "LEU",
        },
        density=19.1,
    )

    mass = leu.openmc(indent_lvl=0)
    mass_exp = (
        '<material id="2" name="LEU" >\n'
        '  <density value="19.1" units="g/cc" />\n'
        '  <nuclide name="U235" wo="4.0000e-02" />\n'
        '  <nuclide name="U238_m1" wo="9.6000e-01" />\n'
        "</material>\n"
    )
    assert mass == mass_exp


def test_openmc_sab():
    leu = Material(
        nucvec={"H1": 0.66, "O16": 0.33},
        metadata={
            "mat_number": 2,
            "sab": "c_H_in_H2O",
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "Water",
        },
        density=1.001,
    )

    mass = leu.openmc(indent_lvl=0)
    mass_exp = (
        '<material id="2" name="Water" >\n'
        '  <density value="1.001" units="g/cc" />\n'
        '  <nuclide name="H1" wo="6.6667e-01" />\n'
        '  <nuclide name="O16" wo="3.3333e-01" />\n'
        '  <sab name="c_H_in_H2O" />\n'
        "</material>\n"
    )
    assert mass == mass_exp


def test_openmc_c():
    csi = Material()
    csi.from_atom_frac({"C": 0.5, "Si": 0.5})
    csi.metadata = {"mat_number": 2, "name": "silicon carbide"}
    csi.density = 3.16

    atom = csi.openmc(frac_type="atom", indent_lvl=0)
    atom_exp = (
        '<material id="2" name="silicon carbide" >\n'
        '  <density value="3.16" units="g/cc" />\n'
        '  <nuclide name="C0" ao="5.0000e-01" />\n'
        '  <nuclide name="Si28" ao="4.6111e-01" />\n'
        '  <nuclide name="Si29" ao="2.3425e-02" />\n'
        '  <nuclide name="Si30" ao="1.5460e-02" />\n'
        "</material>\n"
    )
    assert atom == atom_exp


def test_phits():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "mat_name": "LEU",
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
    )

    mass = leu.phits()
    mass_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c -4.0000e-02\n"
        "     92238.25c -9.6000e-01\n"
    )
    assert mass == mass_exp

    atom = leu.phits(frac_type="atom")
    atom_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c 4.0491e-02\n"
        "     92238.25c 9.5951e-01\n"
    )
    assert atom == atom_exp

    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "mat_name": "LEU",
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
        density=19.1,
    )

    mass = leu.phits()
    mass_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c -7.6400e-01\n"
        "     92238.25c -1.8336e+01\n"
    )
    assert mass == mass_exp

    mass_frac = leu.phits(mult_den=False)
    mass_frac_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c -4.0000e-02\n"
        "     92238.25c -9.6000e-01\n"
    )
    assert mass_frac == mass_frac_exp

    atom = leu.phits(frac_type="atom")
    atom_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c 1.9575e-03\n"
        "     92238.25c 4.6386e-02\n"
    )
    assert atom == atom_exp

    atom_frac = leu.phits(frac_type="atom", mult_den=False)
    atom_frac_exp = (
        "C name: leu\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "M[ 2 ]\n"
        "     92235.15c 4.0491e-02\n"
        "     92238.25c 9.5951e-01\n"
    )
    assert atom_frac == atom_frac_exp


def test_mcnp():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
        density=19.1,
    )

    mass = leu.mcnp()
    mass_exp = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c -7.6400e-01\n"
        "     92238.25c -1.8336e+01\n"
    )
    assert mass == mass_exp

    mass_frac = leu.mcnp(mult_den=False)
    mass_frac_exp = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c -4.0000e-02\n"
        "     92238.25c -9.6000e-01\n"
    )
    assert mass_frac == mass_frac_exp

    atom = leu.mcnp(frac_type="atom")
    atom_exp = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c 1.9575e-03\n"
        "     92238.25c 4.6386e-02\n"
    )
    assert atom == atom_exp

    atom_frac = leu.mcnp(frac_type="atom", mult_den=False)
    atom_frac_exp = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c 4.0491e-02\n"
        "     92238.25c 9.5951e-01\n"
    )
    assert atom_frac == atom_frac_exp


def test_mcnp_mat0():
    leu = Material(
        nucvec={"U235": 0.04, "U236": 0.0, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92236": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
        density=19.1,
    )

    mass = leu.mcnp()
    mass_exp = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c -7.6400e-01\n"
        "     92238.25c -1.8336e+01\n"
    )
    assert mass == mass_exp


def test_uwuw():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
        density=19.1,
    )

    uwuw_name = leu.get_uwuw_name()
    name_exp = "mat:leu/rho:19.1"
    assert uwuw_name == name_exp

    leu2 = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
    )
    uwuw_name = leu2.get_uwuw_name()
    name_exp = "mat:leu"
    assert uwuw_name == name_exp

    no_name = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
        },
        density=19.1,
    )
    uwuw_name = no_name.get_uwuw_name()
    name_exp = ""
    assert uwuw_name == name_exp


def test_alara():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"922350": "15c", "922380": "25c"},
            "name": "LEU",
            "source": "Some URL",
            "comments": "this is a long comment that will definitly go over the 80 character limit, for science",
        },
        density=19.1,
    )
    leu2 = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
        },
        density=19.1,
    )
    leu3 = Material(nucvec={"U235": 0.04, "U238": 0.96})

    written = leu.alara()
    written += leu2.alara()
    written += leu3.alara()

    expected = (
        "# mat number: 2\n"
        "# source: Some URL\n"
        "# comments: this is a long comment that will definitly go over the 80 character\n"
        "#  limit, for science\n"
        "LEU 19.1 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
        "# mat number: 2\n"
        "mat2_rho-19.1 19.1 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
        "mat<mat_num>_rho-<rho> <rho> 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
    )
    assert written == expected


def test_write_openmc():
    if "openmc_mass_fracs.txt" in os.listdir("."):
        os.remove("openmc_mass_fracs.txt")

    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "LEU",
        },
        density=19.1,
    )

    leu.write_openmc("openmc_mass_fracs.txt", indent_lvl=0)
    leu.write_openmc("openmc_mass_fracs.txt", frac_type="atom", indent_lvl=0)

    with open("openmc_mass_fracs.txt") as f:
        written = f.read()
    expected = (
        '<material id="2" name="LEU" >\n'
        '  <density value="19.1" units="g/cc" />\n'
        '  <nuclide name="U235" wo="4.0000e-02" />\n'
        '  <nuclide name="U238" wo="9.6000e-01" />\n'
        "</material>\n"
        '<material id="2" name="LEU" >\n'
        '  <density value="19.1" units="g/cc" />\n'
        '  <nuclide name="U235" ao="4.0491e-02" />\n'
        '  <nuclide name="U238" ao="9.5951e-01" />\n'
        "</material>\n"
    )
    assert written == expected
    os.remove("openmc_mass_fracs.txt")


def test_write_mcnp():
    if "mcnp_mass_fracs.txt" in os.listdir("."):
        os.remove("mcnp_mass_fracs.txt")

    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "source": "Some URL",
            "comments": (
                "this is a long comment that will definitly "
                "go over the 80 character limit, for science"
            ),
            "name": "leu",
        },
        density=19.1,
    )

    leu.write_mcnp("mcnp_mass_fracs.txt")
    leu.write_mcnp("mcnp_mass_fracs.txt", frac_type="atom")

    with open("mcnp_mass_fracs.txt") as f:
        written = f.read()
    expected = (
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c -7.6400e-01\n"
        "     92238.25c -1.8336e+01\n"
        "C name: leu\n"
        "C density = 19.10000\n"
        "C source: Some URL\n"
        "C comments: this is a long comment that will definitly go over the 80 character\n"
        "C  limit, for science\n"
        "m2\n"
        "     92235.15c 1.9575e-03\n"
        "     92238.25c 4.6386e-02\n"
    )
    assert written == expected
    os.remove("mcnp_mass_fracs.txt")


def test_fluka():
    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"92235": "15c", "92238": "25c"},
            "name": "LEU",
            "fluka_name": "URANIUM",
            "fluka_material_index": "35",
            "source": "Some URL",
            "comments": ("Fluka Compound "),
        },
        density=19.1,
    )
    ########################################
    # Part I:  Do not collapse the materials
    id = 25
    matlines = []
    # call fluka() on a material made up of each component
    for key in leu.comp:
        element = Material(nucvec={key: 1})
        matlines.append(element.fluka(id, "atom"))
        id = id + 1
    compound = leu.fluka(id, "atom")
    matlines.append(compound)
    written = "".join(matlines)

    exp = "MATERIAL         92.   235.044        1.       25.                    235-U     \n"
    exp += "MATERIAL         92.   238.051        1.       26.                    238-U     \n"
    exp += "* Fluka Compound \n"
    exp += "MATERIAL          1.        1.      19.1       27.                    URANIUM   \n"
    exp += "COMPOUND   4.000e-02     235-U 9.600e-01     238-U                    URANIUM   \n"

    assert exp == written

    #####################################
    # Part II:  Test a collapsed material
    coll = leu.collapse_elements({920000000})
    coll.metadata["comments"] = "Fluka Element "

    exp = "* Fluka Element"
    exp += " \n"
    exp += (
        "MATERIAL         92.   238.029      19.1       25.                    URANIUM"
    )
    exp += "   \n"

    written = coll.fluka(25, "atom")
    assert exp == written

    ##################################
    # Repeat Part I for mass frac_type
    id = 25
    matlines = []
    # call fluka() on a material made up of each component
    for key in leu.comp:
        element = Material(nucvec={key: 1})
        matlines.append(element.fluka(id, "mass"))
        id = id + 1
    compound = leu.fluka(id, "mass")
    matlines.append(compound)
    written = "".join(matlines)

    exp = "MATERIAL         92.   235.044        1.       25.                    235-U     \n"
    exp += "MATERIAL         92.   238.051        1.       26.                    238-U     \n"
    exp += "* Fluka Compound \n"
    exp += "MATERIAL          1.        1.      19.1       27.                    URANIUM   \n"
    exp += "COMPOUND  -4.000e-02     235-U-9.600e-01     238-U                    URANIUM   \n"
    assert exp == written


def test_fluka_scientific():
    # baseline test for scientific formatting
    mat = Material({"H": 0.1, "O": 0.8, "C": 0.1})
    mat.density = 1.0
    mat.metadata["fluka_name"] = "ORGPOLYM"
    written = mat.fluka(25)

    exp = "MATERIAL          1.        1.        1.       25.                    ORGPOLYM  \n"
    exp += "COMPOUND  -1.000e-01  HYDROGEN-1.000e-01    CARBON-8.000e-01    OXYGENORGPOLYM  \n"
    assert exp == written

    mat = Material({"H": 0.01, "O": 0.8, "C": 0.19})
    mat.density = 1.0
    mat.metadata["fluka_name"] = "ORGPOLYM"
    written = mat.fluka(25)

    exp = "MATERIAL          1.        1.        1.       25.                    ORGPOLYM  \n"
    exp += "COMPOUND  -1.000e-02  HYDROGEN-1.900e-01    CARBON-8.000e-01    OXYGENORGPOLYM  \n"
    assert exp == written


def test_write_alara():
    if "alara.txt" in os.listdir("."):
        os.remove("alara.txt")

    leu = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
            "table_ids": {"922350": "15c", "922380": "25c"},
            "name": "LEU",
            "source": "Some URL",
            "comments": "this is a long comment that will definitly go over the 80 character limit, for science",
        },
        density=19.1,
    )
    leu2 = Material(
        nucvec={"U235": 0.04, "U238": 0.96},
        metadata={
            "mat_number": 2,
        },
        density=19.1,
    )
    leu3 = Material(nucvec={"U235": 0.04, "U238": 0.96})

    leu.write_alara("alara.txt")
    leu2.write_alara("alara.txt")
    leu3.write_alara("alara.txt")

    with open("alara.txt") as f:
        written = f.read()
    expected = (
        "# mat number: 2\n"
        "# source: Some URL\n"
        "# comments: this is a long comment that will definitly go over the 80 character\n"
        "#  limit, for science\n"
        "LEU 19.1 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
        "# mat number: 2\n"
        "mat2_rho-19.1 19.1 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
        "mat<mat_num>_rho-<rho> <rho> 2\n"
        "     u:235 4.0000E+00 92\n"
        "     u:238 9.6000E+01 92\n"
    )
    assert written == expected
    os.remove("alara.txt")


def test_natural_elements():
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    expected_comp = {10000000: 0.1118983878322976, 80000000: 0.8881016121677024}
    for key in expected_comp:
        assert water.comp[key] == pytest.approx(expected_comp[key])


def test_load_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = Material(leu)
    obs = Material()
    json = jsoncpp.Value(
        {
            "mass": 1.0,
            "comp": leu,
            "density": -1.0,
            "metadata": {},
            "atoms_per_molecule": -1.0,
        }
    )
    obs.load_json(json)
    assert exp == obs


def test_dump_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = jsoncpp.Value(
        {
            "mass": 1.0,
            "comp": leu,
            "density": -1.0,
            "metadata": {},
            "atoms_per_molecule": -1.0,
        }
    )
    obs = Material(leu).dump_json()
    assert exp == obs


def test_rw_json():
    filename = "leu.json"
    wmat = Material(leu)
    wmat.write_json(filename)
    rmat = Material()
    rmat.from_json(filename)
    assert wmat == rmat
    os.remove(filename)


def test_material_gammas():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert_array_equal(
        mat.gammas(),
        [
            (19.55, np.nan),
            (31.6, 2.147613653785435e-22),
            (34.7, 4.674217952356534e-22),
            (41.4, 3.7899064478566484e-22),
            (41.96, 7.579812895713297e-22),
            (51.21, 4.29522730757087e-22),
            (54.1, 1.2633021492855498e-23),
            (54.25, 3.7899064478566484e-22),
            (60.5, np.nan),
            (64.45, np.nan),
            (72.7, 1.5159625791426594e-21),
            (74.94, 6.442840961356303e-22),
            (76.2, np.nan),
            (94.0, np.nan),
            (95.7, np.nan),
            (96.09, 1.1496049558498501e-21),
            (109.19, 2.0970815678140124e-20),
            (115.45, 3.7899064478566484e-22),
            (120.35, 3.284585588142429e-22),
            (136.55, 1.5159625791426599e-22),
            (140.76, 2.5266042985710997e-21),
            (142.4, 6.316510746427748e-23),
            (143.76, 1.3845791556169625e-19),
            (147.0, np.nan),
            (150.93, 1.1369719343569947e-21),
            (163.356, 6.417574918370591e-20),
            (173.3, 7.579812895713299e-23),
            (182.1, np.nan),
            (182.62, 4.9268783822136435e-21),
            (185.715, 7.200822250927634e-19),
            (194.94, 7.958803540498964e-21),
            (198.9, 4.547887737427978e-22),
            (202.12, 1.3643663212283936e-20),
            (205.316, 6.341776789413459e-20),
            (215.28, 3.6635762329280945e-22),
            (221.386, 1.4906965361569485e-21),
            (228.78, 8.843115044998849e-23),
            (233.5, 4.800548167285089e-22),
            (240.88, 9.348435904713068e-22),
            (246.83, 6.948161821070524e-22),
            (266.45, 7.579812895713299e-23),
            (275.35, 6.442840961356303e-22),
            (275.49, 4.0425668777137595e-22),
            (281.42, 7.579812895713299e-23),
            (282.92, 7.579812895713299e-23),
            (289.56, 8.843115044998849e-23),
            (291.2, np.nan),
            (291.65, 5.053208597142199e-22),
            (301.7, 6.316510746427748e-23),
            (317.1, 1.2633021492855498e-23),
            (343.5, 3.7899064478566497e-23),
            (345.92, 5.053208597142199e-22),
            (356.03, 6.316510746427748e-23),
            (387.84, 5.053208597142199e-22),
            (410.29, 3.7899064478566497e-23),
            (428.7, np.nan),
            (448.4, 1.2633021492855498e-23),
            (49.55, 3.018820986575914e-21),
            (113.5, 4.811245947355363e-22),
        ],
    )


def test_material_xrays():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert mat.xrays() == [
        (93.35, 7.135646069465136e-18),
        (89.953, 4.4112564001433475e-18),
        (105.0, 3.394789326064896e-18),
        (13.0, 2.9652250902501514e-17),
        (93.35, 5.2767097753586244e-21),
        (89.953, 3.2620619831267015e-21),
        (105.0, 2.510398896994686e-21),
        (13.0, 3.4433858446293e-17),
    ]


def test_material_photons():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert_array_equal(
        mat.photons(),
        [
            (19.55, np.nan),
            (31.6, 2.147613653785435e-22),
            (34.7, 4.674217952356534e-22),
            (41.4, 3.7899064478566484e-22),
            (41.96, 7.579812895713297e-22),
            (51.21, 4.29522730757087e-22),
            (54.1, 1.2633021492855498e-23),
            (54.25, 3.7899064478566484e-22),
            (60.5, np.nan),
            (64.45, np.nan),
            (72.7, 1.5159625791426594e-21),
            (74.94, 6.442840961356303e-22),
            (76.2, np.nan),
            (94.0, np.nan),
            (95.7, np.nan),
            (96.09, 1.1496049558498501e-21),
            (109.19, 2.0970815678140124e-20),
            (115.45, 3.7899064478566484e-22),
            (120.35, 3.284585588142429e-22),
            (136.55, 1.5159625791426599e-22),
            (140.76, 2.5266042985710997e-21),
            (142.4, 6.316510746427748e-23),
            (143.76, 1.3845791556169625e-19),
            (147.0, np.nan),
            (150.93, 1.1369719343569947e-21),
            (163.356, 6.417574918370591e-20),
            (173.3, 7.579812895713299e-23),
            (182.1, np.nan),
            (182.62, 4.9268783822136435e-21),
            (185.715, 7.200822250927634e-19),
            (194.94, 7.958803540498964e-21),
            (198.9, 4.547887737427978e-22),
            (202.12, 1.3643663212283936e-20),
            (205.316, 6.341776789413459e-20),
            (215.28, 3.6635762329280945e-22),
            (221.386, 1.4906965361569485e-21),
            (228.78, 8.843115044998849e-23),
            (233.5, 4.800548167285089e-22),
            (240.88, 9.348435904713068e-22),
            (246.83, 6.948161821070524e-22),
            (266.45, 7.579812895713299e-23),
            (275.35, 6.442840961356303e-22),
            (275.49, 4.0425668777137595e-22),
            (281.42, 7.579812895713299e-23),
            (282.92, 7.579812895713299e-23),
            (289.56, 8.843115044998849e-23),
            (291.2, np.nan),
            (291.65, 5.053208597142199e-22),
            (301.7, 6.316510746427748e-23),
            (317.1, 1.2633021492855498e-23),
            (343.5, 3.7899064478566497e-23),
            (345.92, 5.053208597142199e-22),
            (356.03, 6.316510746427748e-23),
            (387.84, 5.053208597142199e-22),
            (410.29, 3.7899064478566497e-23),
            (428.7, np.nan),
            (448.4, 1.2633021492855498e-23),
            (49.55, 3.018820986575914e-21),
            (113.5, 4.811245947355363e-22),
            (93.35, 7.135646069465136e-18),
            (89.953, 4.4112564001433475e-18),
            (105.0, 3.394789326064896e-18),
            (13.0, 2.9652250902501514e-17),
            (93.35, 5.2767097753586244e-21),
            (89.953, 3.2620619831267015e-21),
            (105.0, 2.510398896994686e-21),
            (13.0, 3.4433858446293e-17),
        ],
    )
    assert mat.photons(True) == [
        (31.6, 2.681359313180288e-06),
        (34.7, 5.835899681627685e-06),
        (41.4, 4.731810552671095e-06),
        (41.96, 9.46362110534219e-06),
        (51.21, 5.362718626360576e-06),
        (54.1, 1.577270184223699e-07),
        (54.25, 4.731810552671095e-06),
        (72.7, 1.892724221068438e-05),
        (74.94, 8.044077939540861e-06),
        (96.09, 1.4353158676435657e-05),
        (109.19, 0.00026182685058113395),
        (115.45, 4.731810552671095e-06),
        (120.35, 4.100902478981616e-06),
        (136.55, 1.8927242210684386e-06),
        (140.76, 3.1545403684473975e-05),
        (142.4, 7.886350921118493e-07),
        (143.76, 0.0017286881219091737),
        (150.93, 1.4195431658013287e-05),
        (163.356, 0.0008012532535856387),
        (173.3, 9.463621105342193e-07),
        (182.62, 6.151353718472424e-05),
        (185.715, 0.008990440050075084),
        (194.94, 9.936802160609301e-05),
        (198.9, 5.678172663205315e-06),
        (202.12, 0.00017034517989615943),
        (205.316, 0.0007917896324802966),
        (215.28, 4.574083534248726e-06),
        (221.386, 1.8611788173839642e-05),
        (228.78, 1.1040891289565891e-06),
        (233.5, 5.993626700050055e-06),
        (240.88, 1.167179936325537e-05),
        (246.83, 8.674986013230342e-06),
        (266.45, 9.463621105342193e-07),
        (275.35, 8.044077939540861e-06),
        (275.49, 5.047264589515836e-06),
        (281.42, 9.463621105342193e-07),
        (282.92, 9.463621105342193e-07),
        (289.56, 1.1040891289565891e-06),
        (291.65, 6.3090807368947945e-06),
        (301.7, 7.886350921118493e-07),
        (317.1, 1.577270184223699e-07),
        (343.5, 4.7318105526710964e-07),
        (345.92, 6.3090807368947945e-06),
        (356.03, 7.886350921118493e-07),
        (387.84, 6.3090807368947945e-06),
        (410.29, 4.7318105526710964e-07),
        (448.4, 1.577270184223699e-07),
        (49.55, 3.7690874952819364e-05),
        (113.5, 6.006983195605586e-06),
        (93.35, 0.08909065655358435),
        (89.953, 0.055075843881425855),
        (105.0, 0.04238495112789303),
        (13.0, 0.37021714298583336),
        (93.35, 6.58812858363439e-05),
        (89.953, 4.0727810904027796e-05),
        (105.0, 3.134307444166928e-05),
        (13.0, 0.42991693068703124),
    ]


def test_decay_h3():
    mat = Material({"H3": 1.0})
    obs = mat.decay(data.half_life("H3"))
    obs = obs.to_atom_frac()
    assert 2 == len(obs)
    assert 0.5 == pytest.approx(obs[nucname.id("H3")])
    assert 0.5 == pytest.approx(obs[nucname.id("He3")])


def test_decay_u235_h3():
    mat = Material({"U235": 1.0, "H3": 1.0})
    obs = mat.decay(365.25 * 24.0 * 3600.0)
    if len(obs) < 4:
        # full decay is not installed
        pytest.skip("full decay is not installed")
    exp = Material(
        {
            10030000: 0.472645829730143,
            20030000: 0.027354079574566214,
            812050000: 4.508489735920749e-26,
            812070000: 9.08083992078195e-22,
            822090000: 5.318134090224469e-29,
            822110000: 1.2900842350157843e-20,
            832110000: 8.383482900183342e-22,
            832150000: 4.5950843264546854e-27,
            842110000: 6.3727159025095244e-27,
            842150000: 1.086256809210682e-26,
            852150000: 1.2293001051164733e-33,
            852190000: 4.0236546470124826e-28,
            862190000: 3.37645671770566e-24,
            872230000: 1.5415521899415466e-22,
            882230000: 4.443454725452303e-18,
            882270000: 2.2016578572254302e-26,
            892270000: 4.766912006620537e-15,
            902270000: 1.1717834795114714e-17,
            902310000: 2.0323933624775356e-12,
            912310000: 4.838922882478526e-10,
            922350000: 0.5000000898212907,
        },
        1.9999996387457337,
        -1.0,
        1.000000000011328,
        {},
    )
    exp_sf = Material(
        {
            10030000: 0.472645829399245,
            20030000: 0.027354080100513937,
            360830000: 6.510656251575286e-23,
            420950000: 9.099745814288571e-22,
            430990000: 8.907321880441823e-22,
            440990000: 1.5296296426796811e-27,
            441010000: 7.6508902214678565e-22,
            441030000: 7.087470445136458e-23,
            441060000: 4.55321359080199e-23,
            451030000: 3.870394241355748e-22,
            451050000: 8.731929824567198e-25,
            451060000: 4.262140221110404e-29,
            461050000: 1.4917209641282872e-22,
            461060000: 1.7247438294978485e-23,
            461070000: 2.340390792626567e-23,
            461080000: 8.797869681724177e-24,
            471070000: 2.8703829980977768e-30,
            471090000: 5.137807148513997e-24,
            481130000: 2.37290541344998e-24,
            491130000: 3.7893178733956063e-31,
            491150000: 2.086255676429189e-24,
            511250000: 5.618864598703359e-24,
            521250000: 7.354580944203714e-25,
            521270000: 8.028983133775622e-27,
            531270000: 2.4596893809913248e-23,
            531350000: 1.3473046087881196e-24,
            541310000: 5.560602432935402e-22,
            541340000: 1.5401102695220697e-21,
            541350000: 1.948090703621348e-24,
            541360000: 1.21369662116544e-21,
            551330000: 1.3048948184817563e-21,
            551340000: 1.2768243455987543e-27,
            551350000: 1.2918499510197182e-21,
            551370000: 1.2586108509091848e-21,
            561340000: 2.2623061311731704e-28,
            561350000: 7.027198085245094e-28,
            561370000: 1.4556789341729914e-23,
            601430000: 1.2493083964840656e-21,
            601440000: 8.996484819074119e-40,
            601450000: 8.370089059982299e-22,
            611470000: 2.251790189422784e-22,
            611480000: 2.1211513872058244e-32,
            611490000: 2.0629596934137826e-24,
            621470000: 3.1056525228624105e-23,
            621480000: 9.791735952953553e-31,
            621490000: 2.340845364493762e-22,
            621500000: 6.574595962369667e-27,
            621510000: 9.227361596246668e-23,
            621520000: 5.946769139576938e-23,
            621540000: 3.0493735453482623e-34,
            631510000: 3.5578068109202716e-25,
            631520000: 3.8851842932960325e-32,
            631530000: 3.5524825148041183e-23,
            631540000: 4.148574990903499e-29,
            631550000: 6.823053915339241e-24,
            641520000: 2.8050385431903425e-34,
            641540000: 1.696263607824536e-30,
            641550000: 5.096040322905076e-25,
            641560000: 3.420111353119681e-24,
            641570000: 1.4355487077868398e-24,
            641580000: 7.755331100244439e-25,
            661600000: 6.578079403862082e-30,
            661610000: 2.1296608420806177e-26,
            661620000: 3.964908141740868e-27,
            661630000: 1.474404469483802e-27,
            661640000: 4.709965621984415e-28,
            671650000: 2.366463632641999e-28,
            681660000: 9.748274822324674e-30,
            681670000: 2.6380928977497267e-34,
            812050000: 4.578934888050694e-26,
            812070000: 1.4679523752742383e-21,
            822090000: 5.375945096401122e-29,
            822110000: 1.1356771668957228e-20,
            832110000: 6.732148973436892e-22,
            832150000: 2.4040023541449982e-27,
            842110000: 7.467000602480984e-27,
            842150000: 9.517786596239864e-27,
            852150000: 1.2291314459098075e-33,
            852190000: 3.2131021040718674e-28,
            862190000: 2.1557203796457774e-23,
            872230000: 1.2854004582605104e-22,
            882230000: 5.474385066152288e-18,
            882270000: 2.201657636910653e-26,
            892270000: 4.936006649788287e-15,
            902270000: 9.88011346084414e-18,
            902310000: 2.0323935890826227e-12,
            912310000: 4.818599281474088e-10,
            922350000: 0.5000000900163438,
        },
        1.9999996379655214,
        -1.0,
        1.0,
        {},
    )
    if len(exp) == len(obs):
        # no spontaneous fission
        assert_mat_almost_equal(exp, obs)
    elif len(exp_sf) == len(obs):
        # with spontaneous fission
        assert_mat_almost_equal(exp_sf, obs)
    else:
        assert False, "Observed material does not have correct length"


def test_cram_h3():
    mat = Material({"H3": 1.0})
    A = -cram.DECAY_MATRIX * data.half_life("H3")
    obs = mat.cram(A, order=16)
    obs = obs.to_atom_frac()
    assert 2 == len(obs)
    assert 0.5 == pytest.approx(obs[nucname.id("H3")])
    assert 0.5 == pytest.approx(obs[nucname.id("He3")])
