"""Material tests"""
import os
from copy import deepcopy
import warnings

from unittest import TestCase
import nose
from nose.plugins.skip import SkipTest
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import nuc_data
from pyne.material import Material, from_atom_frac, from_hdf5, from_text, \
    MapStrMaterial, MultiMaterial, MaterialLibrary
from pyne import jsoncpp
from pyne import data
from pyne import nucname
from pyne import utils
import numpy as np
from numpy.testing import assert_array_equal
import tables as tb

if utils.use_warnings():
    utils.toggle_warnings()

nclides = 9
nucvec = {10010000:  1.0,
          80160000:  1.0,
          691690000: 1.0,
          922350000: 1.0,
          922380000: 1.0,
          942390000: 1.0,
          942410000: 1.0,
          952420000: 1.0,
          962440000: 1.0,
          }

leu = {922380000: 0.96, 922350000: 0.04}


def assert_mat_almost_equal(first, second, places=7):
    assert_almost_equal(first.mass, second.mass, places=places)
    assert_almost_equal(first.density, second.density, places=places)
    assert_almost_equal(first.atoms_per_molecule, second.atoms_per_molecule, places=places)
    assert_equal(first.metadata, second.metadata)
    nucs = set(second.comp)
    assert_equal(set(first.comp), nucs)
    for nuc in nucs:
        assert_almost_equal(first.comp[nuc], second.comp[nuc], places=places)


def make_mat_txt():
    """Helper for mat.txt"""
    with open('mat.txt', 'w') as f:
        f.write('U235  0.05\nU238  0.95')


def make_mat_h5():
    """Helper for mat.h5"""
    f = tb.openFile("mat.h5", "w")
    f.createGroup("/", "mat", "Mass Material Test")
    f.createArray("/mat", "Mass",  np.array([1.0, 0.5,  0.0]), "Mass Test")
    f.createArray("/mat", "U235",  np.array([1.0, 0.75, 0.0]), "U235 Test")
    f.createArray("/mat", "PU239", np.array([0.0, 0.25, 0.0]), "PU239 Test")
    f.close()


###############################################################################

def test_mat1():
    mat = Material("mat.txt")
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.mass, 1.0)

def test_mat2():
    mat = Material("mat.txt", 42)
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.mass, 42.0)

def test_mat3():
    mat = Material("mat.txt", -42)
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.mass, 1.0)

def test_mat4():
    mat = Material({922350000: 0.05, 922380000: 0.95}, 15, metadata={'units': 'kg'})
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.mass, 15.0)
    assert_equal(mat.metadata['units'], 'kg')

def test_from_text():
    mat = Material(metadata={'units': 'kg'})
    mat.from_text("mat.txt")
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.metadata['units'], 'kg')


def test_write_text():
    if 'leu.txt' in os.listdir('.'):
        os.remove('leu.txt')

    leu = Material({'U235': 0.04, 'U238': 0.96}, 42.0, 1.0, 1.0)
    leu.write_text('leu.txt')

    with open('leu.txt') as f:
        written = f.read()
    expected = ("Mass    42\n"
                "Density 1\n"
                "APerM   1\n"
                "U235    0.04\n"
                "U238    0.96\n")
    assert_equal(written, expected)

    read_leu = from_text('leu.txt')
    assert_equal(leu.mass, read_leu.mass)
    assert_equal(leu.atoms_per_molecule, read_leu.atoms_per_molecule)
    assert_equal(leu.comp, read_leu.comp)

    os.remove('leu.txt')



def test_from_hdf5_protocol_0():
    mat = Material()
    mat.from_hdf5("mat.h5", "/mat", protocol=0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat.from_hdf5("mat.h5", "/mat", 0, 0)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat.comp, {922350000: 1.0, 942390000: 0.0})

    mat.from_hdf5("mat.h5", "/mat", 1, 0)
    assert_equal(mat.mass, 0.5)
    assert_equal(mat.comp, {922350000: 0.75, 942390000: 0.25})

    mat.from_hdf5("mat.h5", "/mat", 2, 0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat.from_hdf5("mat.h5", "/mat", -1, 0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat.from_hdf5("mat.h5", "/mat", -2, 0)
    assert_equal(mat.mass, 0.5)
    assert_equal(mat.comp, {922350000: 0.75, 942390000: 0.25})

    mat.from_hdf5("mat.h5", "/mat", -3, 0)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat.comp, {922350000: 1.0, 942390000: 0.0})


def test_hdf5_protocol_1():
    if 'proto1.h5' in os.listdir('.'):
        os.remove('proto1.h5')

    # Test material writing
    leu = Material({'U235': 0.04, 'U238': 0.96}, 4.2, 2.72, 1.0)
    leu.metadata['comment'] = 'first light'
    leu.write_hdf5('proto1.h5', chunksize=10)

    for i in range(2, 11):
        leu = Material({'U235': 0.04, 'U238': 0.96}, i*4.2, 2.72, 1.0*i)
        leu.metadata['comment'] = 'fire in the disco - {0}'.format(i)
        leu.write_hdf5('proto1.h5')

    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5('proto1.h5', '/material', -3, 1)
    assert_equal(m.density, 2.72)
    assert_equal(m.atoms_per_molecule, 8.0)
    assert_equal(m.mass, 33.6)
    assert_equal(m.comp, {922350000: 0.04, 922380000: 0.96})
    assert_equal(m.metadata['comment'], 'fire in the disco - 8')

    m = from_hdf5('proto1.h5', '/material', 3, 1)
    assert_equal(m.density, 2.72)
    assert_equal(m.atoms_per_molecule, 4.0)
    assert_equal(m.mass, 16.8)
    assert_equal(m.comp, {922350000: 0.04, 922380000: 0.96})
    assert_equal(m.metadata['comment'], 'fire in the disco - 4')

    os.remove('proto1.h5')


class TestMaterialMethods(TestCase):
    "Tests that the Material member functions work."

    def test_normalize(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        mat.normalize()
        assert_equal(mat.mass, 1.0)


    def test_mult_by_mass(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        nucvec = mat.mult_by_mass()
        assert_equal(nucvec, {922350000: 0.75, 922380000: 14.25})


    def test_activity(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        obs = mat.activity()
        exp = {922350000: 59953.15101810882, 922380000: 177216.65112976026}       
        assert_equal(set(obs), set(exp))
        assert_equal(set(obs.values()), set(exp.values()))


    def test_decay_heat(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        obs = mat.decay_heat()
        exp = {922350000: 4.48963565256e-14, 922380000: 1.2123912039e-13}
        assert_equal(set(obs), set(exp))
        for key in exp:
            assert_almost_equal(obs[key], exp[key])


    def test_dose_per_g(self):
        mat = Material({922350000: 0.05, 922380000: 0.95}, 15)
        # testing for default source
        obs1 = mat.dose_per_g("ext_air")
        exp1 = {922350000: 1.11264406283e-14, 922380000: 5.01315571163e-15}
        assert_equal(set(obs1), set(exp1))
        for key in exp1:
            assert_almost_equal(obs1[key], exp1[key])
        # testing for non-default source
        obs2 = mat.dose_per_g("ingest", 1)
        exp2 = {922350000: 27.1139475504, 922380000: 77.5921552819}
        assert_equal(set(obs2), set(exp2))
        for key in exp2:
            assert_almost_equal(obs2[key], exp2[key])


    def test_molecular_mass(self):
        mat_empty = Material({})
        assert_equal(mat_empty.molecular_mass(), 0.0)

        mat_u238 = Material({922380000: 1.0})
        mw_u238 = mat_u238.molecular_mass()
        try:
            assert_almost_equal(mw_u238, 238.050788423)
        except AssertionError:
            assert_almost_equal(mw_u238, 238.0)

        mat_mixed = Material({922350000: 0.5, 922380000: 0.5})
        mw_mixed = mat_mixed.molecular_mass()
        try:
            assert_almost_equal(mw_mixed/236.547360417, 1.0, 4)
        except AssertionError:
            assert_almost_equal(mw_mixed/236.5, 1.0, 4)


def test_expand_elements1():
    natmat = Material({'C': 1.0, 902320000: 0.5, 'PU': 4.0, 'U': 3.0},
                       metadata={'y': 1.0})
    expmat = natmat.expand_elements()
    assert_true(60120000 in expmat.comp)
    assert_false(60000000 in expmat.comp)
    assert_true(natmat.metadata == expmat.metadata)
    assert_false(natmat.metadata is expmat.metadata)

def test_expand_elements2():
    """Inspired by #86"""
    natmat = Material({'C': 1.0})
    expmat = natmat.expand_elements()
    afrac = expmat.to_atom_frac()
    assert_almost_equal(data.natural_abund(60120000), afrac[60120000])
    assert_almost_equal(data.natural_abund(60130000), afrac[60130000])

def test_collapse_elements1():
    """ Very simple test to combine nucids"""
    nucvec = {10010000:  1.0,
      80160000: 1.0,
      80160001: 1.0,
      691690000: 1.0,
      922350000: 1.0,
      922380000: 1.0,
      942390000: 1.0,
      952420000: 1.0,
      962440000: 1.0 }

    exception_ids = {nucname.id(1001),
                     nucname.id("U-235"),
                     nucname.id("U-238"),
                     nucname.id("Pu-239"),
                     nucname.id("Pu-241"),
                     }

    mat  = Material(nucvec)

    print("Original")
    print(mat)

    cmat = mat.collapse_elements(exception_ids)
    print("Collapsed")
    print(cmat)

    assert_equal(cmat.comp[80000000],  mat.comp[80160000] + mat.comp[80160001])
    assert_equal(cmat.comp[922350000], mat.comp[922350000])
    assert_equal(cmat.comp[942390000], mat.comp[942390000])
    assert_equal(cmat.comp[950000000], mat.comp[952420000])
    assert_equal(cmat.comp[960000000], mat.comp[952420000])


def test_mass_density():
    ethanol = from_atom_frac({'C':2, 'H':6, 'O':1})
    atom_density_ethanol = 9.282542841E22  # atom density not molecule density
    mass_density = ethanol.mass_density(atom_density_ethanol)
    expected_mass_density = 0.78900
    assert_almost_equal(mass_density, expected_mass_density, 4)


def test_number_density():
    ethanol = from_atom_frac({'C':2, 'H':6, 'O':1}, density=0.78900)
    obs = ethanol.number_density()
    exp = 9.2825E22
    assert_almost_equal(obs / exp, 1.0, 4)


def test_set_mat_int_1():
    mat = Material(nucvec, -1)
    mat1 = mat.set_mat([922350000, 922380000, 80160000], 2)
    comp = 2 / (nclides + 3.)
    assert_almost_equal(mat1.comp[80160000],  comp)
    assert_almost_equal(mat1.comp[922350000], comp)
    assert_almost_equal(mat1.comp[922380000], comp)
    assert_equal(mat1.mass, nclides + 3)


def test_set_mat_int_2():
    mat = Material(nucvec)
    mat1 = mat.set_mat([922350000, 922380000, 80160000], 2)
    comp = 2. / (nclides + 3.)
    assert_almost_equal(mat1.comp[80160000],  comp)
    assert_almost_equal(mat1.comp[922350000], comp)
    assert_almost_equal(mat1.comp[922380000], comp)
    assert_equal(mat1.mass, nclides + 3)


class TestMassSubMaterialMethods(TestCase):
    "Tests that the Material sub-Material ter member functions work."

    def test_sub_mat_int_1(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_mat([922350000, 922380000, 80160000])
        assert_almost_equal(mat1.comp[80160000],  0.3333333333333)
        assert_almost_equal(mat1.comp[922350000], 0.3333333333333)
        assert_almost_equal(mat1.comp[922380000], 0.3333333333333)
        assert_equal(mat1.mass, 3.0)

    def test_sub_mat_int_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_mat([922350000, 922380000, 80160000])
        assert_almost_equal(mat1.comp[80160000],  0.3333333333333)
        assert_almost_equal(mat1.comp[922350000], 0.3333333333333)
        assert_almost_equal(mat1.comp[922380000], 0.3333333333333)
        assert_equal(mat1.mass, 3.0)

    def test_sub_mat_attr_1(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_mat(["U235", "U238", "80160", "H1"])
        assert_almost_equal(mat1.comp[10010000],  0.25)
        assert_almost_equal(mat1.comp[80160000],  0.25)
        assert_almost_equal(mat1.comp[922350000], 0.25)
        assert_almost_equal(mat1.comp[922380000], 0.25)
        assert_equal(mat1.mass, 4.0)

    def test_sub_mat_attr_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_mat(["U235", "U238", "80160", "H1"])
        assert_almost_equal(mat1.comp[10010000],  0.25)
        assert_almost_equal(mat1.comp[80160000],  0.25)
        assert_almost_equal(mat1.comp[922350000], 0.25)
        assert_almost_equal(mat1.comp[922380000], 0.25)
        assert_equal(mat1.mass, 4.0)

    def test_sub_elem_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_elem(8)
        assert_equal(mat1.comp, {80160000: 1.0})
        assert_equal(mat1.mass, 1.0)

    def test_sub_elem_2(self):
        mat = Material(nucvec, -1)
        mat1 = mat.sub_elem(8)
        assert_equal(mat1.comp, {80160000: 1.0})
        assert_equal(mat1.mass, 1.0)

    def test_sub_u_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_u()
        assert_equal(mat1.comp, {922350000: 0.5, 922380000: 0.5})
        assert_equal(mat1.mass, 2.0)

    def test_sub_u_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_u()
        assert_equal(mat1.comp, {922350000: 0.5, 922380000: 0.5})
        assert_equal(mat1.mass, 2.0)

    def test_pu_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_pu()
        assert_equal(mat1.comp, {942390000: 0.5, 942410000: 0.5})
        assert_equal(mat1.mass, 2.0)

    def test_pu_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_pu()
        assert_equal(mat1.comp, {942390000: 0.5, 942410000: 0.5})
        assert_equal(mat1.mass, 2.0)

    def test_lan_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_lan()
        assert_equal(mat1.comp, {691690000: 1.0})
        assert_equal(mat1.mass, 1.0)

    def test_lan_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_lan()
        assert_equal(mat1.comp, {691690000: 1.0})
        assert_equal(mat1.mass, 1.0)

    def test_act_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_act()
        assert_equal(mat1.comp[922350000], 1.0/6.0)
        assert_equal(mat1.comp[922380000], 1.0/6.0)
        assert_equal(mat1.comp[942390000], 1.0/6.0)
        assert_equal(mat1.comp[942410000], 1.0/6.0)
        assert_equal(mat1.comp[952420000], 1.0/6.0)
        assert_equal(mat1.comp[962440000], 1.0/6.0)
        assert_equal(mat1.mass, 6.0)

    def test_act_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_act()
        assert_equal(mat1.comp[922350000], 1.0/6.0)
        assert_equal(mat1.comp[922380000], 1.0/6.0)
        assert_equal(mat1.comp[942390000], 1.0/6.0)
        assert_equal(mat1.comp[942410000], 1.0/6.0)
        assert_equal(mat1.comp[952420000], 1.0/6.0)
        assert_equal(mat1.comp[962440000], 1.0/6.0)
        assert_equal(mat1.mass, 6.0)

    def test_tru_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_tru()
        assert_equal(mat1.comp[942390000], 1.0/4.0)
        assert_equal(mat1.comp[942410000], 1.0/4.0)
        assert_equal(mat1.comp[952420000], 1.0/4.0)
        assert_equal(mat1.comp[962440000], 1.0/4.0)
        assert_equal(mat1.mass, 4.0)

    def test_tru_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_tru()
        assert_equal(mat1.comp[942390000], 1.0/4.0)
        assert_equal(mat1.comp[942410000], 1.0/4.0)
        assert_equal(mat1.comp[952420000], 1.0/4.0)
        assert_equal(mat1.comp[962440000], 1.0/4.0)
        assert_equal(mat1.mass, 4.0)

    def test_ma_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_ma()
        assert_equal(mat1.comp[952420000], 1.0/2.0)
        assert_equal(mat1.comp[962440000], 1.0/2.0)
        assert_equal(mat1.mass, 2.0)

    def test_ma_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_ma()
        assert_equal(mat1.comp[952420000], 1.0/2.0)
        assert_equal(mat1.comp[962440000], 1.0/2.0)
        assert_equal(mat1.mass, 2.0)

    def test_fp_1(self):
        mat = Material(nucvec)
        mat1 = mat.sub_fp()
        assert_equal(mat1.comp[10010000],  1.0/3.0)
        assert_equal(mat1.comp[80160000],  1.0/3.0)
        assert_equal(mat1.comp[691690000], 1.0/3.0)
        assert_equal(mat1.mass, 3.0)

    def test_fp_2(self):
        mat = Material(nucvec)
        mat1 = mat.sub_fp()
        assert_equal(mat1.comp[10010000],  1.0/3.0)
        assert_equal(mat1.comp[80160000],  1.0/3.0)
        assert_equal(mat1.comp[691690000], 1.0/3.0)
        assert_equal(mat1.mass, 3.0)

    def test_sub_range(self):
        mat = Material(nucvec)
        mat1 = mat.sub_range(920000000, 930000000)
        assert_equal(mat1.mass, 2.0)
        for nuc in mat1:
            assert_true(920000000 <= nuc < 930000000)

        mat1 = mat.sub_range(upper="U238")
        assert_equal(mat1.mass, 4.0)
        for nuc in mat1:
            assert_true(nuc < 922380000)


class TestMaterialOperatorOverloading(TestCase):
    "Tests that the Material operator overloads work."
    u235 = Material({922350000: 1.0})
    u238 = Material({922380000: 1.0})

    def test_add_num(self):
        mat = self.u235 + 30.0
        assert_equal(mat.mass, 31.0)

    def test_radd_num(self):
        mat = 90 + self.u235
        assert_equal(mat.mass, 91.0)

    def test_add_mat(self):
        mat = self.u235 + self.u238
        assert_equal(mat.comp, {922350000: 0.5, 922380000: 0.5})
        assert_equal(mat.mass, 2.0)

    def test_mul_num(self):
        mat = self.u235 * 2.0
        assert_equal(mat.mass, 2.0)

    def test_rmul_num(self):
        mat = 150 * self.u235
        assert_equal(mat.mass, 150.0)

    def test_div_num(self):
        mat = self.u235 / 10
        assert_equal(mat.mass, 0.1)


#
# Test atom fraction functions
#
def test_to_atom_frac():
    h2o = {10010000: 0.11191487328808077, 80160000: 0.8880851267119192}
    mat = Material(h2o, atoms_per_molecule=3.0)
    af = mat.to_atom_frac()
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_equal(af[10010000], 2.0)
    assert_equal(af[80160000], 1.0)
    assert_equal(mat.molecular_mass(), 18.01056468403)


def test_from_atom_frac_meth():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.mass, 18.01056468403)
    assert_equal(mat.molecular_mass(), 18.01056468403)

    h2 = Material({10010000: 1.0}, atoms_per_molecule=2.0)
    h2o = {'O16': 1.0, h2: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.molecular_mass(), 18.01056468403)

    mt1 = from_atom_frac({1001: 0.1, 6000: 0.8, 8016: 0.1})
    assert_equal(mt1.comp[10010000], 0.008911815984674479)
    assert_equal(mt1.comp[60000000], 0.849651197215362)
    assert_equal(mt1.comp[80160000], 0.14143698679996367)
    assert_equal(mt1.molecular_mass(), 11.3088626825682)

    ihm = Material()
    ihm.from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, 'O16': 2.0}
    mat = Material()
    mat.from_atom_frac(uox)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_almost_equal(mat.comp[80160000], 0.11912625367051276, 16)
    assert_almost_equal(mat.comp[922350000], 0.43763757904405304, 15)
    assert_almost_equal(mat.comp[922380000], 0.44323616728543414, 15)
    assert_almost_equal(mat.molecular_mass()/268.53718851614, 1.0, 15)


def test_to_atom_dens():
    h2o = {10010000: 0.11191487328808077, 80160000: 0.8880851267119192}
    mat = Material(h2o, density=1.0)
    ad = mat.to_atom_dens()
    assert_almost_equal(ad[10010000]/(10.**22), 6.68734335169385)
    assert_almost_equal(ad[80160000]/(10.**22), 3.34367167584692)
    
#
# Test mapping functions
#

def test_len():
    mat = Material(nucvec)
    assert_equal(len(mat), 9)


def test_contains():
    mat = Material(nucvec)
    assert_true(10010000 in mat)
    assert_true(922350000 in mat)
    assert_false(92000000 in mat)
    assert_raises(TypeError, lambda: 'word' in mat)


def test_getitem_int():
    mat = Material(nucvec)
    assert_equal(mat[922350000], 1.0)
    assert_raises(KeyError, lambda: mat[42])

    mat = Material(leu)
    assert_equal(mat[922350000], 0.04)
    assert_equal(mat[922380000], 0.96)
    assert_raises(KeyError, lambda: mat[922340000])


def test_getitem_str():
    mat = Material(nucvec)
    assert_equal(mat['U235'], 1.0)
    assert_raises(RuntimeError, lambda: mat['word'])

    mat = Material(leu)
    assert_equal(mat['U235'], 0.04)
    assert_equal(mat['U238'], 0.96)
    assert_raises(KeyError, lambda: mat['U234'])


def test_getitem_slice_int():
    mat = Material(nucvec)
    mat1 = mat[920000000:930000000]
    assert_equal(mat1.mass, 2.0)
    for nuc in mat1:
        assert_true(920000000 <= nuc < 930000000)

    mat1 = mat[:922380000]
    assert_equal(mat1.mass, 4.0)
    for nuc in mat1:
        assert_true(nuc < 922380000)

    mat1 = mat[922350000:]
    assert_equal(mat1.mass, 6.0)
    for nuc in mat1:
        assert_true(922350000 <= nuc)


def test_getitem_slice_str():
    mat = Material(nucvec)
    mat1 = mat['U':'NP']
    assert_equal(mat1.mass, 2.0)
    for nuc in mat1:
        assert_true(920000000 <= nuc < 930000000)

    mat1 = mat[:'U238']
    assert_equal(mat1.mass, 4.0)
    for nuc in mat1:
        assert_true(nuc < 922380000)

    mat1 = mat['U235':]
    assert_equal(mat1.mass, 6.0)
    for nuc in mat1:
        assert_true(922350000 <= nuc)


def test_getitem_sequence():
    mat = Material(nucvec)
    mat1 = mat[922380000, 922350000]
    assert_equal(mat1.mass, 2.0)
    assert_equal(set(mat1), set([922380000, 922350000]))

    mat1 = mat[922380, 'H2', 'h1']
    assert_equal(mat1.mass, 2.0)
    assert_equal(set(mat1), set([922380000, 10010000]))


def test_setitem_int():
    mat = Material(nucvec)
    assert_equal(mat.mass, 9.0)
    assert_equal(mat[922350000], 1.0)
    mat[922350000] = 2.0
    assert_equal(mat.mass, 10.0)
    assert_equal(mat[922350000], 2.0)

    mat = Material(leu)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat[922350000], 0.04)
    assert_equal(mat[922380000], 0.96)
    assert_raises(KeyError, lambda: mat[922340000])
    mat[922340000] = 17.0
    assert_equal(mat.mass, 18.0)
    assert_equal(mat[922340000], 17.0)
    assert_equal(mat[922350000], 0.04)
    assert_equal(mat[922380000], 0.96)



def test_setitem_str():
    mat = Material(nucvec)
    assert_equal(mat.mass, 9.0)
    assert_equal(mat[922350000], 1.0)
    mat['U235'] = 2.0
    assert_equal(mat.mass, 10.0)
    assert_equal(mat[922350000], 2.0)

    mat = Material(leu)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat[922350000], 0.04)
    assert_equal(mat[922380000], 0.96)
    assert_raises(KeyError, lambda: mat[922340000])
    mat['U234'] = 17.0
    assert_equal(mat.mass, 18.0)
    assert_equal(mat[922340000], 17.0)
    assert_equal(mat[922350000], 0.04)
    assert_equal(mat[922380000], 0.96)



def test_setitem_slice_int():
    mat = Material(nucvec)
    mat_id = id(mat)
    mat[920000000:930000000] = 42
    assert_equal(mat_id, id(mat))
    assert_equal(mat.mass, 91.0)
    assert_equal(mat[10010000], 1.0)
    assert_equal(mat[922350000], 42.0)
    assert_equal(mat[922380000], 42.0)

    mat = Material(nucvec)
    mat[:922380000] = 0.0
    assert_equal(mat.mass, 5.0)
    for nuc in mat:
        if (nuc < 922380000):
            assert_equal(mat[nuc], 0.0)
        else:
            assert_equal(mat[nuc], 1.0)

    mat = Material(nucvec)
    mat[922350000:] = 2
    assert_equal(mat.mass, 15.0)
    for nuc in mat:
        if (922350000 <= nuc):
            assert_equal(mat[nuc], 2.0)
        else:
            assert_equal(mat[nuc], 1.0)



def test_setitem_slice_str():
    mat = Material(nucvec)
    mat['U':'Np'] = 42
    assert_equal(mat.mass, 91.0)
    assert_equal(mat[10010000], 1.0)
    assert_equal(mat[922350000], 42.0)
    assert_equal(mat[922380000], 42.0)

    mat = Material(nucvec)
    mat[:'U238'] = 0.0
    assert_equal(mat.mass, 5.0)
    for nuc in mat:
        if (nuc < 922380000):
            assert_equal(mat[nuc], 0.0)
        else:
            assert_equal(mat[nuc], 1.0)

    mat = Material(nucvec)
    mat['U235':] = 2
    assert_equal(mat.mass, 15.0)
    for nuc in mat:
        if (922350000 <= nuc):
            assert_equal(mat[nuc], 2.0)
        else:
            assert_equal(mat[nuc], 1.0)


def test_setitem_sequence():
    mat = Material(nucvec)
    mat_id = id(mat)
    mat[922380000, 922350000] = 42
    assert_equal(mat_id, id(mat))
    assert_equal(mat.mass, 91.0)
    assert_equal(mat[10010000], 1.0)
    assert_equal(mat[922350000], 42.0)
    assert_equal(mat[922380000], 42.0)

    mat = Material(nucvec)
    mat[922380000, 'H2', 'h1'] = 0.0
    assert_equal(mat.mass, 7.0)
    assert_equal(len(mat), 10)
    for nuc in mat:
        if (nuc in [10010000, 10020000, 922380000]):
            assert_equal(mat[nuc], 0.0)
        else:
            assert_equal(mat[nuc], 1.0)

    mat = Material(nucvec)
    mat['U235', 'H2', 'h1'] = 2
    assert_equal(mat.mass, 13.0)
    assert_equal(len(mat), 10)
    for nuc in mat:
        if (nuc in [10010000, 10020000, 922350000]):
            assert_equal(mat[nuc], 2.0)
        else:
            assert_equal(mat[nuc], 1.0)



def test_delitem_int():
    mat = Material(nucvec)
    assert_equal(mat[922350000], 1.0)
    del mat[922350000]
    assert_raises(KeyError, lambda: mat[922350000])

    mat = Material(leu)
    assert_equal(mat[922350000], 0.04)
    del mat[922350000]
    assert_equal(mat.mass, 0.96)
    assert_equal(mat.comp[922380000], 1.0)
    assert_raises(KeyError, lambda: mat[922350000])


def test_delitem_str():
    mat = Material(nucvec)
    assert_equal(mat[922350000], 1.0)
    del mat['U235']
    assert_raises(KeyError, lambda: mat[922350000])

    mat = Material(leu)
    assert_equal(mat[922350000], 0.04)
    del mat['U235']
    assert_equal(mat.mass, 0.96)
    assert_equal(mat.comp[922380000], 1.0)
    assert_raises(KeyError, lambda: mat[922350000])


def test_delitem_slice_int():
    mat = Material(nucvec)
    del mat[920000000:930000000]
    assert_equal(mat.mass, 7.0)
    assert_equal(mat[10010000], 1.0)
    assert_raises(KeyError, lambda: mat[922350000])
    assert_raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[:922380000]
    assert_equal(mat.mass, 5.0)
    for nuc in mat:
        if (nuc < 922380000):
            assert_raises(KeyError, lambda: mat[nuc])
        else:
            assert_equal(mat[nuc], 1.0)

    mat = Material(nucvec)
    del mat[922350000:]
    assert_equal(mat.mass, 3.0)
    for nuc in mat:
        if (922350000 <= nuc):
            assert_raises(KeyError, lambda: mat[nuc])
        else:
            assert_equal(mat[nuc], 1.0)


def test_delitem_slice_str():
    mat = Material(nucvec)
    del mat['U':'Np']
    assert_equal(mat.mass, 7.0)
    assert_equal(mat[10010000], 1.0)
    assert_raises(KeyError, lambda: mat[922350000])
    assert_raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[:'U238']
    assert_equal(mat.mass, 5.0)
    for nuc in mat:
        if (nuc < 922380000):
            assert_raises(KeyError, lambda: mat[nuc])
        else:
            assert_equal(mat[nuc], 1.0)

    mat = Material(nucvec)
    del mat['U235':]
    assert_equal(mat.mass, 3.0)
    for nuc in mat:
        if (922350000 <= nuc):
            assert_raises(KeyError, lambda: mat[nuc])
        else:
            assert_equal(mat[nuc], 1.0)


def test_delitem_sequence():
    mat = Material(nucvec)
    del mat[922380000, 922350000]
    assert_equal(mat.mass, 7.0)
    assert_equal(mat[10010000], 1.0)
    assert_raises(KeyError, lambda: mat[922350000])
    assert_raises(KeyError, lambda: mat[922380000])

    mat = Material(nucvec)
    del mat[922380000, 'H2', 'h1']
    assert_equal(mat.mass, 7.0)
    assert_raises(KeyError, lambda: mat[10010000])
    assert_raises(KeyError, lambda: mat[10020000])
    assert_raises(KeyError, lambda: mat[922380000])



def test_iter():
    mat = Material(nucvec)
    for nuc in mat:
        assert_equal(mat[nuc], 1.0)

    mat = Material(leu)
    keys = set([922350000, 922380000])
    values = set([0.04, 0.96])
    items = set([(922350000, 0.04), (922380000, 0.96)])

    assert_equal(set(mat), keys)
    assert_equal(set(mat.values()), values)
    assert_equal(set(mat.items()), items)



#
# test material generation functions
#

def test_from_atom_frac_func():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = from_atom_frac(h2o)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.mass, 18.01056468403)
    assert_equal(mat.molecular_mass(), 18.01056468403)

    h2 = Material({10010000: 1.0}, atoms_per_molecule=2.0)
    h2o = {'O16': 1.0, h2: 1.0}
    mat = from_atom_frac(h2o)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.molecular_mass(), 18.01056468403)

    ihm = from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, 'O16': 2.0}
    mat = from_atom_frac(uox)
    assert_equal(mat.atoms_per_molecule, 3.0)
    assert_almost_equal(mat.comp[80160000], 0.11912625367051276, 16)
    assert_almost_equal(mat.comp[922350000], 0.43763757904405304, 15)
    assert_almost_equal(mat.comp[922380000], 0.44323616728543414, 15)
    assert_almost_equal(mat.molecular_mass()/268.53718851614, 1.0, 15)



def test_from_hdf5_func_protocol_0():
    mat = from_hdf5("mat.h5", "/mat", protocol=0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat = from_hdf5("mat.h5", "/mat", 0, 0)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat.comp, {922350000: 1.0, 942390000: 0.0})

    mat = from_hdf5("mat.h5", "/mat", 1, 0)
    assert_equal(mat.mass, 0.5)
    assert_equal(mat.comp, {922350000: 0.75, 942390000: 0.25})

    mat = from_hdf5("mat.h5", "/mat", 2, 0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat = from_hdf5("mat.h5", "/mat", -1, 0)
    assert_equal(mat.mass, 0.0)
    assert_equal(mat.comp, {922350000: 0.0, 942390000: 0.0})

    mat = from_hdf5("mat.h5", "/mat", -2, 0)
    assert_equal(mat.mass, 0.5)
    assert_equal(mat.comp, {922350000: 0.75, 942390000: 0.25})

    mat = from_hdf5("mat.h5", "/mat", -3, 0)
    assert_equal(mat.mass, 1.0)
    assert_equal(mat.comp, {922350000: 1.0, 942390000: 0.0})


def test_from_text_func():
    mat= from_text("mat.txt")
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})



def test_map_str_material():
    m = MapStrMaterial()
    m['leu'] = Material(leu)
    m['heu'] = Material({'U238': 0.01, 'U235': 0.99}, 42.0)
    assert_equal(len(m), 2)
    assert_equal(m['leu'].mass, 1.0)
    assert_equal(m['leu']['U235'], 0.04)
    assert_equal(m['heu'].mass, 42.0)
    assert_equal(m['heu']['U238'], 0.42)

    m = MapStrMaterial({'leu': Material(leu), 'heu': Material({'U238': 0.01, 'U235': 0.99}, 42.0)})
    assert_equal(len(m), 2)
    assert_equal(m['leu'].mass, 1.0)
    assert_equal(m['leu']['U235'], 0.04)
    assert_equal(m['heu'].mass, 42.0)
    assert_equal(m['heu']['U238'], 0.42)

    n = MapStrMaterial(m, False)
    assert_equal(len(n), 2)
    assert_equal(n['leu'].mass, 1.0)
    assert_equal(n['leu']['U235'], 0.04)
    assert_equal(n['heu'].mass, 42.0)
    assert_equal(n['heu']['U238'], 0.42)

    # points to the same underlying map
    n['other'] = Material({'PU239': 15.0})
    assert_equal(m['other'].mass, 15.0)
    assert_equal(m['other']['PU239'], 15.0)

    assert_equal(n['leu'].mass, 1.0)
    assert_equal(n['leu']['U235'], 0.04)
    assert_equal(n['heu'].mass, 42.0)
    assert_equal(n['heu']['U238'], 0.42)


def test_metadata():
    mat = Material(leu)
    assert_equal(len(mat.metadata), 0)
    mat.metadata['units'] = 'kg'
    assert_equal(len(mat.metadata), 1)
    assert_equal(mat.metadata['units'], 'kg')

    mat.metadata = {'comment': 'rawr', 'amount': 42.0}
    assert_equal(mat.metadata.keys(), ['amount', 'comment'])
    assert_true(isinstance(mat.metadata, jsoncpp.Value))

    aview = mat.metadata
    aview['omnomnom'] = [1, 2, 5, 3]
    assert_equal(len(mat.metadata), 3)
    assert_equal(list(mat.metadata['omnomnom']), [1, 2, 5, 3])
#
# Test MultiMaterial
#
def test_multimaterial_mix_composition():
    mat1 = Material(nucvec={120240000:0.3, 300000000:0.2, 10010000:0.1}, density=2.71)
    mat2 = Material(nucvec={60120000:0.2, 280640000:0.5, 10010000:0.12}, density=8.0)
    mix = MultiMaterial({mat1:0.5, mat2:0.21})
    mat3 = mix.mix_by_mass()
    mat4 = mix.mix_by_volume()

    assert_equal(mat3.comp[10010000], 0.16065498683155846)
    assert_equal(mat3.comp[60120000], 0.0721401580212985)
    assert_equal(mat3.comp[120240000], 0.352112676056338)
    assert_equal(mat3.comp[280640000], 0.18035039505324627)
    assert_equal(mat3.comp[300000000], 0.2347417840375587)

    assert_equal(mat4.comp[10010000], 0.15541581280722197)
    assert_equal(mat4.comp[60120000], 0.13501024631333625)
    assert_equal(mat4.comp[120240000], 0.2232289950576606)
    assert_equal(mat4.comp[280640000], 0.33752561578334067)
    assert_equal(mat4.comp[300000000], 0.14881933003844042)


def test_multimaterial_mix_density():
    mat1 = Material(nucvec={120240000:0.3, 300000000:0.2, 10010000:0.1}, density=1.0)
    mat2 = Material(nucvec={60120000:0.2, 280640000:0.5, 10010000:0.12}, density=2.0)
    # mixing by hand to get density 1.5 when 50:50 by volume of density 1 & 2
    mix = MultiMaterial({mat1:0.5, mat2:0.5})
    mat3 = mix.mix_by_volume()
    # calculated mass fracs by hand, same problem as above stated in mass fraction terms
    # rather than volume fraction terms.
    mix = MultiMaterial({mat1:1/3., mat2:2/3.})
    mat4 = mix.mix_by_mass()

    assert_equal(mat3.density, 1.5)
    assert_equal(mat4.density, 1.5)

    assert_equal(mat3.density, mat4.density)

def test_deepcopy():
    x = Material({'H1': 1.0}, mass=2.0, density=3.0, atoms_per_molecule=4.0,
                 metadata={'name': 'loki'})
    y = deepcopy(x)
    assert_equal(x, y)
    y.comp[10010000] = 42.0
    y[80160000] = 21.0
    y.density = 65.0
    y.atoms_per_molecule = 28.0
    y.metadata['wakka'] = 'jawaka'
    y.mass = 48.0
    y.metadata['name'] = 'odin'
    assert_not_equal(x, y)
    assert_equal(x, Material({'H1': 1.0}, mass=2.0, density=3.0, atoms_per_molecule=4.0,
                             metadata={'name': 'loki'}))

def test_mcnp():

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96},
                   metadata={'mat_number': 2,
                          'table_ids': {'92235':'15c', '92238':'25c'},
                          'mat_name':'LEU',
                          'source':'Some URL',
                          'comments': ('this is a long comment that will definitly '
                                       'go over the 80 character limit, for science'),
                          'name':'leu'},
                   density=19.1)

    mass = leu.mcnp()
    mass_exp = ('C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c -4.0000e-02\n'
                '     92238.25c -9.6000e-01\n')
    assert_equal(mass, mass_exp)

    atom = leu.mcnp(frac_type='atom')
    atom_exp = ('C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c 4.0491e-02\n'
                '     92238.25c 9.5951e-01\n')
    assert_equal(atom, atom_exp)

def test_mcnp_mat0():

    leu = Material(nucvec={'U235': 0.04, 'U236': 0.0, 'U238': 0.96},
                   metadata={'mat_number': 2,
                          'table_ids': {'92235':'15c', '92236':'15c', '92238':'25c'},
                          'mat_name':'LEU',
                          'source':'Some URL',
                          'comments': ('this is a long comment that will definitly '
                                       'go over the 80 character limit, for science'),
                          'name':'leu'},
                   density=19.1)

    mass = leu.mcnp()
    mass_exp = ('C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c -4.0000e-02\n'
                '     92238.25c -9.6000e-01\n')
    assert_equal(mass, mass_exp)


def test_alara():

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96}, metadata={\
          'mat_number':2, 'table_ids':{'922350':'15c', '922380':'25c'},\
          'name':'LEU', 'source':'Some URL', \
          'comments': \
'this is a long comment that will definitly go over the 80 character limit, for science', \
            }, density=19.1)
    leu2 = Material(nucvec={'U235': 0.04, 'U238': 0.96}, metadata={\
          'mat_number':2,}, density=19.1)
    leu3 = Material(nucvec={'U235': 0.04, 'U238': 0.96})


    written = leu.alara()
    written += leu2.alara()
    written += leu3.alara()

    expected = ('# mat number: 2\n'
                '# source: Some URL\n'
                '# comments: this is a long comment that will definitly go over the 80 character\n'
                '#  limit, for science\n'
                'LEU 19.1 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n'
                '# mat number: 2\n'
                'mat2_rho-19.1 19.1 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n'
                'mat<mat_num>_rho-<rho> <rho> 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n')
    assert_equal(written, expected)

def test_write_mcnp():
    if 'mcnp_mass_fracs.txt' in os.listdir('.'):
        os.remove('mcnp_mass_fracs.txt')

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96},
                   metadata={'mat_number': 2,
                          'table_ids': {'92235':'15c', '92238':'25c'},
                          'mat_name':'LEU',
                          'source':'Some URL',
                          'comments': ('this is a long comment that will definitly '
                                       'go over the 80 character limit, for science'),
                          'name':'leu'},
                   density=19.1)

    leu.write_mcnp('mcnp_mass_fracs.txt')
    leu.write_mcnp('mcnp_mass_fracs.txt', frac_type='atom')

    with open('mcnp_mass_fracs.txt') as f:
        written = f.read()
    expected = ('C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c -4.0000e-02\n'
                '     92238.25c -9.6000e-01\n'
                'C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c 4.0491e-02\n'
                '     92238.25c 9.5951e-01\n')
    assert_equal(written, expected)
    os.remove('mcnp_mass_fracs.txt')

def test_fluka():
    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96},
                   metadata={'mat_number': 2,
                          'table_ids': {'92235':'15c', '92238':'25c'},
                          'name':'LEU',
                          'fluka_name':'URANIUM',
                          'fluka_material_index': '35',
                          'source':'Some URL',
                          'comments': ('Fluka Compound '),
                          },
                   density=19.1)
    ########################################
    # Part I:  Do not collapse the materials
    id = 25
    matlines = []
    # call fluka() on a material made up of each component
    for key in leu.comp:
        element = Material(nucvec={key:1})
        matlines.append(element.fluka(id,'atom'))
        id=id+1
    compound = leu.fluka(id,'atom')
    matlines.append(compound)
    written = ''.join(matlines)

    exp =  'MATERIAL         92.   235.044        1.       25.                    235-U     \n'
    exp += 'MATERIAL         92.   238.051        1.       26.                    238-U     \n'
    exp += '* Fluka Compound \n'
    exp += 'MATERIAL          1.        1.      19.1       27.                    URANIUM   \n'
    exp += 'COMPOUND   4.000e-02     235-U 9.600e-01     238-U                    URANIUM   \n'

    assert_equal(exp,written)

    #####################################
    # Part II:  Test a collapsed material
    coll = leu.collapse_elements({920000000})
    coll.metadata['comments'] = 'Fluka Element '

    exp = '* Fluka Element'
    exp += ' \n'
    exp += 'MATERIAL         92.   238.029      19.1       25.                    URANIUM'
    exp += '   \n'

    written = coll.fluka(25,'atom')
    assert_equal(exp, written)

    ##################################
    # Repeat Part I for mass frac_type
    id = 25
    matlines = []
    # call fluka() on a material made up of each component
    for key in leu.comp:
        element = Material(nucvec={key:1})
        matlines.append(element.fluka(id,'mass'))
        id=id+1
    compound = leu.fluka(id,'mass')
    matlines.append(compound)
    written = ''.join(matlines)

    exp =  'MATERIAL         92.   235.044        1.       25.                    235-U     \n'
    exp += 'MATERIAL         92.   238.051        1.       26.                    238-U     \n'
    exp += '* Fluka Compound \n'
    exp += 'MATERIAL          1.        1.      19.1       27.                    URANIUM   \n'
    exp += 'COMPOUND  -4.000e-02     235-U-9.600e-01     238-U                    URANIUM   \n'
    assert_equal(exp,written)

def test_fluka_scientific():
    # baseline test for scientific formatting
    mat = Material({'H':0.1,'O':0.8,'C':0.1})
    mat.density=1.0
    mat.metadata['fluka_name'] = 'ORGPOLYM'
    written = mat.fluka(25)

    exp  = 'MATERIAL          1.        1.        1.       25.                    ORGPOLYM  \n'
    exp += 'COMPOUND  -1.000e-01  HYDROGEN-1.000e-01    CARBON-8.000e-01    OXYGENORGPOLYM  \n'
    assert_equal(exp,written)

    mat = Material({'H':0.01,'O':0.8,'C':0.19})
    mat.density=1.0
    mat.metadata['fluka_name'] = 'ORGPOLYM'
    written = mat.fluka(25)
    
    exp =  'MATERIAL          1.        1.        1.       25.                    ORGPOLYM  \n'
    exp += 'COMPOUND  -1.000e-02  HYDROGEN-1.900e-01    CARBON-8.000e-01    OXYGENORGPOLYM  \n'
    assert_equal(exp,written)


    

def test_write_alara():
    if 'alara.txt' in os.listdir('.'):
        os.remove('alara.txt')

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96}, metadata={\
          'mat_number':2, 'table_ids':{'922350':'15c', '922380':'25c'},\
          'name':'LEU', 'source':'Some URL', \
          'comments': \
'this is a long comment that will definitly go over the 80 character limit, for science', \
            }, density=19.1)
    leu2 = Material(nucvec={'U235': 0.04, 'U238': 0.96}, metadata={\
          'mat_number':2,}, density=19.1)
    leu3 = Material(nucvec={'U235': 0.04, 'U238': 0.96})

    leu.write_alara('alara.txt')
    leu2.write_alara('alara.txt')
    leu3.write_alara('alara.txt')

    with open('alara.txt') as f:
        written = f.read()
    expected = ('# mat number: 2\n'
                '# source: Some URL\n'
                '# comments: this is a long comment that will definitly go over the 80 character\n'
                '#  limit, for science\n'
                'LEU 19.1 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n'
                '# mat number: 2\n'
                'mat2_rho-19.1 19.1 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n'
                'mat<mat_num>_rho-<rho> <rho> 2\n'
                '     u:235 4.0000E-02 92\n'
                '     u:238 9.6000E-01 92\n')
    assert_equal(written, expected)
    os.remove('alara.txt')

def test_natural_elements():
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    expected_comp = {10000000: 0.11189838783149784, 80000000: 0.8881016121685023}
    for key in expected_comp:
        assert_almost_equal(water.comp[key], expected_comp[key])


def test_load_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = Material(leu)
    obs = Material()
    json = jsoncpp.Value({"mass": 1.0, "comp": leu, "density": -1.0, "metadata": {},
                         "atoms_per_molecule": -1.0})
    obs.load_json(json)
    assert_equal(exp, obs)

def test_dump_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = jsoncpp.Value({"mass": 1.0, "comp": leu, "density": -1.0, "metadata": {},
                         "atoms_per_molecule": -1.0})
    obs = Material(leu).dump_json()
    assert_equal(exp, obs)

def test_rw_json():
    filename = "leu.json"
    wmat = Material(leu)
    wmat.write_json(filename)
    rmat = Material()
    rmat.from_json(filename)
    assert_equal(wmat, rmat)
    os.remove(filename)



#
#  Material Library
#


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
    assert_equal(set(wmatlib), set(rmatlib))
    for key in rmatlib:
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)

def test_matlib_hdf5_nuc_data():
    matlib = MaterialLibrary()
    matlib.from_hdf5(nuc_data, datapath="/material_library/materials",
                     nucpath="/material_library/nucid")

def test_matlib_hdf5():
    filename = "matlib.h5"
    if filename in os.listdir('.'):
        os.remove(filename)
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.metadata["name"] = "Aqua sera."
    lib = {"leu": Material(leu), "nucvec": nucvec, "aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_hdf5(filename)
    rmatlib = MaterialLibrary()
    rmatlib.from_hdf5(filename)
    os.remove(filename)
    # Round trip!
    rmatlib.write_hdf5(filename)
    wmatlib = MaterialLibrary(filename)
    assert_equal(set(wmatlib), set(rmatlib))
    for key in rmatlib:
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)


def test_material_gammas():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert_array_equal(mat.gammas(), [(19.55, np.nan),
                 (31.6, 2.1476136491211186e-20),
                 (34.7, 4.674217942204787e-20),
                 (41.4, 3.789906439625503e-20),
                 (41.96, 7.579812879251006e-20),
                 (51.21, 4.295227298242237e-20),
                 (54.1, 1.2633021465418345e-21),
                 (54.25, 3.789906439625503e-20),
                 (60.5, np.nan),
                 (64.45, np.nan),
                 (72.7, 1.5159625758502012e-19),
                 (74.94, 6.442840947363355e-20),
                 (76.2, np.nan),
                 (94.0, np.nan),
                 (95.7, np.nan),
                 (96.09, 1.1496049533530693e-19),
                 (109.19, 2.097081563259445e-18),
                 (115.45, 3.789906439625503e-20),
                 (120.35, 3.2845855810087694e-20),
                 (136.55, 1.5159625758502013e-20),
                 (140.76, 2.526604293083669e-19),
                 (142.4, 6.316510732709172e-21),
                 (143.76, 1.3845791526098507e-17),
                 (147.0, np.nan),
                 (150.93, 1.136971931887651e-19),
                 (163.356, 6.4175749044325184e-18),
                 (173.3, 7.579812879251007e-21),
                 (182.1, np.nan),
                 (182.62, 4.926878371513154e-19),
                 (185.715, 7.200822235288456e-17),
                 (194.94, 7.958803523213557e-19),
                 (198.9, 4.5478877275506033e-20),
                 (202.12, 1.3643663182651812e-18),
                 (205.316, 6.3417767756400084e-18),
                 (215.28, 3.66357622497132e-20),
                 (221.386, 1.4906965329193644e-19),
                 (228.78, 8.843115025792842e-21),
                 (233.5, 4.800548156858971e-20),
                 (240.88, 9.348435884409573e-20),
                 (246.83, 6.948161805980088e-20),
                 (266.45, 7.579812879251007e-21),
                 (275.35, 6.442840947363355e-20),
                 (275.49, 4.0425668689338704e-20),
                 (281.42, 7.579812879251007e-21),
                 (282.92, 7.579812879251007e-21),
                 (289.56, 8.843115025792842e-21),
                 (291.2, np.nan),
                 (291.65, 5.0532085861673375e-20),
                 (301.7, 6.316510732709172e-21),
                 (317.1, 1.2633021465418345e-21),
                 (343.5, 3.789906439625503e-21),
                 (345.92, 5.0532085861673375e-20),
                 (356.03, 6.316510732709172e-21),
                 (387.84, 5.0532085861673375e-20),
                 (410.29, 3.789906439625503e-21),
                 (428.7, np.nan),
                 (448.4, 1.2633021465418345e-21),
                 (49.55, 3.0188209868525946e-19),
                 (113.5, 4.811245947796322e-20)])

def test_material_xrays():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert_equal(mat.xrays(), [(93.35, 7.135646053967512e-18),
     (89.953, 4.4112563905627166e-18),
     (105.0, 3.394789318691889e-18),
     (13.0, 2.965225083810098e-17),
     (93.35, 5.2767097758422435e-21),
     (89.953, 3.262061983425675e-21),
     (105.0, 2.5103988972247685e-21),
     (13.0, 3.4433858449448927e-17)])

def test_material_photons():
    leu = {"U238": 0.96, "U235": 0.04}
    mat = Material(leu)
    assert_array_equal(mat.photons(), [(19.55, np.nan),
                 (31.6, 2.1476136491211186e-20),
                 (34.7, 4.674217942204787e-20),
                 (41.4, 3.789906439625503e-20),
                 (41.96, 7.579812879251006e-20),
                 (51.21, 4.295227298242237e-20),
                 (54.1, 1.2633021465418345e-21),
                 (54.25, 3.789906439625503e-20),
                 (60.5, np.nan),
                 (64.45, np.nan),
                 (72.7, 1.5159625758502012e-19),
                 (74.94, 6.442840947363355e-20),
                 (76.2, np.nan),
                 (94.0, np.nan),
                 (95.7, np.nan),
                 (96.09, 1.1496049533530693e-19),
                 (109.19, 2.097081563259445e-18),
                 (115.45, 3.789906439625503e-20),
                 (120.35, 3.2845855810087694e-20),
                 (136.55, 1.5159625758502013e-20),
                 (140.76, 2.526604293083669e-19),
                 (142.4, 6.316510732709172e-21),
                 (143.76, 1.3845791526098507e-17),
                 (147.0, np.nan),
                 (150.93, 1.136971931887651e-19),
                 (163.356, 6.4175749044325184e-18),
                 (173.3, 7.579812879251007e-21),
                 (182.1, np.nan),
                 (182.62, 4.926878371513154e-19),
                 (185.715, 7.200822235288456e-17),
                 (194.94, 7.958803523213557e-19),
                 (198.9, 4.5478877275506033e-20),
                 (202.12, 1.3643663182651812e-18),
                 (205.316, 6.3417767756400084e-18),
                 (215.28, 3.66357622497132e-20),
                 (221.386, 1.4906965329193644e-19),
                 (228.78, 8.843115025792842e-21),
                 (233.5, 4.800548156858971e-20),
                 (240.88, 9.348435884409573e-20),
                 (246.83, 6.948161805980088e-20),
                 (266.45, 7.579812879251007e-21),
                 (275.35, 6.442840947363355e-20),
                 (275.49, 4.0425668689338704e-20),
                 (281.42, 7.579812879251007e-21),
                 (282.92, 7.579812879251007e-21),
                 (289.56, 8.843115025792842e-21),
                 (291.2, np.nan),
                 (291.65, 5.0532085861673375e-20),
                 (301.7, 6.316510732709172e-21),
                 (317.1, 1.2633021465418345e-21),
                 (343.5, 3.789906439625503e-21),
                 (345.92, 5.0532085861673375e-20),
                 (356.03, 6.316510732709172e-21),
                 (387.84, 5.0532085861673375e-20),
                 (410.29, 3.789906439625503e-21),
                 (428.7, np.nan),
                 (448.4, 1.2633021465418345e-21),
                 (49.55, 3.0188209868525946e-19),
                 (113.5, 4.811245947796322e-20),
                 (93.35, 7.135646053967512e-18),
                 (89.953, 4.4112563905627166e-18),
                 (105.0, 3.394789318691889e-18),
                 (13.0, 2.965225083810098e-17),
                 (93.35, 5.2767097758422435e-21),
                 (89.953, 3.262061983425675e-21),
                 (105.0, 2.5103988972247685e-21),
                 (13.0, 3.4433858449448927e-17)])
    assert_equal(mat.photons(True), [(31.6, 0.00011635441715962317),
                 (34.7, 0.0002532419667591798),
                 (41.4, 0.000205331324399335),
                 (41.96, 0.00041066264879867),
                 (51.21, 0.00023270883431924635),
                 (54.1, 6.844377479977834e-06),
                 (54.25, 0.000205331324399335),
                 (72.7, 0.00082132529759734),
                 (74.94, 0.00034906325147886946),
                 (96.09, 0.0006228383506779829),
                 (109.19, 0.011361666616763204),
                 (115.45, 0.000205331324399335),
                 (120.35, 0.00017795381447942367),
                 (136.55, 8.213252975973401e-05),
                 (140.76, 0.0013688754959955667),
                 (142.4, 3.422188739988917e-05),
                 (143.76, 0.07501437718055706),
                 (150.93, 0.0006159939731980049),
                 (163.356, 0.03476943759828739),
                 (173.3, 4.1066264879867005e-05),
                 (182.62, 0.002669307217191355),
                 (185.715, 0.39012951635873655),
                 (194.94, 0.004311957812386035),
                 (198.9, 0.00024639758927920197),
                 (202.12, 0.007391927678376061),
                 (205.316, 0.034358774949488725),
                 (215.28, 0.00019848694691935719),
                 (221.386, 0.0008076365426373843),
                 (228.78, 4.791064235984484e-05),
                 (233.5, 0.00026008634423915766),
                 (240.88, 0.0005064839335183596),
                 (246.83, 0.0003764407613987808),
                 (266.45, 4.1066264879867005e-05),
                 (275.35, 0.00034906325147886946),
                 (275.49, 0.00021902007935929068),
                 (281.42, 4.1066264879867005e-05),
                 (282.92, 4.1066264879867005e-05),
                 (289.56, 4.791064235984484e-05),
                 (291.65, 0.00027377509919911336),
                 (301.7, 3.422188739988917e-05),
                 (317.1, 6.844377479977834e-06),
                 (343.5, 2.0533132439933503e-05),
                 (345.92, 0.00027377509919911336),
                 (356.03, 3.422188739988917e-05),
                 (387.84, 0.00027377509919911336),
                 (410.29, 2.0533132439933503e-05),
                 (448.4, 6.844377479977834e-06),
                 (49.55, 0.0016355509594484913),
                 (113.5, 0.0002606659341621033),
                 (93.35, 0.038659837071091864),
                 (89.953, 0.023899511277348993),
                 (105.0, 0.018392448414441622),
                 (13.0, 0.1606513520320623),
                 (93.35, 2.858840512304754e-05),
                 (89.953, 1.767335204706799e-05),
                 (105.0, 1.3600956608013968e-05),
                 (13.0, 0.18655736948227228)])


def test_decay_h3():
    mat = Material({'H3': 1.0})
    obs = mat.decay(data.half_life('H3'))
    obs = obs.to_atom_frac()
    assert_equal(2, len(obs))
    assert_almost_equal(0.5, obs[nucname.id('H3')])
    assert_almost_equal(0.5, obs[nucname.id('He3')])

def test_decay_u235_h3():
    mat = Material({'U235': 1.0, 'H3': 1.0})
    obs = mat.decay(365.25 * 24.0 * 3600.0)
    if len(obs) < 4:
        # full decay is not installed
        raise SkipTest
    exp = Material({10030000: 0.472645829913563, 20030000: 0.02735407958518153, 
            812070000: 8.609275765376378e-22, 822090000: 5.262160149924099e-29,
            822110000: 1.2876726344327517e-20, 832110000: 8.369159135473628e-22,
            832150000: 4.595084328237903e-27, 842110000: 6.3568211371206766e-27, 
            842150000: 1.0842337907192214e-26, 852190000: 4.023654648573945e-28,
            862190000: 3.796661905033895e-24, 872230000: 1.5415521905397778e-22,
            882230000: 4.431826308123027e-18, 882270000: 2.2016578580798292e-26,
            892270000: 4.7669120084704354e-15, 902270000: 1.1717834799662062e-17,
            902310000: 2.0323933632662482e-12, 912310000: 4.83892288435637e-10,
            922350000: 0.5000000900153261}, 1.99999963797, -1.0, 1.0, {})
    assert_mat_almost_equal(exp, obs)

# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
