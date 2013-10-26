"""Material tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

import os
from pyne import nuc_data
from pyne.material import Material, from_atom_frac, from_hdf5, from_text, \
    MapStrMaterial, MultiMaterial, MaterialLibrary
from pyne import jsoncpp 
from pyne import data
import numpy  as np
import tables as tb

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
    assert_almost_equal(first.atoms_per_mol, second.atoms_per_mol, places=places)
    assert_equal(first.attrs, second.attrs)
    nucs = set(second.comp.keys())
    assert_equal(set(first.comp.keys()), nucs)
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
    mat = Material({922350000: 0.05, 922380000: 0.95}, 15, attrs={'units': 'kg'})
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.mass, 15.0)
    assert_equal(mat.attrs['units'], 'kg')

def test_from_text():
    mat = Material(attrs={'units': 'kg'})
    mat.from_text("mat.txt")
    assert_equal(mat.comp, {922350000: 0.05, 922380000: 0.95})
    assert_equal(mat.attrs['units'], 'kg')


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
    assert_equal(leu.atoms_per_mol, read_leu.atoms_per_mol)
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
    leu.attrs['comment'] = 'first light'
    leu.write_hdf5('proto1.h5', chunksize=10)

    for i in range(2, 11):
        leu = Material({'U235': 0.04, 'U238': 0.96}, i*4.2, 2.72, 1.0*i)
        leu.attrs['comment'] = 'fire in the disco - {0}'.format(i)
        leu.write_hdf5('proto1.h5')

    # Loads with protocol 1 now.
    m = Material()
    m.from_hdf5('proto1.h5', '/material', -3, 1)
    assert_equal(m.density, 2.72)
    assert_equal(m.atoms_per_mol, 8.0)
    assert_equal(m.mass, 33.6)
    assert_equal(m.comp, {922350000: 0.04, 922380000: 0.96})
    assert_equal(m.attrs['comment'], 'fire in the disco - 8')

    m = from_hdf5('proto1.h5', '/material', 3, 1)
    assert_equal(m.density, 2.72)
    assert_equal(m.atoms_per_mol, 4.0)
    assert_equal(m.mass, 16.8)
    assert_equal(m.comp, {922350000: 0.04, 922380000: 0.96})
    assert_equal(m.attrs['comment'], 'fire in the disco - 4')

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


    def test_molecular_weight(self):
        mat_empty = Material({})
        assert_equal(mat_empty.molecular_weight(), 0.0)

        mat_u238 = Material({922380000: 1.0})
        mw_u238 = mat_u238.molecular_weight()
        try:
            assert_almost_equal(mw_u238, 238.050788423)
        except AssertionError:
            assert_almost_equal(mw_u238, 238.0)            

        mat_mixed = Material({922350000: 0.5, 922380000: 0.5})
        mw_mixed = mat_mixed.molecular_weight()
        try:
            assert_almost_equal(mw_mixed/236.547360417, 1.0, 4)            
        except AssertionError:
            assert_almost_equal(mw_mixed/236.5, 1.0, 4)


def test_expand_elements1():
    natmat = Material({'C': 1.0, 902320000: 0.5, 'PU': 4.0, 'U': 3.0}, 
                       attrs={'y': 1.0})
    expmat = natmat.expand_elements()
    assert_true(60120000 in expmat.comp)
    assert_false(60000000 in expmat.comp)
    assert_true(natmat.attrs == expmat.attrs)
    assert_false(natmat.attrs is expmat.attrs)

def test_expand_elements2():
    """Inspired by #86"""
    natmat = Material({'C': 1.0})
    expmat = natmat.expand_elements()
    afrac = expmat.to_atom_frac()
    assert_almost_equal(data.natural_abund(60120000), afrac[60120000])
    assert_almost_equal(data.natural_abund(60130000), afrac[60130000])


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
    mat = Material(h2o, atoms_per_mol=3.0)
    af = mat.to_atom_frac()
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_equal(af[10010000], 2.0)
    assert_equal(af[80160000], 1.0)
    assert_equal(mat.molecular_weight(), 18.01056468403)    


def test_from_atom_frac_meth():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.mass, 18.01056468403)    
    assert_equal(mat.molecular_weight(), 18.01056468403)    

    h2 = Material({10010000: 1.0}, atoms_per_mol=2.0)
    h2o = {'O16': 1.0, h2: 1.0}
    mat = Material()
    mat.from_atom_frac(h2o)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.molecular_weight(), 18.01056468403)    

    ihm = Material()
    ihm.from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, 'O16': 2.0}
    mat = Material()
    mat.from_atom_frac(uox)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_almost_equal(mat.comp[80160000], 0.11912625367051276, 16)
    assert_almost_equal(mat.comp[922350000], 0.43763757904405304, 15)
    assert_almost_equal(mat.comp[922380000], 0.44323616728543414, 15)
    assert_almost_equal(mat.molecular_weight()/268.53718851614, 1.0, 15)


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
    assert_equal(set(mat1.keys()), set([922380000, 922350000]))

    mat1 = mat[922380, 'H2', 'h1']
    assert_equal(mat1.mass, 2.0)
    assert_equal(set(mat1.keys()), set([922380000, 10010000]))


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

    assert_equal(set(mat.keys()), keys)
    assert_equal(set(mat.values()), values)
    assert_equal(set(mat.items()), items)



# 
# test material generation functions
#

def test_from_atom_frac_func():
    h2o = {10010000: 2.0, 80160000: 1.0}
    mat = from_atom_frac(h2o)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.mass, 18.01056468403)    
    assert_equal(mat.molecular_weight(), 18.01056468403)    

    h2 = Material({10010000: 1.0}, atoms_per_mol=2.0)
    h2o = {'O16': 1.0, h2: 1.0}
    mat = from_atom_frac(h2o)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_equal(mat.comp[10010000], 0.11191487328808077)
    assert_equal(mat.comp[80160000], 0.8880851267119192)
    assert_equal(mat.molecular_weight(), 18.01056468403)    

    ihm = from_atom_frac({922350000: 0.5, 922380000: 0.5})
    uox = {ihm: 1.0, 'O16': 2.0}
    mat = from_atom_frac(uox)
    assert_equal(mat.atoms_per_mol, 3.0)
    assert_almost_equal(mat.comp[80160000], 0.11912625367051276, 16)
    assert_almost_equal(mat.comp[922350000], 0.43763757904405304, 15)
    assert_almost_equal(mat.comp[922380000], 0.44323616728543414, 15)
    assert_almost_equal(mat.molecular_weight()/268.53718851614, 1.0, 15)



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


def test_attrs():
    mat = Material(leu)
    assert_equal(len(mat.attrs), 0)
    mat.attrs['units'] = 'kg'
    assert_equal(len(mat.attrs), 1)
    assert_equal(mat.attrs['units'], 'kg')
    
    mat.attrs = {'comment': 'rawr', 'amount': 42.0}
    assert_equal(mat.attrs.keys(), ['amount', 'comment'])
    assert_true(isinstance(mat.attrs, jsoncpp.Value))

    aview = mat.attrs
    aview['omnomnom'] = [1, 2, 5, 3]
    assert_equal(len(mat.attrs), 3)
    assert_equal(list(mat.attrs['omnomnom']), [1, 2, 5, 3])
#
# Test MultiMaterial
#
def test_multimaterial():
    mat1 = Material(nucvec={120240000:0.3, 300000000:0.2, 10010000:0.1}, density=2.71)
    mat2 = Material(nucvec={60120000:0.2, 280640000:0.5, 10010000:0.12}, density=8.0)
    mix = MultiMaterial({mat1:0.5, mat2:0.21})
    mat3 = mix.mix_by_mass()
    mat4 = mix.mix_by_volume()

    assert_equal(mat3.density, -1.0)
    assert_equal(mat3.comp[10010000], 0.16065498683155846)
    assert_equal(mat3.comp[60120000], 0.0721401580212985)
    assert_equal(mat3.comp[120240000], 0.352112676056338)
    assert_equal(mat3.comp[280640000], 0.18035039505324627)
    assert_equal(mat3.comp[300000000], 0.2347417840375587)

    assert_equal(mat4.density, -1.0)
    assert_equal(mat4.comp[10010000], 0.15541581280722197)
    assert_equal(mat4.comp[60120000], 0.13501024631333625)
    assert_equal(mat4.comp[120240000], 0.2232289950576606)
    assert_equal(mat4.comp[280640000], 0.33752561578334067)
    assert_equal(mat4.comp[300000000], 0.14881933003844042)

def test_write_mcnp():
    if 'mcnp_mass_fracs.txt' in os.listdir('.'):
        os.remove('mcnp_mass_fracs.txt')

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96}, 
                   attrs={'mat_number': 2, 
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
                '     92235.15c -4.0000E-02\n'
                '     92238.25c -9.6000E-01\n'
                'C name: leu\n'
                'C density = 19.1\n'
                'C source: Some URL\n'
                'C comments: this is a long comment that will definitly go over the 80 character\n'
                'C  limit, for science\n'
                'm2\n'
                '     92235.15c 4.0491E-02\n'
                '     92238.25c 9.5951E-01\n')
    assert_equal(written, expected)
    os.remove('mcnp_mass_fracs.txt')


def test_write_alara():
    if 'alara.txt' in os.listdir('.'):
        os.remove('alara.txt')

    leu = Material(nucvec={'U235': 0.04, 'U238': 0.96}, attrs={\
          'mat_number':2, 'table_ids':{'922350':'15c', '922380':'25c'},\
          'name':'LEU', 'source':'Some URL', \
          'comments': \
'this is a long comment that will definitly go over the 80 character limit, for science', \
            }, density=19.1)
    leu2 = Material(nucvec={'U235': 0.04, 'U238': 0.96}, attrs={\
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
    for key in expected_comp.keys():
        assert_almost_equal(water.comp[key], expected_comp[key])


def test_load_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = Material(leu)
    obs = Material()
    json = jsoncpp.Value({"mass": 1.0, "comp": leu, "density": -1.0, "attrs": {}, 
                         "atoms_per_mol": -1.0})
    obs.load_json(json)
    assert_equal(exp, obs)

def test_dump_json():
    leu = {"U238": 0.96, "U235": 0.04}
    exp = jsoncpp.Value({"mass": 1.0, "comp": leu, "density": -1.0, "attrs": {}, 
                         "atoms_per_mol": -1.0})
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
    water.attrs["name"] = "Aqua sera."
    lib = {"leu": Material(leu), "nucvec": nucvec, "aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_json(filename)
    rmatlib = MaterialLibrary()
    rmatlib.from_json(filename)
    assert_equal(set(wmatlib.keys()), set(rmatlib.keys()))
    for key in rmatlib.keys():
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)

def test_matlib_hdf5_nuc_data():
    matlib = MaterialLibrary()
    matlib.from_hdf5(nuc_data, datapath="/material_library/materials", 
                     nucpath="/material_library/nucid")

def test_matlib_hdf5():
    filename = "matlib.h5"
    water = Material()
    water.from_atom_frac({10000000: 2.0, 80000000: 1.0})
    water.attrs["name"] = "Aqua sera."
    lib = {"leu": Material(leu), "nucvec": nucvec, "aqua": water}
    wmatlib = MaterialLibrary(lib)
    wmatlib.write_hdf5(filename)
    rmatlib = MaterialLibrary()
    rmatlib.from_hdf5(filename)
    os.remove(filename)
    # Round trip!
    rmatlib.write_hdf5(filename)
    wmatlib = MaterialLibrary(filename)
    assert_equal(set(wmatlib.keys()), set(rmatlib.keys()))
    for key in rmatlib.keys():
        assert_mat_almost_equal(wmatlib[key], rmatlib[key])
    os.remove(filename)

# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
