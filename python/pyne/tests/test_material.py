"""Material tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal

import os
from pyne.material import Material
import numpy  as np
import tables as tb


#def setup_class(TestMaterialConstructor):
#    with open('mat.txt', 'w') as f:
#        f.write('U235  0.05\nU238  0.95')


"""\
class TestMaterialConstructor(TestCase):
    "Tests that the Material constructors works."

    @classmethod
    def setup_class(cls):
        "Make temporary file for constructors to read in."
        with open('mat.txt', 'w') as f:
            f.write('U235  0.05\nU238  0.95')

    @classmethod
    def teardown_class(cls):
        "Remove temporary file so that we don't clutter up the filesystem."
        os.remove('mat.txt')

    def test_mat1(self):
        mat = Material("mat.txt")
        assert_equal(mat.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(mat.mass, 1.0)
        assert_equal(mat.name, '')

    def test_mat2(self):
        mat = Material("mat.txt", 42)
        assert_equal(mat.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(mat.mass, 42.0)
        assert_equal(mat.name, '')

    def test_mat3(self):
        mat = Material("mat.txt", -42, "My Stream")
        assert_equal(mat.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(mat.mass, 1.0)
        assert_equal(mat.name, 'My Stream')

    def test_mat4(self):
        mat = Material({922350: 0.05, 922380: 0.95}, 15, "Dict Try")
        assert_equal(mat.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(mat.mass, 15.0)
        assert_equal(mat.name, 'Dict Try')

    def test_load_from_hdf5(self):
        #First make a temp file
        f = tb.openFile("mat.h5", "w")
        f.createGroup("/", "mat", "Mass Stream Test")
        f.createArray("/mat", "Mass",  np.array([1.0, 0.5,  0.0]), "Mass Test")
        f.createArray("/mat", "U235",  np.array([1.0, 0.75, 0.0]), "U235 Test")
        f.createArray("/mat", "PU239", np.array([0.0, 0.25, 0.0]), "PU239 Test")
        f.close()

        #perform tests
        mat = Material()
        mat.load_from_hdf5("mat.h5", "/mat")
        assert_equal(mat.mass, 0.0)
        assert_equal(mat.comp, {922350: 0.0, 942390: 0.0})

        mat.load_from_hdf5("mat.h5", "/mat", 0)
        assert_equal(mat.mass, 1.0)
        assert_equal(mat.comp, {922350: 1.0, 942390: 0.0})

        mat.load_from_hdf5("mat.h5", "/mat", 1)
        assert_equal(mat.mass, 0.5)
        assert_equal(mat.comp, {922350: 0.75, 942390: 0.25})

        mat.load_from_hdf5("mat.h5", "/mat", 2)
        assert_equal(mat.mass, 0.0)
        assert_equal(mat.comp, {922350: 0.0, 942390: 0.0})

        mat.load_from_hdf5("mat.h5", "/mat", -1)
        assert_equal(mat.mass, 0.0)
        assert_equal(mat.comp, {922350: 0.0, 942390: 0.0})

        mat.load_from_hdf5("mat.h5", "/mat", -2)
        assert_equal(mat.mass, 0.5)
        assert_equal(mat.comp, {922350: 0.75, 942390: 0.25})

        mat.load_from_hdf5("mat.h5", "/mat", -3)
        assert_equal(mat.mass, 1.0)
        assert_equal(mat.comp, {922350: 1.0, 942390: 0.0})

        #clean up
        os.remove('mat.h5')

    def test_load_from_text(self):
        mat = Material()
        mat.load_from_text("mat.txt")
        assert_equal(mat.comp, {922350: 0.05, 922380: 0.95})
"""



class TestMaterialMethods(TestCase):
    "Tests that the Material member functions work."

    #def test_normalize(self):
    #    mat = Material({922350: 0.05, 922380: 0.95}, 15)
    #    mat.normalize()
    #    assert_equal(mat.mass, 1.0)


    def test_mult_by_mass(self):
        mat = Material({922350: 0.05, 922380: 0.95}, 15)
        nucvec = mat.mult_by_mass()
        assert_equal(nucvec, {922350: 0.75, 922380: 14.25})

"""\

    def test_atomic_weight(self):
        mat_empty = Material({})
        assert_equal(mat_empty.atomic_weight(), 0.0)

        mat_u238 = Material({922380: 1.0})
        assert_equal(mat_u238.atomic_weight(), 238.0)

        mat_mixed = Material({922350: 0.5, 922380: 0.5})
        assert_almost_equal(mat_mixed.atomic_weight()/236.5, 1.0, 4)



class TestMassSubStreamMethods(TestCase):
    "Tests that the Material sub-stream getter member functions work."

    isovec = {
        10010:  1.0,   
        80160:  1.0,   
        691690: 1.0,
        922350: 1.0,
        922380: 1.0,
        942390: 1.0,
        942410: 1.0,
        952420: 1.0,
        962440: 1.0,
        }

    def test_get_sub_streamInt_1(self):
        mat = Material(self.isovec, -1, "Old Stream")
        mat1 = mat.sub_stream([92, 80160])
        assert_almost_equal(mat1.comp[80160],  0.3333333333333)
        assert_almost_equal(mat1.comp[922350], 0.3333333333333)
        assert_almost_equal(mat1.comp[922380], 0.3333333333333)
        assert_equal(mat1.mass, 3.0)
        assert_equal(mat1.name, '')

    def test_get_sub_streamInt_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_stream([92, 80160], "New Stream")
        assert_almost_equal(mat1.comp[80160],  0.3333333333333)
        assert_almost_equal(mat1.comp[922350], 0.3333333333333)
        assert_almost_equal(mat1.comp[922380], 0.3333333333333)
        assert_equal(mat1.mass, 3.0)
        assert_equal(mat1.name, 'New Stream')

    def test_get_sub_streamattr_1(self):
        mat = Material(self.isovec, -1, "Old Stream")
        mat1 = mat.sub_stream(["U", "80160", "H1"])
        assert_almost_equal(mat1.comp[10010],  0.25)
        assert_almost_equal(mat1.comp[80160],  0.25)
        assert_almost_equal(mat1.comp[922350], 0.25)
        assert_almost_equal(mat1.comp[922380], 0.25)
        assert_equal(mat1.mass, 4.0)
        assert_equal(mat1.name, '')

    def test_get_sub_streamattr_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_stream(["U", "80160", "H1"], "New Stream")
        assert_almost_equal(mat1.comp[10010],  0.25)
        assert_almost_equal(mat1.comp[80160],  0.25)
        assert_almost_equal(mat1.comp[922350], 0.25)
        assert_almost_equal(mat1.comp[922380], 0.25)
        assert_equal(mat1.mass, 4.0)
        assert_equal(mat1.name, 'New Stream')

    def test_get_u_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_u()
        assert_equal(mat1.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, '')

    def test_get_u_2(self):
        mat = Material(self.isovec)
        mat1 = mat.get_u("U Stream")
        assert_equal(mat1.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, 'U Stream')

    def test_get_pu_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_pu()
        assert_equal(mat1.comp, {942390: 0.5, 942410: 0.5})
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, '')

    def test_get_pu_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_pu("PU Stream")
        assert_equal(mat1.comp, {942390: 0.5, 942410: 0.5})
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, 'PU Stream')

    def test_get_lan_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_lan()
        assert_equal(mat1.comp, {691690: 1.0})
        assert_equal(mat1.mass, 1.0)
        assert_equal(mat1.name, '')

    def test_get_lan_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_lan("LAN Stream")
        assert_equal(mat1.comp, {691690: 1.0})
        assert_equal(mat1.mass, 1.0)
        assert_equal(mat1.name, 'LAN Stream')

    def test_get_act_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_act()
        assert_equal(mat1.comp[922350], 1.0/6.0)
        assert_equal(mat1.comp[922380], 1.0/6.0)
        assert_equal(mat1.comp[942390], 1.0/6.0)
        assert_equal(mat1.comp[942410], 1.0/6.0)
        assert_equal(mat1.comp[952420], 1.0/6.0)
        assert_equal(mat1.comp[962440], 1.0/6.0)
        assert_equal(mat1.mass, 6.0)
        assert_equal(mat1.name, '')

    def test_get_act_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_act("ACT Stream")
        assert_equal(mat1.comp[922350], 1.0/6.0)
        assert_equal(mat1.comp[922380], 1.0/6.0)
        assert_equal(mat1.comp[942390], 1.0/6.0)
        assert_equal(mat1.comp[942410], 1.0/6.0)
        assert_equal(mat1.comp[952420], 1.0/6.0)
        assert_equal(mat1.comp[962440], 1.0/6.0)
        assert_equal(mat1.mass, 6.0)
        assert_equal(mat1.name, 'ACT Stream')

    def test_get_tru_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_tru()
        assert_equal(mat1.comp[942390], 1.0/4.0)
        assert_equal(mat1.comp[942410], 1.0/4.0)
        assert_equal(mat1.comp[952420], 1.0/4.0)
        assert_equal(mat1.comp[962440], 1.0/4.0)
        assert_equal(mat1.mass, 4.0)
        assert_equal(mat1.name, '')

    def test_get_tru_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_tru("TRU Stream")
        assert_equal(mat1.comp[942390], 1.0/4.0)
        assert_equal(mat1.comp[942410], 1.0/4.0)
        assert_equal(mat1.comp[952420], 1.0/4.0)
        assert_equal(mat1.comp[962440], 1.0/4.0)
        assert_equal(mat1.mass, 4.0)
        assert_equal(mat1.name, 'TRU Stream')

    def test_get_ma_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_ma()
        assert_equal(mat1.comp[952420], 1.0/2.0)
        assert_equal(mat1.comp[962440], 1.0/2.0)
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, '')

    def test_get_ma_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_ma("MA Stream")
        assert_equal(mat1.comp[952420], 1.0/2.0)
        assert_equal(mat1.comp[962440], 1.0/2.0)
        assert_equal(mat1.mass, 2.0)
        assert_equal(mat1.name, 'MA Stream')

    def test_get_fp_1(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_fp()
        assert_equal(mat1.comp[10010],  1.0/3.0)
        assert_equal(mat1.comp[80160],  1.0/3.0)
        assert_equal(mat1.comp[691690], 1.0/3.0)
        assert_equal(mat1.mass, 3.0)
        assert_equal(mat1.name, '')

    def test_get_fp_2(self):
        mat = Material(self.isovec)
        mat1 = mat.sub_fp("FP Stream")
        assert_equal(mat1.comp[10010],  1.0/3.0)
        assert_equal(mat1.comp[80160],  1.0/3.0)
        assert_equal(mat1.comp[691690], 1.0/3.0)
        assert_equal(mat1.mass, 3.0)
        assert_equal(mat1.name, 'FP Stream')

        
class TestMaterialOperatorOverloading(TestCase):
    "Tests that the Material operator overloads work."
    u235 = Material({922350: 1.0})
    u238 = Material({922380: 1.0})

    def test_add_num(self):
        mat = self.u235 + 30.0
        assert_equal(mat.mass, 31.0)

    def test_radd_num(self):
        mat = 90 + self.u235
        assert_equal(mat.mass, 91.0)

    def test_add_mat(self):
        mat = self.u235 + self.u238
        assert_equal(mat.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(mat.mass, 2.0)
        assert_equal(mat.name, '')

    def test_mul_num(self):
        mat = self.u235 * 2.0
        assert_equal(mat.mass, 2.0)

    def test_rmul_num(self):
        mat = 150 * self.u235
        assert_equal(mat.mass, 150.0)

    def test_div_num(self):
        mat = self.u235 / 10
        assert_equal(mat.mass, 0.1)
"""\


if __name__ == "__main__":
    nose.main()
