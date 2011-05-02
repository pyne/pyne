"""Mass stream tests"""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal

import os
import mass_stream
import numpy  as np
import tables as tb

#def setup_class(TestMassStreamConstructor):
#    with open('ms.txt', 'w') as f:
#        f.write('U235  0.05\nU238  0.95')


class TestMassStreamConstructor(TestCase):
    """Tests that the MassStream constructors works."""

    @classmethod
    def setup_class(cls):
        """Make temporary file for constructors to read in."""
        with open('ms.txt', 'w') as f:
            f.write('U235  0.05\nU238  0.95')

    @classmethod
    def teardown_class(cls):
        """Remove temporary file so that we don't clutter up the filesystem."""
        os.remove('ms.txt')

    def test_ms1(self):
        ms = mass_stream.MassStream("ms.txt")
        assert_equal(ms.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(ms.mass, 1.0)
        assert_equal(ms.name, '')

    def test_ms2(self):
        ms = mass_stream.MassStream("ms.txt", 42)
        assert_equal(ms.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(ms.mass, 42.0)
        assert_equal(ms.name, '')

    def test_ms3(self):
        ms = mass_stream.MassStream("ms.txt", -42, "My Stream")
        assert_equal(ms.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(ms.mass, 1.0)
        assert_equal(ms.name, 'My Stream')

    def test_ms4(self):
        ms = mass_stream.MassStream({922350: 0.05, 922380: 0.95}, 15, "Dict Try")
        assert_equal(ms.comp, {922350: 0.05, 922380: 0.95})
        assert_equal(ms.mass, 15.0)
        assert_equal(ms.name, 'Dict Try')

    def test_load_from_hdf5(self):
        #First make a temp file
        f = tb.openFile("ms.h5", "w")
        f.createGroup("/", "ms", "Mass Stream Test")
        f.createArray("/ms", "Mass",  np.array([1.0, 0.5,  0.0]), "Mass Test")
        f.createArray("/ms", "U235",  np.array([1.0, 0.75, 0.0]), "U235 Test")
        f.createArray("/ms", "PU239", np.array([0.0, 0.25, 0.0]), "PU239 Test")
        f.close()

        #perform tests
        ms = mass_stream.MassStream()
        ms.load_from_hdf5("ms.h5", "/ms")
        assert_equal(ms.mass, 0.0)
        assert_equal(ms.comp, {922350: 0.0, 942390: 0.0})

        ms.load_from_hdf5("ms.h5", "/ms", 0)
        assert_equal(ms.mass, 1.0)
        assert_equal(ms.comp, {922350: 1.0, 942390: 0.0})

        ms.load_from_hdf5("ms.h5", "/ms", 1)
        assert_equal(ms.mass, 0.5)
        assert_equal(ms.comp, {922350: 0.75, 942390: 0.25})

        ms.load_from_hdf5("ms.h5", "/ms", 2)
        assert_equal(ms.mass, 0.0)
        assert_equal(ms.comp, {922350: 0.0, 942390: 0.0})

        ms.load_from_hdf5("ms.h5", "/ms", -1)
        assert_equal(ms.mass, 0.0)
        assert_equal(ms.comp, {922350: 0.0, 942390: 0.0})

        ms.load_from_hdf5("ms.h5", "/ms", -2)
        assert_equal(ms.mass, 0.5)
        assert_equal(ms.comp, {922350: 0.75, 942390: 0.25})

        ms.load_from_hdf5("ms.h5", "/ms", -3)
        assert_equal(ms.mass, 1.0)
        assert_equal(ms.comp, {922350: 1.0, 942390: 0.0})

        #clean up
        os.remove('ms.h5')

    def test_load_from_text(self):
        ms = mass_stream.MassStream()
        ms.load_from_text("ms.txt")
        assert_equal(ms.comp, {922350: 0.05, 922380: 0.95})


class TestMassStreamMethods(TestCase):
    """Tests that the MassStream member functions work."""

    def test_normalize(self):
        ms = mass_stream.MassStream({922350: 0.05, 922380: 0.95}, 15)
        ms.normalize()
        assert_equal(ms.mass, 1.0)


    def test_mult_by_mass(self):
        ms = mass_stream.MassStream({922350: 0.05, 922380: 0.95}, 15)
        isovec = ms.mult_by_mass()
        assert_equal(isovec, {922350: 0.75, 922380: 14.25})


    def test_atomic_weight(self):
        ms_empty = mass_stream.MassStream({})
        assert_equal(ms_empty.atomic_weight(), 0.0)

        ms_u238 = mass_stream.MassStream({922380: 1.0})
        assert_equal(ms_u238.atomic_weight(), 238.0)

        ms_mixed = mass_stream.MassStream({922350: 0.5, 922380: 0.5})
        assert_almost_equal(ms_mixed.atomic_weight()/236.5, 1.0, 4)


class TestMassSubStreamMethods(TestCase):
    """Tests that the MassStream sub-stream getter member functions work."""

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
        ms = mass_stream.MassStream(self.isovec, -1, "Old Stream")
        ms1 = ms.get_sub_stream([92, 80160])
        assert_almost_equal(ms1.comp[80160],  0.3333333333333)
        assert_almost_equal(ms1.comp[922350], 0.3333333333333)
        assert_almost_equal(ms1.comp[922380], 0.3333333333333)
        assert_equal(ms1.mass, 3.0)
        assert_equal(ms1.name, '')

    def test_get_sub_streamInt_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_sub_stream([92, 80160], "New Stream")
        assert_almost_equal(ms1.comp[80160],  0.3333333333333)
        assert_almost_equal(ms1.comp[922350], 0.3333333333333)
        assert_almost_equal(ms1.comp[922380], 0.3333333333333)
        assert_equal(ms1.mass, 3.0)
        assert_equal(ms1.name, 'New Stream')

    def test_get_sub_streamStr_1(self):
        ms = mass_stream.MassStream(self.isovec, -1, "Old Stream")
        ms1 = ms.get_sub_stream(["U", "80160", "H1"])
        assert_almost_equal(ms1.comp[10010],  0.25)
        assert_almost_equal(ms1.comp[80160],  0.25)
        assert_almost_equal(ms1.comp[922350], 0.25)
        assert_almost_equal(ms1.comp[922380], 0.25)
        assert_equal(ms1.mass, 4.0)
        assert_equal(ms1.name, '')

    def test_get_sub_streamStr_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_sub_stream(["U", "80160", "H1"], "New Stream")
        assert_almost_equal(ms1.comp[10010],  0.25)
        assert_almost_equal(ms1.comp[80160],  0.25)
        assert_almost_equal(ms1.comp[922350], 0.25)
        assert_almost_equal(ms1.comp[922380], 0.25)
        assert_equal(ms1.mass, 4.0)
        assert_equal(ms1.name, 'New Stream')

    def test_get_u_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_u()
        assert_equal(ms1.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, '')

    def test_get_u_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_u("U Stream")
        assert_equal(ms1.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, 'U Stream')

    def test_get_pu_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_pu()
        assert_equal(ms1.comp, {942390: 0.5, 942410: 0.5})
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, '')

    def test_get_pu_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_pu("PU Stream")
        assert_equal(ms1.comp, {942390: 0.5, 942410: 0.5})
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, 'PU Stream')

    def test_get_lan_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_lan()
        assert_equal(ms1.comp, {691690: 1.0})
        assert_equal(ms1.mass, 1.0)
        assert_equal(ms1.name, '')

    def test_get_lan_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_lan("LAN Stream")
        assert_equal(ms1.comp, {691690: 1.0})
        assert_equal(ms1.mass, 1.0)
        assert_equal(ms1.name, 'LAN Stream')

    def test_get_act_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_act()
        assert_equal(ms1.comp[922350], 1.0/6.0)
        assert_equal(ms1.comp[922380], 1.0/6.0)
        assert_equal(ms1.comp[942390], 1.0/6.0)
        assert_equal(ms1.comp[942410], 1.0/6.0)
        assert_equal(ms1.comp[952420], 1.0/6.0)
        assert_equal(ms1.comp[962440], 1.0/6.0)
        assert_equal(ms1.mass, 6.0)
        assert_equal(ms1.name, '')

    def test_get_act_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_act("ACT Stream")
        assert_equal(ms1.comp[922350], 1.0/6.0)
        assert_equal(ms1.comp[922380], 1.0/6.0)
        assert_equal(ms1.comp[942390], 1.0/6.0)
        assert_equal(ms1.comp[942410], 1.0/6.0)
        assert_equal(ms1.comp[952420], 1.0/6.0)
        assert_equal(ms1.comp[962440], 1.0/6.0)
        assert_equal(ms1.mass, 6.0)
        assert_equal(ms1.name, 'ACT Stream')

    def test_get_tru_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_tru()
        assert_equal(ms1.comp[942390], 1.0/4.0)
        assert_equal(ms1.comp[942410], 1.0/4.0)
        assert_equal(ms1.comp[952420], 1.0/4.0)
        assert_equal(ms1.comp[962440], 1.0/4.0)
        assert_equal(ms1.mass, 4.0)
        assert_equal(ms1.name, '')

    def test_get_tru_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_tru("TRU Stream")
        assert_equal(ms1.comp[942390], 1.0/4.0)
        assert_equal(ms1.comp[942410], 1.0/4.0)
        assert_equal(ms1.comp[952420], 1.0/4.0)
        assert_equal(ms1.comp[962440], 1.0/4.0)
        assert_equal(ms1.mass, 4.0)
        assert_equal(ms1.name, 'TRU Stream')

    def test_get_ma_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_ma()
        assert_equal(ms1.comp[952420], 1.0/2.0)
        assert_equal(ms1.comp[962440], 1.0/2.0)
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, '')

    def test_get_ma_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_ma("MA Stream")
        assert_equal(ms1.comp[952420], 1.0/2.0)
        assert_equal(ms1.comp[962440], 1.0/2.0)
        assert_equal(ms1.mass, 2.0)
        assert_equal(ms1.name, 'MA Stream')

    def test_get_fp_1(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_fp()
        assert_equal(ms1.comp[10010],  1.0/3.0)
        assert_equal(ms1.comp[80160],  1.0/3.0)
        assert_equal(ms1.comp[691690], 1.0/3.0)
        assert_equal(ms1.mass, 3.0)
        assert_equal(ms1.name, '')

    def test_get_fp_2(self):
        ms = mass_stream.MassStream(self.isovec)
        ms1 = ms.get_fp("FP Stream")
        assert_equal(ms1.comp[10010],  1.0/3.0)
        assert_equal(ms1.comp[80160],  1.0/3.0)
        assert_equal(ms1.comp[691690], 1.0/3.0)
        assert_equal(ms1.mass, 3.0)
        assert_equal(ms1.name, 'FP Stream')

        
class TestMassStreamOperatorOverloading(TestCase):
    """Tests that the MassStream operator overloads work."""
    u235 = mass_stream.MassStream({922350: 1.0})
    u238 = mass_stream.MassStream({922380: 1.0})

    def test_add_num(self):
        ms = self.u235 + 30.0
        assert_equal(ms.mass, 31.0)

    def test_radd_num(self):
        ms = 90 + self.u235
        assert_equal(ms.mass, 91.0)

    def test_add_ms(self):
        ms = self.u235 + self.u238
        assert_equal(ms.comp, {922350: 0.5, 922380: 0.5})
        assert_equal(ms.mass, 2.0)
        assert_equal(ms.name, '')

    def test_mul_num(self):
        ms = self.u235 * 2.0
        assert_equal(ms.mass, 2.0)

    def test_rmul_num(self):
        ms = 150 * self.u235
        assert_equal(ms.mass, 150.0)

    def test_div_num(self):
        ms = self.u235 / 10
        assert_equal(ms.mass, 0.1)

if __name__ == "__main__":
    nose.main()
