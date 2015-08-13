#!/usr/bin/env python

from unittest import TestCase
import warnings
from nose.tools import assert_equal, assert_raises
from nose.plugins.skip import SkipTest
import numpy as np

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne.cccc import Isotxs, Rtflux

try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    from nose.plugins.skip import SkipTest
    HAVE_PYTAPS = False

class TestIsotxs(TestCase):
    def setUp(self):
        self.iso = Isotxs('ISOTXS')
        try:
            self.iso.read()
        except:
            raise SkipTest

    def test_isotxs_data(self):
        assert self.iso.emax[0] == 10000000.0
        assert self.iso.emax[4] == 0.625
        assert self.iso.emin == 9.999999747378752e-06
        assert self.iso.chi[0] == 0.7654485106468201
        assert self.iso.chi[3] == 0.0
        assert self.iso.vel[2] == 7650702.0
        assert self.iso.label[:6] == 'ISOTXS'

        assert self.iso.fc == {'ichidst': 1, 'maxdown': 6, 'maxord': 2,
                               'maxup': 6, 'ngroup': 7, 'niso': 39,
                               'nsblok': 1, 'nscmax': 4}

        assert len(self.iso.nucNames) == 39
        assert self.iso.nucNames[0] == 'H1gp 0'
        assert self.iso.nucNames[20] == 'MNbl 0'
        assert self.iso.nucNames[-1] == 'GPHl 0'

    def test_nuclide_data(self):
        assert len(self.iso.nuclides) == 39

        nuc = self.iso.nuclides[0]
        assert nuc.name == 'H1gp 0'
        assert nuc.libParams['libName'] == 'ENDF/B-6'
        assert nuc.libParams['ords'] == [2, 2, 0, 0]
        assert nuc.libParams['temp'] == 300.0
        assert nuc.libParams['amass'] == 1.0078279972076416
        assert nuc.libParams['fisFlag'] == 0
        assert nuc.libParams['jband'][3, 3] == 7

    def test_methods(self):
        nuc = self.iso.nuclides[20]
        nuc2 = self.iso.find_nuclide('MNbl 0')
        assert nuc == nuc2


def test_rtflux_basics():
    rt = Rtflux("files_test_cccc/rtflux_3D")
    assert_equal(rt.hname, "rtflux")
    assert_equal(rt.huse, "155339")
    assert_equal(rt.ivers,"081220")
    assert_equal(rt.ndim, 3)
    assert_equal(rt.ngroup, 4)
    assert_equal(rt.ninti, 4)
    assert_equal(rt.nintj, 4)
    assert_equal(rt.nintk, 4)
    assert_equal(rt.niter, 23)
    assert_equal(rt.effk, 0.892322838306427)
    assert_equal(rt.nblok, 1)
    assert_equal(rt.adjoint, False)

def test_rtflux_3D():

    if not HAVE_PYTAPS:
        raise SkipTest
    from pyne.mesh import Mesh, IMeshTag

    rt = Rtflux("files_test_cccc/rtflux_3D")
    structured_coords=[[0.0,  30.0, 40.0, 50.0, 70.0],
                       [0.0, 15.0, 40.0,  50.0, 70.0],
                       [0.0, 20.0, 75.0, 130.0, 150.0]]
    m = Mesh(structured=True, structured_coords=structured_coords)
    rt.to_mesh(m, "flux")
    m.tag = IMeshTag(4, float, name="flux")
    flux = m.tag[:]
    # test energy ordering
    np.allclose(flux[0], [2.34609601e-05, 1.20151503e-04, 
                          7.48262958e-05, 2.66320088e-06])
    # test spatial ordering
    exp = [2.34609601e-05, 6.06449525e-05, 6.06449525e-05, 2.34609601e-05, 1.93416323e-05,
           4.99965586e-05, 4.99965586e-05, 1.93416323e-05, 1.37740815e-05, 3.55703333e-05,
           3.55703333e-05, 1.37740815e-05, 6.93809615e-06, 1.78326363e-05, 1.78326363e-05,
           6.93809615e-06, 7.28244260e-06, 1.82284016e-05, 1.82284016e-05, 7.28244260e-06,
           6.00219150e-06, 1.50229626e-05, 1.50229626e-05, 6.00219150e-06, 4.28577852e-06,
           1.07128320e-05, 1.07128320e-05, 4.28577852e-06, 2.22228994e-06, 5.52725252e-06,
           5.52725252e-06, 2.22228994e-06, 2.09818192e-06, 4.99869927e-06, 4.99869927e-06,
           2.09818192e-06, 1.72912792e-06, 4.11851949e-06, 4.11851949e-06, 1.72912792e-06,
           1.24560384e-06, 2.96332430e-06, 2.96332430e-06, 1.24560384e-06, 6.68628056e-07,
           1.58531510e-06, 1.58531510e-06, 6.68628056e-07, 5.91257887e-07, 1.38688199e-06,
           1.38688199e-06, 5.91257887e-07, 4.88278620e-07, 1.14497300e-06, 1.14497300e-06,
           4.88278620e-07, 3.55008243e-07, 8.31776241e-07, 8.31776241e-07, 3.55008243e-07,
           1.91609508e-07, 4.46100170e-07, 4.46100170e-07, 1.91609508e-07]
    np.allclose(exp, flux[:,0])


#def test_rt_flux_1D():
#    if not HAVE_PYTAPS:
#        raise SkipTest
#    from pyne.mesh import Mesh, IMeshTag
#
#    rt = Rtflux("files_test_cccc/rtflux_1D")
#    structured_coords=[[10*x for x in range(8)], [0.0, 1.0], [0.0, 1.0]]
#    m = Mesh(structured=True, structured_coords=structured_coords)
#    rt.to_mesh(m, "flux")
#    m.tag = IMeshTag(4, float, name="flux")
#    flux = m.tag[:]

def test_atflux():
    rt = Rtflux("files_test_cccc/atflux_3D")
    assert_equal(rt.adjoint, True)
