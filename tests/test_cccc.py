#!/usr/bin/env python
from unittest import TestCase
import warnings
from nose.tools import assert_equal, assert_raises
from nose.plugins.skip import SkipTest
import numpy as np
from numpy.testing import assert_array_almost_equal

from pyne.mesh import HAVE_PYMOAB
from pyne.cccc import Isotxs, Rtflux, Atflux
from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)


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
    assert_equal(rt.ivers, "081220")
    assert_equal(rt.ndim, 3)
    assert_equal(rt.ngroup, 4)
    assert_equal(rt.ninti, 4)
    assert_equal(rt.nintj, 4)
    assert_equal(rt.nintk, 4)
    assert_equal(rt.niter, 23)
    assert_equal(rt.effk, 0.892322838306427)
    assert_equal(rt.power, 1.0)
    assert_equal(rt.nblok, 1)
    assert_equal(rt.adjoint, False)


def test_rtflux_3D():

    if not HAVE_PYMOAB:
        raise SkipTest
    from pyne.mesh import Mesh, NativeMeshTag

    rt = Rtflux("files_test_cccc/rtflux_3D")
    structured_coords = [[0.0,  30.0, 40.0, 50.0, 70.0],
                         [0.0, 15.0, 40.0,  50.0, 70.0],
                         [0.0, 20.0, 75.0, 130.0, 150.0]]
    m = Mesh(structured=True, structured_coords=structured_coords)
    rt.to_mesh(m, "flux")
    m.tag = NativeMeshTag(4, float, name="flux")
    flux = m.tag[:]
    # test energy ordering
    assert(np.allclose(flux[0], [2.66320088e-06, 7.48262958e-05,
                                 1.20151503e-04, 2.34609601e-05]))
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
    assert(np.allclose(exp, flux[:, 3]))


def test_rtflux_2D():
    if not HAVE_PYMOAB:
        raise SkipTest
    from pyne.mesh import Mesh, NativeMeshTag

    rt = Rtflux("files_test_cccc/rtflux_2D")
    structured_coords = [
        [10*x for x in range(8)], [0.0, 10, 20, 30, 40], [0.0, 1.0]]
    m = Mesh(structured=True, structured_coords=structured_coords)
    rt.to_mesh(m, "flux")
    m.tag = NativeMeshTag(4, float, name="flux")
    flux = m.tag[:]
    # test energy ordering
    assert(np.allclose(flux[0], [1.66764798e-03, 4.59426961e-02,
                                 7.35252284e-02, 1.54202809e-02]))
    # test spatial ordering
    exp = [1.54202809e-02, 1.22833140e-02, 8.24652761e-03, 4.10247328e-03, 1.29812236e-02,
           7.31613464e-03, 4.27769488e-03, 2.96312777e-03, 9.98971577e-03, 4.57188750e-03,
           2.15025301e-03, 1.86188954e-03, 7.46111984e-03, 3.45483912e-03, 1.54143733e-03,
           1.24946029e-03, 5.24386070e-03, 2.36004487e-03, 1.01537828e-03, 8.43692879e-04,
           3.31488925e-03, 1.52156795e-03, 6.57790710e-04, 5.31460286e-04, 1.52433295e-03,
           7.17349305e-04, 3.12441756e-04, 2.43518040e-04]
    assert(np.allclose(exp, flux[:, 3]))


def test_rt_flux_1D():
    if not HAVE_PYMOAB:
        raise SkipTest
    from pyne.mesh import Mesh, NativeMeshTag

    rt = Rtflux("files_test_cccc/rtflux_1D")
    structured_coords = [[10*x for x in range(8)], [0.0, 1.0], [0.0, 1.0]]
    m = Mesh(structured=True, structured_coords=structured_coords)
    rt.to_mesh(m, "flux")
    m.tag = NativeMeshTag(4, float, name="flux")
    flux = m.tag[:]
    exp = [[1.13102481e-03, 2.48423595e-02, 4.07499865e-02, 1.12382315e-02],
           [1.22305619e-03, 2.06497203e-02, 2.49217686e-02, 4.70910079e-03],
           [9.99735815e-04, 1.35176066e-02, 1.08674872e-02, 1.10359953e-03],
           [6.18705350e-04, 7.42586666e-03, 4.91338778e-03, 5.60093921e-04],
           [3.40591572e-04, 3.81165736e-03, 2.18964779e-03, 1.75256136e-04],
           [1.65193192e-04, 1.78217782e-03, 9.88094844e-04, 1.16255047e-04],
           [5.67537762e-05, 6.01543637e-04, 3.08311260e-04, 1.42177473e-05]]
    assert(np.allclose(exp, flux))


def test_rtflux_raises():
    if not HAVE_PYMOAB:
        raise SkipTest
    from pyne.mesh import Mesh, NativeMeshTag

    rt = Rtflux("files_test_cccc/rtflux_1D")
    structured_coords = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
    m = Mesh(structured=True, structured_coords=structured_coords)
    with assert_raises(ValueError):
        rt.to_mesh(m, "flux")


def test_atflux_adjoint():
    at = Atflux("files_test_cccc/atflux_3D")
    assert_equal(at.adjoint, True)


def test_atflux_eng_order():
    """Ensure the energy order is read in reverse for atflux file."""
    if not HAVE_PYMOAB:
        raise SkipTest
    from pyne.mesh import Mesh, NativeMeshTag

    # This file is created with: source=1 174R 0 0 1 40R 0
    # So only the 2 highest energy photon groups and the 1 highest energy
    # neutron group should have non-zero fluxes.
    at = Atflux("files_test_cccc/atflux_eng_order")
    sc = [np.linspace(-100, 100, 5),
          np.linspace(-100, 100, 5),
          np.linspace(0, 330, 5)]
    m = Mesh(structured=True, structured_coords=sc, mats=None)
    at.to_mesh(m, "flux")
    m.flux = NativeMeshTag(217, float)
    assert_array_almost_equal(m.flux[0],
                              np.array([0]*40 + [57.3204927667, 1.16690395827] + [0]*174 + [14.2312186922]))
