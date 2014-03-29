#!/usr/bin/env python
from __future__ import unicode_literals
from unittest import TestCase

from nose.tools import assert_equal, assert_raises

from pyne.cccc import Isotxs

class TestIsotxs(TestCase):

    def setUp(self):
        self.iso = Isotxs('ISOTXS')
        self.iso.read()

    def test_isotxs_data(self):
        assert self.iso.emax[0] == 10000000.0
        assert self.iso.emax[4] == 0.625
        assert self.iso.emin == 9.999999747378752e-06
        assert self.iso.chi[0] == 0.7654485106468201
        assert self.iso.chi[3] == 0.0
        assert self.iso.vel[2] == 7650702.0
        assert self.iso.label[:6] == b'ISOTXS'

        assert self.iso.fc == {'ichidst': 1, 'maxdown': 6, 'maxord': 2,
                               'maxup': 6, 'ngroup': 7, 'niso': 39,
                               'nsblok': 1, 'nscmax': 4}

        assert len(self.iso.nucNames) == 39
        assert self.iso.nucNames[0] == b'H1gp 0'
        assert self.iso.nucNames[20] == b'MNbl 0'
        assert self.iso.nucNames[-1] == b'GPHl 0'

    def test_nuclide_data(self):
        assert len(self.iso.nuclides) == 39

        nuc = self.iso.nuclides[0]
        assert nuc.name == b'H1gp 0'
        assert nuc.libParams['libName'] == b'ENDF/B-6'
        assert nuc.libParams['ords'] == [2, 2, 0, 0]
        assert nuc.libParams['temp'] == 300.0
        assert nuc.libParams['amass'] == 1.0078279972076416
        assert nuc.libParams['fisFlag'] == 0
        assert nuc.libParams['jband'][3,3] == 7

    def test_methods(self):
        nuc = self.iso.nuclides[20]
        nuc2 = self.iso.find_nuclide(b'MNbl 0')
        assert nuc == nuc2
