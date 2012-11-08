"""Transmute tests"""

import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false

from numpy.testing import dec

import os
from pyne.material import Material
import numpy  as np
import tables as tb

from pyne import transmute

@dec.skipif(True)
def test_decay1():
    mat = Material({'C14': 1.0, 'N14': 0.0, 'C13': 0.0})
    obs = transmute.decay(mat, 1.0)

    assert False
