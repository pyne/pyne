"""PyNE Material expand elements test under the presence of data or no data tests"""
import os
import math
import warnings

import nose
from nose.tools import assert_equal, assert_in, assert_true, assert_almost_equal
import numpy as np
import numpy.testing as npt

from pyne.utils import QAWarning
from pyne.pyne_config import pyne_conf

warnings.simplefilter("ignore", QAWarning)

import pyne
from pyne import data, nucname
from pyne import utils

if utils.use_warnings():
    utils.toggle_warnings()

nucvec = {'H':  1.0,
          'Fe': 60.0,
          'Mn': 29.0
          }

def test_with_no_data():
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b'thisisanonsensedatapath'

    mat = Material(nucvec,density=7.8)
    assert_equal(mat.comp.has_key('He3'),False)
    assert_equal(mat.comp.has_key('Co58'),False)
    assert_equal(mat.comp.has_key('H1'),True)
    assert_equal(mat.comp.has_key('Mn55'),True)
    assert_equal(mat.comp.has_key('Fe56'),True)

    pyne_conf.NUC_DATA_PATH = orig

def test_with_data():
    mat = Material(nucvec,density=7.8)
    assert_equal(mat.comp.has_key('He3'),False)
    assert_equal(mat.comp.has_key('Co58'),False)
    assert_equal(mat.comp.has_key('H1'),True)
    assert_equal(mat.comp.has_key('Mn55'),True)
    assert_equal(mat.comp.has_key('Fe56'),True)

# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
