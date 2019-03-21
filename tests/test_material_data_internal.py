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
from pyne.material import Material
from pyne import data, nucname
from pyne import utils

from pyne.nucname import id

if utils.use_warnings():
    utils.toggle_warnings()

nucvec = {'H':  1.0,
          'Fe': 60.0,
          'Mn': 39.0
          }

def test_with_internal_data():
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b'thisisanonsensedatapath'

    mat = Material(nucvec,density=7.8)
    mat = mat.expand_elements()
    assert_equal(id('He3') in mat.comp,False)
    assert_equal(id('Co58') in mat.comp,False)
    assert_equal(id('Ni58') in mat.comp,False)

    assert_equal(id('H1') in mat.comp,True)
    assert_equal(id('Fe56') in mat.comp,True)
    assert_equal(id('Mn55') in mat.comp,True)
   


# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
