"""PyNE nuclear data tests"""
import os
import math
import warnings

import nose
from nose.tools import assert_equal, assert_in, assert_true
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

def test_atomic_mass():
    # set the datapath to nonsense so we call, the cpp data version
    pyne_conf.NUC_DATA_PATH = b'bobbobhonkeytonk'
    o16 = [15.99491461957, 16.0]
    u235 = [235.043930131, 235.0]
    am242m = [242.059549364, 242.0]

    # zzaam form
    assert_in(data.atomic_mass(80160), o16)
    assert_in(data.atomic_mass(922350), u235)
    assert_in(data.atomic_mass(952421), am242m)


def test_natural_abund_excited_state():
    # set the datapath to nonsense so we call, the cpp data version
    pyne_conf.NUC_DATA_PATH = b'bobbobhonkeytonk'
    # initialize natural_abund_map
    gnd = 902320000
    excited = gnd + 1
    data.natural_abund(gnd)
    # excited state should not be in the map yet
    assert_equal(data.natural_abund_map.get(excited), None)
    nabund = data.natural_abund(excited)
    assert_equal(nabund, data.natural_abund_map.get(excited))



if __name__ == "__main__":
    nose.runmodule()
