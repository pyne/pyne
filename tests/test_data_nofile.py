"""PyNE nuclear data tests"""
import os
import math
import warnings

import pytest
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
    # set the datapath to nonsense so we call  the cpp data version
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b"bobbobhonkeytonk"
    o16 = [15.9949146196, 16.0]
    u235 = [235.04392819, 235.0]
    am242m = [242.059547428, 242.0]

    # zzaam form
    assert data.atomic_mass(80160) in o16
    assert data.atomic_mass(922350) in u235
    assert data.atomic_mass(952421) in am242m
    pyne_conf.NUC_DATA_PATH = orig


def test_natural_abund_excited_state():
    # set the datapath to nonsense so we call, the cpp data version
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b"bobbobhonkeytonk"
    # initialize natural_abund_map
    gnd = 902320000
    excited = gnd + 42
    data.natural_abund(gnd)
    # excited state should not be in the map yet
    assert excited not in data.natural_abund_map
    nabund = data.natural_abund(excited)
    assert nabund == data.natural_abund_map.get(excited)
    pyne_conf.NUC_DATA_PATH = orig


def test_elements():
    # set the datapath to nonsense so we call, the cpp data version
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b"bobbobhonkeytonk"
    # initialize natural_abund_map
    # test a series of elements
    assert data.atomic_mass("H") == pytest.approx(1.007940754065775)
    assert data.atomic_mass("Li") == pytest.approx(6.940036602972684)
    assert data.atomic_mass("U") == pytest.approx(238.02890905398374)
    # NB any elements beyond z = 92 are not natural
    # and therefore have 0.0 atomic mass
    assert data.atomic_mass("Pu") == pytest.approx(0.0)
    # note if you use the nuc_data.h5 file it
    # has the same behaviour
    pyne_conf.NUC_DATA_PATH = orig

