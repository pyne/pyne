"""PyNE Material expand elements test under the presence of data or no data tests"""
import os
import math
import warnings

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

nucvec = {"H": 1.0, "Fe": 60.0, "Mn": 39.0}


def test_with_internal_data():
    orig = pyne_conf.NUC_DATA_PATH
    pyne_conf.NUC_DATA_PATH = b"thisisanonsensedatapath"

    mat = Material(nucvec, density=7.8)
    mat = mat.expand_elements()
    assert (id("He3") in mat.comp) == False
    assert (id("Co58") in mat.comp) == False
    assert (id("Ni58") in mat.comp) == False

    assert (id("H1") in mat.comp) == True
    assert (id("Fe56") in mat.comp) == True
    assert (id("Mn55") in mat.comp) == True

    pyne_conf.NUC_DATA_PATH = orig

