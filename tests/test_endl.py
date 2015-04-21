import os
import io
import warnings
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from math import e

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose, \
    assert_array_almost_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)

from pyne.endl import Library

import nose
from nose.tools import assert_equal

def ignore_future_warnings(func):
    """This is a decorator which can be used to ignore FutureWarnings
    occurring in a function."""
    def new_func(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def test_loadfile():
    try:
        assert(os.path.isfile('epdl97_eedl_Pb'))
    except AssertionError:
        try:
            import urllib.request as urllib
        except ImportError:
            import urllib
        urllib.urlretrieve("https://www-nds.iaea.org/epdl97/data/endl/eedl/za082000",
                    "epdl97_eedl_Pb")
    from hashlib import md5
    with open("epdl97_eedl_Pb", "rb") as f:
        obs_hash = md5(f.read()).hexdigest()
    exp_hash = "502105669e0c0ad917301d219c649aaf"
    try:
        assert_equal(obs_hash, exp_hash)
    except AssertionError:
        raise AssertionError("epdl97_eedl_Pb hash check failed; please try redownloading the epdl97_eedl_Pb data file.")
    testlib = Library("epdl97_eedl_Pb")
#    exp_nuclides = [10010000, 30060000, 40090000, 50110000, 30070000, 60000000, 50100000]
#    obs_nuclides = list(map(int, testlib.structure.keys()))
#    assert_array_equal(exp_nuclides, obs_nuclides)

if __name__ == "__main__":
    nose.runmodule()
