"""PyNE decay tests and benchmarks."""
from __future__ import print_function, unicode_literals
import os
import sys
import warnings
from functools import lru_cache
if sys.version_info[0] >= 3:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve

import nose
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in, assert_true

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import nucname
from pyne import data
from pyne import origen22
from pyne.material import Material, MaterialLibrary

# import decay gen
srcdir = os.path.join(os.path.dirname(__file__), '..', 'src')
srcdir = os.path.abspath(srcdir)
sys.path.insert(0, srcdir)
import decaygen

MATS = None

def setup():
    global MATS
    o2benchfile = 'origen-benchmark.h5'
    o2benchurl = 'http://data.pyne.io/origen-benchmark.h5'
    if not os.path.exists(o2benchfile):
        sys.stderr.write("\nDownloading " + o2benchurl + ": ")
        sys.stderr.flush()
        urlretrieve(o2benchurl, o2benchfile)
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    MATS = MaterialLibrary(o2benchfile)


#
# helper functions
#

def t9_half_life(nuc):
    nuczz = nucname.zzaaam(nuc)
    hl = None
    for i in [1, 2, 3]:
        hl = t9[i]['half_life'].get(nuczz, None)
        if hl is not None:
            break
    return hl

@lru_cache()
def hl_relerr(nuc):
    """Half-life relative error."""
    dhl = data.half_life(nuc) or np.inf
    ohl = t9_half_life(nuc) or np.inf
    if np.isinf(ohl) and np.isinf(dhl):
        return 0.0
    hlre = np.abs(ohl - dhl) * 2 / (ohl + dhl)
    return np.inf if np.isnan(hlre) else hlre

if __name__ == "__main__":
    nose.runmodule()
