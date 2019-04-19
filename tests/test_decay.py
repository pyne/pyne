"""PyNE decay tests and benchmarks."""
from __future__ import print_function, unicode_literals
import os
import sys
import json
import warnings
if sys.version_info[0] >= 3:
    from urllib.request import urlretrieve
    from functools import lru_cache
else:
    from urllib import urlretrieve
    lru_cache = lambda *args, **kwargs: (lambda f: f)

import nose
from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_in, assert_true, assert_less

import numpy as np
import tables as tb
from numpy.testing import assert_array_equal, assert_array_almost_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne import nucname
from pyne import data
from pyne import origen22
from pyne.material import Material
from pyne.material_library import MaterialLibrary

# import decay gen
srcdir = os.path.join(os.path.dirname(__file__), '..', 'src')
srcdir = os.path.abspath(srcdir)
sys.path.insert(0, srcdir)
import decaygen

h5ver = tuple(map(int, tb.hdf5_version.split('-', 1)[0].split('.')))
if h5ver == (1, 8, 13):
    H5NAME = 'origen-benchmark-hdf5-1.8.13.h5'
else:
    H5NAME = 'origen-benchmark-hdf5-1.8.14.h5'
MATS = None
O2HLS = None  # Origen Half-lives

def setup():
    global MATS, O2HLS
    o2benchurl = 'https://github.com/pyne/data/raw/master/' + H5NAME
    if not os.path.exists(H5NAME):
        sys.stderr.write("\nDownloading " + o2benchurl + ": ")
        sys.stderr.flush()
        urlretrieve(o2benchurl, H5NAME)
        sys.stderr.write("done.\n")
        sys.stderr.flush()
    MATS = MaterialLibrary(H5NAME)
    with open('o2hls.json', 'r') as f:
        O2HLS = json.load(f)
    O2HLS = {int(nuc): v for nuc, v in O2HLS.items()}


METASTABLE_BLACKLIST = {
    771940001,  # have 2+ metastables that ORIGEN lumps together, or
    340830001,  # ...
    350830000,  # decay to metastables without reason.
    340830000,  # ...
    481170000,  # ...
    501130001,  # missing branch in origen
    }


#
# helper functions
#

def t9_half_life(nuc):
    return O2HLS.get(nuc, None)


@lru_cache()
def hl_relerr(nuc):
    """Half-life relative error."""
    dhl = data.half_life(nuc) or np.inf
    ohl = t9_half_life(nuc) or np.inf
    if np.isinf(ohl) and np.isinf(dhl):
        return 0.0
    hlre = np.abs(ohl - dhl) * 2 / (ohl + dhl)
    return np.inf if np.isnan(hlre) else hlre


def pivot_mat_keys():
    """This puts the material keys into a dict mapping nuclides to the data
    set names (eg 'origen_922350000_3'). This filters out any data sets that
    contain high relative errors in the half-lives. It also filters out
    species for which Origen has weird metastable behavior that seems
    unphysical.
    """
    nuc_keys = {}
    for key in MATS.keys():
        key = str(key, 'UTF-8')
        _, nuc, t = key.split('_')
        nuc = int(nuc)
        if nuc in METASTABLE_BLACKLIST:
            continue
        chains = decaygen.genchains([(nuc,)])
        maxrelerr = max([max(list(map(hl_relerr, c))) for c in chains])
        if maxrelerr > 0.1:
            continue
        t = int(t)
        if nuc not in nuc_keys:
            nuc_keys[nuc] = []
        nuc_keys[nuc].append(key)
    sfunc = lambda x: int(x.rsplit('_', 1)[1])
    for keys in nuc_keys.values():
        keys.sort(key=sfunc)
    return nuc_keys


def matdiff(x, y, threshold=1e-3):
    """Takes the difference between to materials, returns diff dict,
    the maximum relative error, the child which has the max difference,
    that child's mass, and the mass-weighted error.
    Skips nuclides that are not in both materials.
    """
    diff = {}
    maxrelerr = -1.0
    child = None
    childmass = None
    weightederr = None
    for nuc in x:
        if nuc not in y:
            continue
        xcomp = x[nuc]
        ycomp = y[nuc]
        if (xcomp <= threshold) and (ycomp <= threshold):
            continue
        diff[nuc] = d = np.abs(xcomp - ycomp) * (2 / (xcomp + ycomp))
        if d > maxrelerr:
            maxrelerr = d
            child = nuc
            childmass = max(xcomp, ycomp)
            weightederr = childmass * maxrelerr
    return diff, maxrelerr, child, childmass, weightederr


def mat_err_compare(nuc_keys, threshold=1e-3):
    """This returns a generator that compares the origen decayed
    material to the pyne material decay. The comparison stores:

    * the material key, matkey
    * the parent nuclide, nuc
    * the decay time, t
    * the maximum relative error over all nuclides, relerr
    * the nuclide with the maximum relative error, child
    * The mass of the child, childmass
    * the child's mass weighted relative error, weightederr

    Note that this also filters by when Origen itself has messed up by
    gaining or losing too much mass (1%). Furthermore, the maximum relative
    error is only computed for species that have a certain threshold
    unit mass (default 1e-3).
    """
    for nuc, keys in list(nuc_keys.items()):
        fresh = MATS[keys[0]]
        for matkey in keys[1:]:
            mat = MATS[matkey]
            if (mat.mass < 0.99) or (mat.mass > 1.1):
                continue  # Origen lost too much mass
            t = mat.metadata['decay_time']
            decmat = fresh.decay(t)
            row = matdiff(mat, decmat)
            if row[1] < 0.0 and row[2] is None:
                continue  # uncomparable materials, no common species
            row = (matkey, nuc, t) + row
            yield row

#
# Tests
#

def check_materr(row):
    maxrelerr = row[4]
    if maxrelerr < 0.01:
        return
    weightederr = row[7]
    assert_less(weightederr, 0.1)


def test_benchmark():
    nuc_keys = pivot_mat_keys()
    for row in mat_err_compare(nuc_keys):
        yield check_materr, row


if __name__ == "__main__":
    nose.runmodule()
