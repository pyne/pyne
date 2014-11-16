import os
import io
import warnings
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose, \
    assert_array_almost_equal

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)

from pyne.endf2 import library

import nose
from nose.tools import assert_equal


def test_endcpp_basics():
    try:
        assert(os.path.isfile('U235.txt'))
    except AssertionError:
        try:
            import urllib.request as urllib
        except ImportError:
            import urllib
        urllib.urlretrieve("http://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/U/235",
                    "U235.txt")
    from hashlib import md5
    with open("U235.txt", "rb") as f:
        obs_hash = md5(f.read()).hexdigest()
    exp_hash = "1b71da3769d8b1e675c3c579ba5cb2d3"
    try:
        assert_equal(obs_hash, exp_hash)
    except AssertionError:
        raise AssertionError("U235.txt hash check failed; please try redownloading the U235 data file.")

    test_library = library()
    test_library.read_endf('U235.txt')
    contents = [[9228,    1,  451],
                [9228,    1,  452],
                [9228,    1,  455],
                [9228,    1,  456],
                [9228,    1,  458],
                [9228,    1,  460]]
    for index, item in enumerate(test_library.content_list):
        assert_array_equal(item, contents[index])
    allcontent = [[1, 451, 934, 7],
     [1, 452, 30, 6],
     [1, 455, 7, 7],
     [1, 456, 30, 7],
     [1, 458, 11, 4],
     [1, 460, 253745, 0],
     [2, 151, 3298, 7],
     [3, 1, 124, 7],
     [3, 2, 124, 7],
     [3, 3, 124, 7],
     [3, 4, 106, 7],
     [3, 16, 21, 6],
     [3, 17, 12, 6],
     [3, 18, 108, 7],
     [3, 19, 108, 7],
     [3, 20, 24, 6],
     [3, 21, 16, 6],
     [3, 37, 5, 6],
     [3, 38, 9, 6],
     [3, 51, 94, 6],
     [3, 52, 77, 6],
     [3, 53, 72, 6],
     [3, 54, 57, 6],
     [3, 55, 43, 6],
     [3, 56, 47, 6],
     [3, 57, 34, 7],
     [3, 58, 45, 7],
     [3, 59, 33, 7],
     [3, 60, 43, 7],
     [3, 61, 30, 7],
     [3, 62, 29, 7],
     [3, 63, 29, 7],
     [3, 64, 28, 7],
     [3, 65, 39, 7],
     [3, 66, 27, 7],
     [3, 67, 37, 7],
     [3, 68, 36, 7],
     [3, 69, 35, 7],
     [3, 70, 34, 7],
     [3, 71, 33, 7],
     [3, 72, 33, 7],
     [3, 73, 32, 7],
     [3, 74, 31, 7],
     [3, 75, 30, 7],
     [3, 76, 29, 7],
     [3, 77, 28, 7],
     [3, 78, 28, 7],
     [3, 79, 27, 7],
     [3, 80, 26, 7],
     [3, 81, 26, 7],
     [3, 82, 25, 7],
     [3, 83, 24, 7],
     [3, 84, 23, 7],
     [3, 85, 22, 7],
     [3, 86, 22, 7],
     [3, 87, 21, 7],
     [3, 88, 21, 7],
     [3, 89, 20, 7],
     [3, 90, 20, 7],
     [3, 91, 37, 7],
     [3, 102, 106, 6],
     [4, 2, 184, 6],
     [4, 18, 2, 6],
     [4, 51, 65, 6],
     [4, 52, 59, 6],
     [4, 53, 178, 6],
     [4, 54, 57, 6],
     [4, 55, 57, 6],
     [4, 56, 176, 6],
     [4, 57, 53, 7],
     [4, 58, 172, 7],
     [4, 59, 55, 7],
     [4, 60, 172, 7],
     [4, 61, 53, 7],
     [4, 62, 53, 7],
     [4, 63, 53, 7],
     [4, 64, 51, 7],
     [4, 65, 69, 7],
     [4, 66, 51, 7],
     [4, 67, 69, 7],
     [4, 68, 69, 7],
     [4, 69, 69, 7],
     [4, 70, 66, 7],
     [4, 71, 69, 7],
     [4, 72, 72, 7],
     [4, 73, 72, 7],
     [4, 74, 69, 7],
     [4, 75, 69, 7],
     [4, 76, 68, 7],
     [4, 77, 68, 7],
     [4, 78, 65, 7],
     [4, 79, 64, 7],
     [4, 80, 64, 7],
     [4, 81, 60, 7],
     [4, 82, 60, 7],
     [4, 83, 56, 7],
     [4, 84, 55, 7],
     [4, 85, 52, 7],
     [4, 86, 52, 7],
     [4, 87, 52, 7],
     [4, 88, 52, 7],
     [4, 89, 52, 7],
     [4, 90, 49, 7],
     [5, 18, 4346, 7],
     [5, 455, 589, 7],
     [6, 16, 699, 6],
     [6, 17, 210, 6],
     [6, 37, 42, 6],
     [6, 91, 1320, 6],
     [12, 4, 253, 6],
     [12, 18, 5, 6],
     [12, 102, 5, 6],
     [12, 460, 12928, 0],
     [13, 3, 8, 6],
     [14, 3, 1, 6],
     [14, 4, 1, 6],
     [14, 18, 1, 6],
     [14, 102, 1, 6],
     [14, 460, 1, 0],
     [15, 3, 127, 6],
     [15, 18, 54, 6],
     [15, 102, 58, 6],
     [31, 452, 56, 7],
     [31, 456, 199, 7],
     [33, 1, 158903, 7],
     [33, 2, 113418, 7],
     [33, 4, 319, 7],
     [33, 16, 86, 7],
     [33, 17, 29, 7],
     [33, 18, 81910, 7],
     [33, 102, 23060, 7],
     [35, 18, 6606, 7]]
    for index, item in enumerate(test_library.get_mt451(9228, 1, 451).mt_list):
        assert_array_equal(allcontent[index], item)


if __name__ == "__main__":
    nose.runmodule()
