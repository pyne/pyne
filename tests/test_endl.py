from __future__ import print_function
import os
import io
import warnings
import sys
from hashlib import md5

import numpy as np
from numpy.testing import assert_array_equal
import nose
from nose.tools import assert_equal

from pyne.utils import QAWarning
warnings.simplefilter("ignore", QAWarning)
from pyne.endl import Library

if sys.version_info[0] > 2:
    from io import StringIO
else:
    from StringIO import StringIO

if sys.version_info[0] > 2:
    import urllib.request as urllib
else:
    import urllib


def test_loadfile():
    if not os.path.isfile('epdl97_eedl_Pb'):
        urllib.urlretrieve(
                "https://www-nds.iaea.org/epdl97/data/endl/eedl/za082000",
                "epdl97_eedl_Pb"
                )
    with open("epdl97_eedl_Pb", "rb") as f:
        obs_hash = md5(f.read()).hexdigest()
    exp_hash = "502105669e0c0ad917301d219c649aaf"
    assert_equal(
            obs_hash, exp_hash,
            msg="epdl97_eedl_Pb hash check failed; please try redownloading"
            " the epdl97_eedl_Pb data file."
            )
    testlib = Library("epdl97_eedl_Pb")

    # test the nuclides
    pb_nuclide = 820000000
    exp_nuclides = [pb_nuclide]
    obs_nuclides = list(map(int, testlib.structure.keys()))
    assert_array_equal(exp_nuclides, obs_nuclides)

    # test the incoming particles
    exp_pin = [9]
    obs_pin = testlib.structure[pb_nuclide]['pin']
    assert_array_equal(sorted(exp_pin), sorted(obs_pin))

    # test the reaction properties
    exp_rdesc = [7, 8, 10, 81, 82, 83]
    obs_rdesc = testlib.structure[pb_nuclide]['rdesc']
    assert_array_equal(sorted(exp_rdesc), sorted(obs_rdesc))

    # test the reaction descriptors
    exp_rprop = [0, 10, 11, 21, 22]
    obs_rprop = testlib.structure[pb_nuclide]['rprop']
    assert_array_equal(sorted(exp_rprop), sorted(obs_rprop))

if __name__ == "__main__":
    nose.runmodule()
