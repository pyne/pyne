import os
import io
import warnings
import sys
from hashlib import md5

import numpy as np
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
#    exp_nuclides = [10010000, 30060000, 40090000, 50110000, 30070000,
#                    60000000, 50100000]
#    obs_nuclides = list(map(int, testlib.structure.keys()))
#    assert_array_equal(exp_nuclides, obs_nuclides)

if __name__ == "__main__":
    nose.runmodule()
