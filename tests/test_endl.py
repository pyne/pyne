from __future__ import print_function
import os
import io
import warnings
import sys
from hashlib import md5

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import nose
from nose.tools import assert_equal

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
from pyne.endl import Library

if sys.version_info[0] > 2:
    from io import StringIO
else:
    from StringIO import StringIO


from utils import download_file


def test_loadfile():
    download_file(
        "https://www-nds.iaea.org/epdl97/data/endl/eedl/za082000",
        "epdl97_eedl_Pb",
        "502105669e0c0ad917301d219c649aaf",
    )
    testlib = Library("epdl97_eedl_Pb")

    # test the nuclides
    pb_nuclide = 820000000
    exp_nuclides = [pb_nuclide]
    obs_nuclides = list(map(int, testlib.structure.keys()))
    assert_array_equal(exp_nuclides, obs_nuclides)

    # test the incoming particles
    exp_pin = [9]
    obs_pin = testlib.structure[pb_nuclide]["pin"]
    assert_array_equal(sorted(exp_pin), sorted(obs_pin))

    # test the reaction properties
    exp_rdesc = [7, 8, 10, 81, 82, 83]
    obs_rdesc = testlib.structure[pb_nuclide]["rdesc"]
    assert_array_equal(sorted(exp_rdesc), sorted(obs_rdesc))

    # test the reaction descriptors
    exp_rprop = [0, 10, 11, 21, 22]
    obs_rprop = testlib.structure[pb_nuclide]["rprop"]
    assert_array_equal(sorted(exp_rprop), sorted(obs_rprop))


def test_elastic_cross_section():
    testlib = Library("files_test_endl/testfile")
    # elastic-scattering cross section (2 columns)
    pb_nuclide = "Pb"
    p_in = 9  # electron
    rdesc = 10  # elastic scattering
    rprop = 0  # integrated cross section
    data = testlib.get_rx(pb_nuclide, p_in, rdesc, rprop)
    # test the min and max values
    min_values = [1e-5, 3.63530e5]
    max_values = [1e5, 8.42443e9]
    assert_allclose(np.min(data, axis=0), min_values)
    assert_allclose(np.max(data, axis=0), max_values)


def test_bremsstrahlung_spectra():
    testlib = Library("files_test_endl/testfile")
    # bremsstrahlung spectra (3 columns)
    pb_nuclide = "Pb"
    p_in = 9  # electron
    rdesc = 82  # bremsstrahlung
    rprop = 21  # spectra of secondary photon
    data = testlib.get_rx(pb_nuclide, p_in, rdesc, rprop)
    # test the min and max values
    min_values = [1e-5, 1e-7, 4.61320e-8]
    max_values = [1e5, 1e5, 9.15055e6]
    assert_allclose(np.min(data, axis=0), min_values)
    assert_allclose(np.max(data, axis=0), max_values)


def test_ionization_cross_section():
    testlib = Library("files_test_endl/testfile")
    # ionization cross section (2 columns)
    # this is non-trivial because this table is indexed by subshell
    pb_nuclide = "Pb"
    p_in = 9  # electron
    rdesc = 81  # ionization
    rprop = 0  # integrated cross section
    x1 = 1  # subshell
    data = testlib.get_rx(pb_nuclide, p_in, rdesc, rprop, x1=x1)
    # test the min and max values
    min_values = [8.82900e-2, 5.66158e-1]
    max_values = [1e5, 7.52177e1]
    assert_allclose(np.min(data, axis=0), min_values)
    assert_allclose(np.max(data, axis=0), max_values)


def test_ionization_spectra():
    testlib = Library("files_test_endl/testfile")
    # ionization spectra of recoil electron (3 columns)
    # this is non-trivial because this table is indexed by subshell and
    # outgoing particle
    pb_nuclide = "Pb"
    p_in = 9  # electron
    rdesc = 81  # ionization
    rprop = 21  # spectra of secondary electron
    x1 = 1  # subshell
    p_out = 19  # electron as recoil
    data = testlib.get_rx(pb_nuclide, p_in, rdesc, rprop, x1=x1, p_out=p_out)
    # test the min and max values
    min_values = [8.82900e-2, 1e-8, 1.33833e-11]
    max_values = [1e5, 5e4, 7.06473e7]
    assert_allclose(np.min(data, axis=0), min_values)
    assert_allclose(np.max(data, axis=0), max_values)


if __name__ == "__main__":
    nose.runmodule()
