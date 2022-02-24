""" test data_checksums and hashing functions"""
import os
import warnings
from shutil import copyfile

import nose
from nose.tools import assert_equal, assert_true

from pyne.utils import QAWarning

warnings.simplefilter("ignore", QAWarning)
import pyne

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")


def test_data_checksums():
    from pyne.data import data_checksums

    assert_equal(len(data_checksums), 7)
    assert_equal(data_checksums["/atomic_mass"], "10edfdc662e35bdfab91beb89285efff")
    assert_equal(
        data_checksums["/material_library"], "8b10864378fbd88538434679acf908cc"
    )
    assert_equal(data_checksums["/neutron/eaf_xs"], "29622c636c4a3a46802207b934f9516c")
    assert_equal(
        data_checksums["/neutron/scattering_lengths"],
        "a24d391cc9dc0fc146392740bb97ead4",
    )
    assert_equal(
        data_checksums["/neutron/simple_xs"], "3d6e086977783dcdf07e5c6b0c2416be"
    )
    assert_equal(data_checksums["/decay"], "4f41f3e46f4306cc44449f08a20922e0")
    assert_equal(data_checksums["/dose_factors"], "dafa32c24b2303850a0bebdf3e6b122e")


def test_internal_hashes():
    from pyne.dbgen import hashtools

    remove_data = False
    if os.access(pyne.nuc_data, os.W_OK):
        test_data = pyne.nuc_data
    else:
        # Create a copy so we don't try to modify a file we don't have
        # permissions for
        test_data = "test_nuc_data.h5"
        copyfile(pyne.nuc_data, test_data)
        remove_data = True
    hashtools.set_internal_hashes(test_data)
    for item, val in hashtools.check_internal_hashes(test_data):
        assert_true(val)
    if remove_data:
        os.remove(test_data)


if __name__ == "__main__":
    nose.runmodule()
