""" test data_checksums and hashing functions"""
import os
import warnings

from nose.tools import assert_equal, assert_true

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
import pyne

# These tests require nuc_data
if not os.path.isfile(pyne.nuc_data):
    raise RuntimeError("Tests require nuc_data.h5.  Please run nuc_data_make.")

def test_data_checksums():
    from pyne.data import data_checksums
    assert_equal(len(data_checksums), 5)
    assert_equal(data_checksums['/neutron/simple_xs'],
                 '3d6e086977783dcdf07e5c6b0c2416be')

def test_internal_hashes():
    from pyne.dbgen import hashtools
    hashtools.set_internal_hashes(pyne.nuc_data)
    for item, val in hashtools.check_internal_hashes(pyne.nuc_data):
        assert_true(val)

