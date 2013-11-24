""" test data_checksums"""
from nose.tools import assert_equal


def test_data_checksums():
    from pyne.data import data_checksums
    assert_equal(len(data_checksums), 6)
    assert_equal(data_checksums['/neutron/simple_xs'], '3d6e086977783dcdf07e5c6b0c2416be')