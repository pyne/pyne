"""Transmuter tests"""
from nose.tools import assert_equal, assert_almost_equal

from pyne import data
from pyne import cram
from pyne import nucname
from pyne import transmuters


def test_transmuters_cram():
    n0 = {"H3": 1.0}
    A = -cram.DECAY_MATRIX * data.half_life("H3")
    n1 = transmuters.cram(A, n0, order=16)
    assert_equal(2, len(n1))
    assert_almost_equal(0.5, n1[nucname.id("H3")])
    assert_almost_equal(0.5, n1[nucname.id("He3")])


# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
