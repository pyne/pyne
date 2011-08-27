"""PyNE nuclear data tests"""
from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, assert_in

from pyne import data

def test_nuc_weight():
    o16 = [15.99491461957, 16.0]
    u235 = [235.043931368, 235.0]
    am242m = [242.059550625, 242.0]

    # zzaam form
    assert_in(data.nuc_weight(80160),  o16)
    assert_in(data.nuc_weight(922350), u235)
    assert_in(data.nuc_weight(952421), am242m)


if __name__ == "__main__":
    nose.main()

