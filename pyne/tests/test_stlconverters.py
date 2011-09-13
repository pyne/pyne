"""Tests the part of stlconverters that is accessible from Python."""

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

import os
import numpy  as np
import tables as tb

import pyne.stlconverters as conv


#
# Test set proxies
# 

def test_set_int():
    s = conv.SetInt()
    s.add(7)
    assert_true(7 in s)
    assert_true(11 not in s)

    s = conv.SetInt([1, 43, 16])
    assert_true(1 in s)
    assert_true(13 not in s)

    s = conv.SetInt({i**2 for i in range(10)})
    assert_true(4 in s)
    assert_true(42 not in s)



def test_set_str():
    s = conv.SetStr()
    s.add("old")
    assert_true("old" in s)
    assert_true("new" not in s)

    s = conv.SetStr(["Aha", "Take", "Me", "On"])
    assert_true("Aha" in s)
    assert_true("Captain Hammer" not in s)

    s = conv.SetStr("Dr. Horrible")
    assert_true('r' in s)
    assert_true('Sweetness' not in s)



#
# Test map proxies
# 


def test_map_str_int():
    m = conv.MapStrInt()
    m['worm'] = 69
    m['drum'] = 100
    assert_equal(len(m), 2)
    assert_equal(m['worm'], 69)

    m = conv.MapStrInt({'yes': 1, 'no': 0})
    assert_equal(len(m), 2)
    assert_equal(m['no'], 0)

    n = conv.MapStrInt(m, False)
    assert_equal(len(m), 2)
    assert_equal(m['yes'], 1)

    # points to the same underlying map
    n['maybe'] = -1 
    assert_equal(m['maybe'], -1)


def test_map_int_str():
    m = conv.MapIntStr()
    m[88] = 'two eight'
    m[100] = 'cent'
    assert_equal(len(m), 2)
    assert_equal(m[100], 'cent')

    m = conv.MapIntStr({1: 'yes', 0: 'no'})
    assert_equal(len(m), 2)
    assert_equal(m[0], 'no')

    n = conv.MapIntStr(m, False)
    assert_equal(len(m), 2)
    assert_equal(m[1], 'yes')

    # points to the same underlying map
    n[-1] = 'maybe' 
    assert_equal(m[-1], 'maybe')



def test_map_int_double():
    m = conv.MapIntDouble()
    m[88] = 88.88
    m[100] = 0.01
    assert_equal(len(m), 2)
    assert_equal(m[100], 0.01)

    m = conv.MapIntDouble({1: 4.2, 0: 3.14})
    assert_equal(len(m), 2)
    assert_equal(m[0], 3.14)

    n = conv.MapIntDouble(m, False)
    assert_equal(len(m), 2)
    assert_equal(m[1], 4.2)

    # points to the same underlying map
    n[-1] = 10.1
    assert_equal(m[-1], 10.1)



def test_map_int_complex():
    m = conv.MapIntComplex()
    m[88] = 88.88 + 0j
    m[100] = 0.01 - 66.5j
    assert_equal(len(m), 2)
    assert_equal(m[100], 0.01 - 66.5j)

    m = conv.MapIntComplex({1: 4.2 + 42.0j, 0: 3.14 - 2.2j})
    assert_equal(len(m), 2)
    assert_equal(m[0], 3.14 - 2.2j)

    n = conv.MapIntComplex(m, False)
    assert_equal(len(m), 2)
    assert_equal(m[1], 4.2 + 42.0j)

    # points to the same underlying map
    n[-1] = 10.1 + 3.3j
    assert_equal(m[-1], 10.1 + 3.3j)




#
# Run as script
#
if __name__ == "__main__":
    nose.main()
