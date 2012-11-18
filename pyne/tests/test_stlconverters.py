"""Tests the part of stlconverters that is accessible from Python."""
###################
###  WARNING!!! ###
###################
# This file has been autogenerated

from unittest import TestCase
import nose

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

import os
import numpy  as np
import tables as tb

import pyne.stlconverters as conv
















# SetInt
def test_set_int():
    s = conv.SetInt()
    s.add(1)
    assert_true(1 in s)
    assert_true(-65 not in s)

    s = conv.SetInt([1, 42, -65])
    assert_true(42 in s)
    assert_true(18 not in s)



# SetStr
def test_set_str():
    s = conv.SetStr()
    s.add('Aha')
    assert_true('Aha' in s)
    assert_true('Me' not in s)

    s = conv.SetStr(['Aha', 'Take', 'Me'])
    assert_true('Take' in s)
    assert_true('On' not in s)



# MapStrStr
def test_map_str_str():
    m = conv.MapStrStr()
    m['Aha'] = 'On'
    m['Take'] = 'Me'
    assert_equal(len(m), 2)
    assert_equal(m['Take'], 'Me')

    m = conv.MapStrStr({'Me': 'Take', 'On': 'Aha'})
    assert_equal(len(m), 2)
    assert_equal(m['Me'], 'Take')

    n = conv.MapStrStr(m, False)
    assert_equal(len(n), 2)
    assert_equal(n['Me'], 'Take')

    # points to the same underlying map
    n['Take'] = 'Me'
    assert_equal(m['Take'], 'Me')



# MapStrInt
def test_map_str_int():
    m = conv.MapStrInt()
    m['Aha'] = 18
    m['Take'] = -65
    assert_equal(len(m), 2)
    assert_equal(m['Take'], -65)

    m = conv.MapStrInt({'Me': 42, 'On': 1})
    assert_equal(len(m), 2)
    assert_equal(m['Me'], 42)

    n = conv.MapStrInt(m, False)
    assert_equal(len(n), 2)
    assert_equal(n['Me'], 42)

    # points to the same underlying map
    n['Take'] = -65
    assert_equal(m['Take'], -65)



# MapIntStr
def test_map_int_str():
    m = conv.MapIntStr()
    m[1] = 'On'
    m[42] = 'Me'
    assert_equal(len(m), 2)
    assert_equal(m[42], 'Me')

    m = conv.MapIntStr({-65: 'Take', 18: 'Aha'})
    assert_equal(len(m), 2)
    assert_equal(m[-65], 'Take')

    n = conv.MapIntStr(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[-65], 'Take')

    # points to the same underlying map
    n[42] = 'Me'
    assert_equal(m[42], 'Me')



# MapStrUInt
def test_map_str_uint():
    m = conv.MapStrUInt()
    m['Aha'] = 18
    m['Take'] = 65
    assert_equal(len(m), 2)
    assert_equal(m['Take'], 65)

    m = conv.MapStrUInt({'Me': 42, 'On': 1})
    assert_equal(len(m), 2)
    assert_equal(m['Me'], 42)

    n = conv.MapStrUInt(m, False)
    assert_equal(len(n), 2)
    assert_equal(n['Me'], 42)

    # points to the same underlying map
    n['Take'] = 65
    assert_equal(m['Take'], 65)



# MapUIntStr
def test_map_uint_str():
    m = conv.MapUIntStr()
    m[1] = 'On'
    m[42] = 'Me'
    assert_equal(len(m), 2)
    assert_equal(m[42], 'Me')

    m = conv.MapUIntStr({65: 'Take', 18: 'Aha'})
    assert_equal(len(m), 2)
    assert_equal(m[65], 'Take')

    n = conv.MapUIntStr(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[65], 'Take')

    # points to the same underlying map
    n[42] = 'Me'
    assert_equal(m[42], 'Me')



# MapStrDouble
def test_map_str_dbl():
    m = conv.MapStrDouble()
    m['Aha'] = 18
    m['Take'] = -65.5555
    assert_equal(len(m), 2)
    assert_equal(m['Take'], -65.5555)

    m = conv.MapStrDouble({'Me': 42.42, 'On': 1.0})
    assert_equal(len(m), 2)
    assert_equal(m['Me'], 42.42)

    n = conv.MapStrDouble(m, False)
    assert_equal(len(n), 2)
    assert_equal(n['Me'], 42.42)

    # points to the same underlying map
    n['Take'] = -65.5555
    assert_equal(m['Take'], -65.5555)



# MapIntInt
def test_map_int_int():
    m = conv.MapIntInt()
    m[1] = 18
    m[42] = -65
    assert_equal(len(m), 2)
    assert_equal(m[42], -65)

    m = conv.MapIntInt({-65: 42, 18: 1})
    assert_equal(len(m), 2)
    assert_equal(m[-65], 42)

    n = conv.MapIntInt(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[-65], 42)

    # points to the same underlying map
    n[42] = -65
    assert_equal(m[42], -65)



# MapIntDouble
def test_map_int_dbl():
    m = conv.MapIntDouble()
    m[1] = 18
    m[42] = -65.5555
    assert_equal(len(m), 2)
    assert_equal(m[42], -65.5555)

    m = conv.MapIntDouble({-65: 42.42, 18: 1.0})
    assert_equal(len(m), 2)
    assert_equal(m[-65], 42.42)

    n = conv.MapIntDouble(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[-65], 42.42)

    # points to the same underlying map
    n[42] = -65.5555
    assert_equal(m[42], -65.5555)



# MapIntComplex
def test_map_int_complex():
    m = conv.MapIntComplex()
    m[1] = 0.18j
    m[42] = (-65.55-1j)
    assert_equal(len(m), 2)
    assert_equal(m[42], (-65.55-1j))

    m = conv.MapIntComplex({-65: (42+42j), 18: 1.0})
    assert_equal(len(m), 2)
    assert_equal(m[-65], (42+42j))

    n = conv.MapIntComplex(m, False)
    assert_equal(len(n), 2)
    assert_equal(n[-65], (42+42j))

    # points to the same underlying map
    n[42] = (-65.55-1j)
    assert_equal(m[42], (-65.55-1j))



