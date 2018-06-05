"""Tests Gamma Effects"""


import nose 
from nose.tools import assert_equal, assert_true, assert_almost_equal, assert_raises
from pyne import gamma_effects as ge

E = 1500.

def test_compton_edge():
    assert_equal(ge.compton_edge(E),1281.6861293078896)

def test_backscatter():
    assert_equal(ge.backscatter(E),218.3138706921105)

def test_single_escape():
    assert_equal(ge.single_escape(E),989.0)

def test_double_escape():
    assert_equal(ge.double_escape(E),478.0)
