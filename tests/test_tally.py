"""Tally tests"""
import os
import warnings

from unittest import TestCase
import nose 

from nose.tools import assert_equal, assert_not_equal, assert_raises, raises, \
    assert_almost_equal, assert_true, assert_false, assert_in

from pyne.utils import VnVWarning
warnings.simplefilter("ignore", VnVWarning)
from pyne.tally import Tally
from pyne import jsoncpp 
from pyne import data
import numpy  as np
import tables as tb


###############################################################################
# tests the default constructor
def test_tally1():
    tally = Tally()
    assert_equal(tally.tally_type,"")

def test_tally2():
    tally = Tally()
    assert_equal(tally.particle_name,"")

def test_tally3():
    tally = Tally()
    assert_equal(tally.entity_type,"")

def test_tally4():
    tally = Tally()
    assert_equal(tally.entity_name,"")

def test_tally5():
    tally = Tally()
    assert_equal(tally.entity_id,-1)

###############################################################################
# tests the constructor
def test_tally6():
    tally = Tally("flux","photon",12,"volume","volume 12")
    assert_equal(tally.tally_type,"flux")

def test_tally7():
    tally = Tally("flux","photon",12,"volume","volume 12")
    assert_equal(tally.particle_name,"photon")

def test_tally8():
    tally = Tally("flux","photon",12,"volume","volume 12")
    assert_equal(tally.entity_type,"volume")

def test_tally9():
    tally = Tally("flux","photon",12,"volume","volume 12")
    assert_equal(tally.entity_name,"volume 12")

def test_tally10():
    tally = Tally("flux","photon",12,"volume","volume 12")
    assert_equal(tally.entity_id,12)

# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
