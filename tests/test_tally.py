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

def test_tally6():
    tally = Tally()
    assert_equal(tally.tally_name,"")

def test_tally7():
    tally = Tally()
    assert_equal(tally.entity_size,-1.0)


###############################################################################
# tests the constructor
def test_tally8():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.tally_type,"Flux")

def test_tally9():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.particle_name,"photon")

def test_tally10():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.entity_type,"Volume")

def test_tally11():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.entity_name,"Volume 12")

def test_tally12():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.entity_id,12)

def test_tally13():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.tally_name,"Photon Flux in Cell 12")

def test_tally14():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    assert_equal(tally.entity_size,35.0)

################################################################################
# tests the save abilty
def test_tally15():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    tally.write_hdf5("test_tally.h5","tally")
    
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.tally_type,new_tally.tally_type)

def test_tally16():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
   
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_equal(tally.particle_name,new_tally.particle_name)

def test_tally17():
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
   
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_equal(tally.entity_type,new_tally.entity_type)



# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
