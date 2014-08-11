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

# helper function to delete files before each test
def clean(files):
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

# writes the neutron tally to the file filename
def write_neutron(file_name):
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    tally.write_hdf5(file_name,"tally")

# writes the photon tally to the file filename
def write_photon(file_name):
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    tally.write_hdf5(file_name,"tally")

# writes the neutron tally to arb file and path
def write_arb_n(file_name,path):
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    tally.write_hdf5(file_name,path)

# writes the photon tally to arb file and path
def write_arb_p(file_name,path):
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    tally.write_hdf5(file_name,path)


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
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")

    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.tally_type,new_tally.tally_type)

def test_tally16():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")
   
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.particle_name,new_tally.particle_name)

def test_tally17():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")
  
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.entity_type,new_tally.entity_type)


def test_tally18():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")

    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.entity_name,new_tally.entity_name)


def test_tally19():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")

    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.entity_id,new_tally.entity_id)

def test_tally20():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")

    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_equal(tally.tally_name,new_tally.tally_name)

def test_tally21():
    clean(["test_tally.h5"])
    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    write_photon("test_tally.h5")

    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_equal(tally.entity_size,new_tally.entity_size)

################################################################################
# tests the append and load abilty
def test_tally22():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)
    assert_equal(tally.tally_type,new_tally.tally_type)

def test_tally23():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.particle_name,new_tally.particle_name)

def test_tally24():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.entity_type,new_tally.entity_type)


def test_tally25():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.entity_name,new_tally.entity_name)


def test_tally26():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.entity_id,new_tally.entity_id)

def test_tally27():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.tally_name,new_tally.tally_name)

def test_tally28():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    # there are now two tallies in the h5 file
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0) 
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",1)

    assert_equal(tally.entity_size,new_tally.entity_size)

################################################################################
# tests the load abilty
def test_tally29():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)
    assert_not_equal(tally.tally_type,new_tally.tally_type)

def test_tally30():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.particle_name,new_tally.particle_name)

def test_tally31():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.entity_type,new_tally.entity_type)


def test_tally32():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.entity_name,new_tally.entity_name)


def test_tally33():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.entity_id,new_tally.entity_id)

def test_tally34():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.tally_name,new_tally.tally_name)

def test_tally35():
    clean(["test_tally.h5"])
    write_photon("test_tally.h5")
    write_neutron("test_tally.h5")

    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","tally",0)

    assert_not_equal(tally.entity_size,new_tally.entity_size)


################################################################################
# test write and retrive from arbitrary location
def test_tally36():
    clean(["test_tally.h5"])
    write_arb_n("test_tally.h5","bob_geldof")
    tally = Tally("Current","neutron",14,"Surface","Surface 14","Neutron Current Across surface 14",100.0)
    
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","bob_geldof",0)
    assert_equal(tally.tally_type,new_tally.tally_type)

# test append and retrive from arbitrary location
def test_tally37():
    clean(["test_tally.h5"])
    write_arb_n("test_tally.h5","bob_geldof")
    write_arb_p("test_tally.h5","bob_geldof")

    tally = Tally("Flux","photon",12,"Volume","Volume 12","Photon Flux in Cell 12",35.0)
    
    new_tally = Tally()
    new_tally.from_hdf5("test_tally.h5","bob_geldof",1)
    assert_equal(tally.tally_type,new_tally.tally_type)



# Run as script
#
if __name__ == "__main__":
    nose.runmodule()
