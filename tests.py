#!/usr/bin/python

import argparse
import nose
from unittest import TestCase
from nose.tools import assert_equal, assert_in, assert_true, assert_almost_equal
from pyne import material
from pyne.material import Material


"""
function to load the material library
"""
def load_mat_lib(filename):
    mat_lib = material.Material()
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(
        filename, datapath="/material_library/materials", nucpath="/material_library/nucid")
    return mat_lib
    

"""
function to load the output h5m file created using the get_tag_values script
"""
def load_output(filename, material_library):
    output_lib= material.MaterialLibrary()
    output_lib.from_hdf5(filename)
    mat=material.Material()   
    for m in output_lib.items() :
        mat=m
        test_existence(material_library, mat)
        test_density(material_library, mat)
        test_composition(material_library, mat)
    return output_lib


"""
test to check the exsitence of materials from the output library in the PyNE material 
library
"""
def test_existence(material_library, material): 
    assert_in(material[0],material_library.iterkeys())
    
"""
test to check that the density of materials exsits and is in a proper format
"""    
               
def test_density(material_library, material):
    assert_true(material[1].density,float)
    
"""
test to check the composition
"""
def test_composition(material_library,material):
    for item in material_library.iteritems():
        if material[0] == item[0]:
            for c in material[1].comp.keys():
                assert_almost_equal(material[1].comp[c], item[1].comp[c], places=4)                 

"""
Parsing
"""
def parse():
    parser=argparse.ArgumentParser()
    parser.add_argument('-m',action='store',dest='model',help='the output file .h5m path')
    args=parser.parse_args()
    if not args.model:
        raise Exception('h5m path needed')
    return args
    
"""
main
"""    
def main():
    args=parse()
    # now load material library
    mat_lib = load_mat_lib("/home/moataz/.local/lib/python2.7/site-packages/pyne/nuc_data.h5")
    #load the output h5m file
    output_lib=load_output(args.model, mat_lib)
    
    

if __name__ == "__main__" :
    main()


