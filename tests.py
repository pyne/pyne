#!/usr/bin/python

import argparse
from  pyne import material
from pyne.material import Material

"""
function to load the material library
"""
def load_mat_lib(filename):
    mat_lib = material.Material()
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(
        filename, datapath="/material_library/materials", nucpath="/material_library/nucid")
    for item in mat_lib.iterkeys():
        print item
    exit()
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
        compare(material_library, mat)
    return output_lib


"""
function to compare items from the output library with that in the material 
library
"""
def compare(material_library, material): 
        counter=0
        for item in material_library.items():
            counter=counter+1
            if material[0] == item[0] :
                print 'material :' , material
                print 'item :' ,item
                if material[1].density == item[1].density:
                    print "complete match found for : ", material[0]
                elif material[1].density != item[1].density :
                    print "match found with a difference in density for : ", material[0]   
                elif counter == len(material_library.items())-1 :
                    print "no match found for %s" %material[0]

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


