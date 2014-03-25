#!/usr/bin/python

from  pyne import material #, MaterialLibrary
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
        compare(material_library, mat)
    return output_lib


"""
function to compare items from the output library with that in the material 
library
"""
def compare(material_library, material): 
        counter=0
        for key in material_library.items():
            counter=counter+1
            if material == key :
                print "OK", material
                print("+++++++++++++++++")
                print key
                break
            elif material[0] == key[0] and material[1].attrs["density"] != key[1].attrs["density"] :
                print "difference in density", material
                #raise Exception('error found') 
                print "///////////////"   
            elif counter == len(material_library.items())-1 :
                print "no match found for %s" %material[0]
                print material

def main():
    # now load material library
    mat_lib = load_mat_lib("/home/moataz/.local/lib/python2.7/site-packages/pyne/nuc_data.h5")
    #load the output h5m file
    output_lib=load_output("output.h5m", mat_lib)
    
    

if __name__ == "__main__" :
    main()


