#!/usr/bin/python 


datafile='test_models/test_4.h5m'

"""
gets all tags on dagmc geometry
------------------------------
filename : the dagmc filename
"""
from itaps import iMesh,iBase
import string

def get_tag_values(filename):
    mesh = iMesh.Mesh()
    mesh.load(filename)

    ents = mesh.getEntities()
    print mesh.getEntType(ents)
    ents = mesh.getEntities(iBase.Type.all, iMesh.Topology.triangle)
    print len(ents)
    mesh_set = mesh.getEntSets()
    print len(mesh_set)

    global tag_values
    tag_values=[]
    found_all_tags=0
    for i in mesh_set:
        if found_all_tags == 1:
            break
        
    # get all the tags
        tags = mesh.getAllTags(i)
    #print tags
    # loop over the tags checking for name
        for t in tags:
            # if we find name
            if t.name == 'NAME':
                # the handle is the tag name
                t_handle = mesh.getTagHandle(t.name)
	    # get the data for the tag, with taghandle in element i
                tag = t_handle[i]
                a=[]
	    # since we have a byte type tag loop over the 32 elements
                for part in tag:
                    # if the byte char code is non 0
                #print part
                    if (part != 0 ):
                        # convert to ascii 
                        a.append(str(unichr(part)))
		    # join to end string
                        test=''.join(a)
                        # the the string we are testing for is not in the list of found
                        # tag values, add to the list of tag_values
                if not any(test in s for s in tag_values):
                    # print test
                    tag_values.append(test)
		    #print test
		    # if the tag is called impl_complement, this is the 
		    # last tag we are done
                if any('impl_complement' in s for s in tag_values):
                    found_all_tags=1 
    print tag_values

""" 
function to print near matches to material name
"""
def print_near_match(material,material_library):
    for item in material_library.iterkeys():  
#        print item
        if ( material.lower() in item.lower()) or (material.upper() in item.upper()) :
	    print "near matches to ", material, " are " 
	    print item
            print material_library.get(item)

"""
function which loads pyne material library
"""
nuc_data='/home/moataz/.local/lib/python2.7/site-packages/pyne/nuc_data.h5'
from pyne import material
from pyne.material import Material
def load_mat_lib(filename):
    mat_lib=material.Material()
    mat_lib=material.MaterialLibrary()
    mat_lib.from_hdf5(filename,datapath="/material_library/materials",nucpath="/material_library/nucid")
    #print mat_lib
    global mat_lib
"""
function to check that materials exist in library
-------------------------------------------------
tag_values - list of tags from dagmc file
mat_lib - pyne material library instance
"""

def check_matname(tag_values,mat_lib):
    # loop over tags 
    mat_list=[]
    global mat_list
    for tag in tag_values :
        # material=steel
        #name=tag.split("=")
        # find material tag
        mat_name=[]
        mat_list_matname=[]
        if "mat" in tag:
            mat_name = tag.split("/")
            # list of material names from tagged geometry
            mat_list_matname.append(mat_name[0]) 
            for matname in mat_list_matname :
                  mat_name=matname.split(':')
                  mat_list.append(mat_name[1])         

    print mat_list
    # list of pyne materials to add to h5m
    materials=[]
    for key in mat_lib.iterkeys() :    
    # if name from geom matches name in lib
        # loop over materials in geometry
        for item in mat_list:
        # loop over materials in library
            if item in key :
                # get the material
                new_mat = mat_lib.get(key)
                materials.append(new_mat)

            else :
                print('material {%s} doesn''t exist in pyne material lib' %item)
                print_near_match(item,mat_lib)
                exit()
                continue 
    # list of pyne material objects
    print materials
    return materials

            
# get list of tag valuesaz
get_tag_values(datafile)
# now load material library
load_mat_lib(nuc_data)
# check that material tags exist in library
material_objects=check_matname(tag_values,mat_lib)
exit() 
