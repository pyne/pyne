#!/usr/bin/python 

from itaps import iMesh,iBase
from pyne import material
from pyne.material import Material
import string
import argparse

"""
gets all tags on dagmc geometry
------------------------------
filename : the dagmc filename
"""
def get_tag_values(filename):
    mesh = iMesh.Mesh()
    mesh.load(filename)

   # ents = mesh.getEntities()
   # print mesh.getEntType(ents)
    ents = mesh.getEntities(iBase.Type.all, iMesh.Topology.triangle)
    #print len(ents)
    mesh_set = mesh.getEntSets()
    #print len(mesh_set)

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
    print('The matrials group names found in the h5m file are: ')
    print tag_values


"""
function which loads pyne material library
"""
#nuc_data='/home/moataz/.local/lib/python2.7/site-packages/pyne/nuc_data.h5'
def load_mat_lib(filename):
    global mat_lib
    mat_lib=material.Material()
    mat_lib=material.MaterialLibrary()
    mat_lib.from_hdf5(filename,datapath="/material_library/materials",nucpath="/material_library/nucid")
   # print mat_lib.keys()
    

"""
function to check that materials exist in library
-------------------------------------------------
tag_values - list of tags from dagmc file
mat_lib - pyne material library instance
"""
def check_matname(tag_values,mat_lib):
    global mat_list, d
    global materials_list
    # loop over tags 
    mat_list=[]   
    d=0
    for tag in tag_values :
        # material=steel
        #name=tag.split("=")
        # find material tag
        mat_name=[]
        mat_list_matname=[]
        if "mat" in tag:
            if "/" in tag :
                mat_name = tag.split("/")
                # list of material names from tagged geometry
                mat_list_matname.append(mat_name[0])
            else :
                mat_list_matname.append(tag) 
            for matname in mat_list_matname :
		  try: 
                       mat_name=matname.split(':')
                       mat_list.append(mat_name[1])         
                  except:
                       print("Could not find group name in approaprite format"), tag
	             
    if len(mat_list) == 0:
	print("no group names found")
	exit()                      
    print mat_list
    # for the sake of testing, fmat_list >>> fluka materials
    global fmat_list
    fmat_list=[]
    # list of pyne materials to add to h5m
    materials_list=[]
    # if name from geom matches name in lib
    # loop over materials in geometry
    for u in range(len(mat_list)):
       item=mat_list[u]
       # loop over materials in library
       for key in mat_lib.iterkeys():  
           if item == key :
                d=d+1
                # get the material
                new_mat = mat_lib.get(key)
                # set the mcnp material number
                set_attrs(new_mat,d, code)
                materials_list.append(new_mat)
                break
           if mat_lib.keys().index(key) == len(mat_lib.keys())-1:	
                print('material {%s} doesn''t exist in pyne material lib' %item)
                print_near_match(item,mat_lib)
                exit()
    print fmat_list
    
    # check that there are as many materials as there are groups
    if d != len(mat_list):
	print "There are insuficient materials"
	exit()
    # list of pyne material objects
    #print materials_list/
    print materials_list
    return materials_list
    
""" 
function to print near matches to material name
"""
def print_near_match(material,material_library):
    #p=open('w.txt','w')
    #m=Material()
    for item in material_library.iterkeys() :  
        if ( material.lower() in item.lower()) or (material.upper() in item.upper()) :
	    print "near matches to ", material, " are " 
	    print item
            return
           # p.write(str(item))
           # p.write('\n')
           # print material_library.get(item)
           # p.write(str(material_library.get(item)))
           # p.write('\n')
           # p.write('\n')    
    #p.close()

"""
function to set the attributes of the materials:
"""
def set_attrs(mat,number,code):
    if code is 'mcnp' or 'both' :     
        mat.attrs['mat_number']=str(number)
    if code == 'fluka' or 'both' :   
        fluka_material_naming(mat,number)
    return
     

"""
Function to prepare fluka material names:
"""
def fluka_material_naming(material,number) :
    matf=material.attrs['name']
    matf=''.join(c for c in matf if c.isalnum())
    if len(matf) <= 8 :
        if matf.upper() in fmat_list :
            if number <= 9 :
                matf=matf.rstrip(matf[-1])
                matf=matf+str(number)
            if number >= 9 and number <= 99 :
                for i in range(2) :
                    matf=matf.rstrip(matf[-1])
                matf=matf+str(number)
            fmat_list.append(matf.upper())    
        else :            
            fmat_list.append(matf.upper())
    else :
        matf=matf[0:8]
        if matf.upper() in fmat_list :
            if number <= 9 :
                matf=matf.rstrip(matf[-1])
                matf=matf+str(number)
            if number >= 9 and number <= 99 :
                for i in 2 :
                    matf=matf.rstrip(matf[-1])
                matf=matf+str(number)
            fmat_list.append(matf.upper())    
        else :            
            fmat_list.append(matf.upper())
    material.attrs['name']=matf.upper()   
    material.attrs['mat_number']=str(number)
    return fmat_list
    return material


"""
Function write_mats, writes material objects to hdf5 file
-------
material_list: list of PyNE Material Objects
filename: filename to write the objects to
"""
def write_mats_h5m(material_list,filename):
    for material in material_list:
	material.write_hdf5(filename)


"""
function to parse the script, adding options:
defining 
-f  : the .h5m file path
-d  : nuc_data path
"""
def parsing(parsescript) :
    parser=argparse.ArgumentParser()
    parser.add_argument('-f', action='store',dest='datafile', help='the path to the .h5m file')
    parser.add_argument('-d', action='store',dest='nuc_data', help='the path to the PyNE materials library   nuc_data.h5')
    parser.add_argument('-c', action='store',dest='code',help='the format of the output h5m file; mcnp or fluka')
    parser.add_argument('-o',action='store',dest='output',help='the name of the output file')
    args=parser.parse_args()
    if args.datafile:
       #print args
       global datafile
       datafile=args.datafile
    else :
       print('h5m file not found!!. [-f] is not set')   
       exit()   
    if args.nuc_data :
       global nuc_data
       nuc_data=args.nuc_data
    else :
       print('nuc_data file not found!!. [-d] is not set') 
       exit()
    if args.code :
       global code
       code=args.code
    else :
       print('output file format is not specified!!. [-c] is not set') 
       exit()
    if args.output :
       global output
       output=args.output
    else :
       global output
       output='output.h5m'


#parse the script
parsing(1)            
# get list of tag valuesaz
get_tag_values(datafile)
# now load material library
load_mat_lib(nuc_data)
# check that material tags exist in library
# material_list is list of pyne objects in problem
material_list=check_matname(tag_values,mat_lib)
# write materials to file
write_mats_h5m(material_list,output)

exit()
