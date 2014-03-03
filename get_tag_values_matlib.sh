<<<<<<< HEAD:get_tag_values_matlib.sh
#!/usr/bin/python
=======
#!/usr/env/python 
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py

from itaps import iMesh, iBase
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
<<<<<<< HEAD:get_tag_values_matlib.sh

=======
    
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    # get all entities
    ents = mesh.getEntities(iBase.Type.all, iMesh.Topology.triangle)
    # get mesh set
    mesh_set = mesh.getEntSets()

<<<<<<< HEAD:get_tag_values_matlib.sh
    tag_values = []  # list of tag values
    found_all_tags = 0
=======
    tag_values=[] # list of tag values
    found_all_tags=0
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    for i in mesh_set:
        if found_all_tags == 1:
            break

    # get all the tags
        tags = mesh.getAllTags(i)
        # loop over the tags checking for name
        for t in tags:
            # look for NAME tag
            if t.name == 'NAME':
                # the handle is the tag name
                t_handle = mesh.getTagHandle(t.name)
                # get the data for the tag, with taghandle in element i
                tag = t_handle[i]
<<<<<<< HEAD:get_tag_values_matlib.sh
                a = []
                # since we have a byte type tag loop over the 32 elements
                for part in tag:
                    # if the byte char code is non 0
                    if (part != 0):
                        # convert to ascii
                        a.append(str(unichr(part)))
                        # join to end string
                        test = ''.join(a)
=======
                a=[]
                # since we have a byte type tag loop over the 32 elements
                for part in tag:
                    # if the byte char code is non 0
                    if (part != 0 ):
                        # convert to ascii 
                        a.append(str(unichr(part)))
                        # join to end string
                        test=''.join(a)
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
                        # the the string we are testing for is not in the list of found
                        # tag values, add to the list of tag_values
                # if not already in list append to lilst
                if not any(test in s for s in tag_values):
                    tag_values.append(test)
                # last tag we are done
                if any('impl_complement' in s for s in tag_values):
<<<<<<< HEAD:get_tag_values_matlib.sh
                    found_all_tags = 1
=======
                    found_all_tags=1 
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py

    print('The matrials group names found in the h5m file are: ')
    print tag_values
    return tag_values

"""
function which loads pyne material library
"""
<<<<<<< HEAD:get_tag_values_matlib.sh


def load_mat_lib(filename):
    mat_lib = material.Material()
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(
        filename, datapath="/material_library/materials", nucpath="/material_library/nucid")
    return mat_lib

=======
def load_mat_lib(filename):
    mat_lib=material.Material()
    mat_lib=material.MaterialLibrary()
    mat_lib.from_hdf5(filename,datapath="/material_library/materials",nucpath="/material_library/nucid")
    return  mat_lib
    
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
"""
function to check that materials exist in library
-------------------------------------------------
tag_values - list of tags from dagmc file
mat_lib - pyne material library instance
"""
<<<<<<< HEAD:get_tag_values_matlib.sh


def check_matname(tag_values):
    # loop over tags
    mat_list = []   # list of materials
    d = 0  # counter of the

    mat_list_matname = []  # list of material names

=======
def check_matname(tag_values):
    # loop over tags 
    mat_list = []   # list of materials
    d = 0  # counter of the 

    mat_list_matname = [] # list of material names 
        
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    # loop over the tags in the file
    for tag in tag_values:
        # look for mat, this id's material in group name
        if "mat" in tag:
            # split on the basis of "/" being delimiter
<<<<<<< HEAD:get_tag_values_matlib.sh
            if "/" in tag:
                mat_name = tag.split("/")
                # list of material name only
                mat_list_matname.append(mat_name[0])
            # otherwise we have only mat:
            else:
                mat_list_matname.append(tag)

    # split colons from name
    for matname in mat_list_matname:
        try:
            mat_name = matname.split(':')
            mat_list.append(mat_name[1])
        except:
            print("Could not find group name in appropriate format"), tag
            exit
    # error conditions, not tags found
    if len(mat_list) == 0:
        print("no group names found")
        exit()
=======
            if "/" in tag: 
                mat_name = tag.split("/")
                # list of material name only 
                mat_list_matname.append(mat_name[0])
            # otherwise we have only mat:
            else :
                mat_list_matname.append(tag) 

    # split colons from name
    for matname in mat_list_matname:
        try: 
            mat_name=matname.split(':')
            mat_list.append(mat_name[1])         
        except:
            print("Could not find group name in appropriate format"), tag	             
            exit
    # error conditions, not tags found
    if len(mat_list) == 0:
	print("no group names found")
	exit()                      
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py

    print mat_list
    return mat_list

<<<<<<< HEAD:get_tag_values_matlib.sh

def check_and_create_materials(material_list, mat_lib):
=======
def check_and_create_materials(material_list,mat_lib):
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    # for the sake of testing, fmat_list >>> fluka materials
    flukamaterials_list = []
    material_object_list = []
    d = 0
    # loop over materials in geometry
    for material in material_list:
<<<<<<< HEAD:get_tag_values_matlib.sh
        # loop over materials in library
        for key in mat_lib.iterkeys():
            if material == key:
                d = d + 1
                # get the material
                new_mat = mat_lib.get(key)[:]
                # copy attrs 'cos python is dumb
                copy_attrs(new_mat, mat_lib.get(key))

                # set the mcnp material number or fluka material name
                set_attrs(new_mat, d, code, flukamaterials_list)
                material_object_list.append(new_mat)
                break
            if mat_lib.keys().index(key) == len(mat_lib.keys()) - 1:
                print(
                    'material {%s} doesn''t exist in pyne material lib' % material)
                print_near_match(material, mat_lib)
=======
       # loop over materials in library
       for key in mat_lib.iterkeys():  
           if material == key :
                d=d+1
                # get the material
                new_mat = mat_lib.get(key)[:]
                # copy attrs 'cos python is dumb
                copy_attrs(new_mat,mat_lib.get(key))
                
                # set the mcnp material number or fluka material name
                set_attrs(new_mat,d,code,flukamaterials_list)
                material_object_list.append(new_mat)
                break
           if mat_lib.keys().index(key) == len(mat_lib.keys())-1:	
                print('material {%s} doesn''t exist in pyne material lib' %material)
                print_near_match(material,mat_lib)
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
                exit()

    # check that there are as many materials as there are groups
    if d != len(mat_list):
<<<<<<< HEAD:get_tag_values_matlib.sh
        print "There are insuficient materials"
        exit()
=======
	print "There are insuficient materials"
	exit()
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py

    # return the list of material objects to write to file
    return material_object_list

<<<<<<< HEAD:get_tag_values_matlib.sh

def copy_attrs(material, material_from_lib):
    # copy attrs from lib to material
    for key in list(material_from_lib.attrs.keys()):
        material.attrs[key] = material_from_lib.attrs[key]

    material.density = material_from_lib.density
    material.mass = material_from_lib.mass
    material.atoms_per_mol = material_from_lib.atoms_per_mol

    return

""" 
function to print near matches to material name
"""


def print_near_match(material, material_library):
    for item in material_library.iterkeys():
        if (material.lower() in item.lower()) or (material.upper() in item.upper()):
            print "near matches to ", material, " are "
            print item
=======
def copy_attrs(material,material_from_lib):
    # copy attrs from lib to material
    for key in list(material_from_lib.attrs.keys()):
        material.attrs[key]=material_from_lib.attrs[key]

    material.density = material_from_lib.density
    material.mass  = material_from_lib.mass
    material.atoms_per_mol = material_from_lib.atoms_per_mol

    return
    
""" 
function to print near matches to material name
"""
def print_near_match(material,material_library):
    for item in material_library.iterkeys() :  
        if ( material.lower() in item.lower()) or (material.upper() in item.upper()) :
	    print "near matches to ", material, " are " 
	    print item
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
            return


"""
function to set the attributes of the materials:
"""
<<<<<<< HEAD:get_tag_values_matlib.sh


def set_attrs(mat, number, code, flukamat_list):
    if code is 'mcnp' or 'both':
        mat.attrs['mat_number'] = str(number)
    if code == 'fluka' or 'both':
        fluka_material_naming(mat, number, flukamat_list)
=======
def set_attrs(mat,number,code,flukamat_list):
    if code is 'mcnp' or 'both' :     
        mat.attrs['mat_number']=str(number)
    if code == 'fluka' or 'both' :   
        fluka_material_naming(mat,number,flukamat_list)
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    return

"""
Function to prepare fluka material names:
"""
<<<<<<< HEAD:get_tag_values_matlib.sh


def fluka_material_naming(matl, number, flukamat_list):
    matf = matl.attrs['name']
    matf = ''.join(c for c in matf if c.isalnum())

    if len(matf) > 8:
        matf = matf[0:7]
=======
def fluka_material_naming(matl,number,flukamat_list):
    matf=matl.attrs['name']
    matf=''.join(c for c in matf if c.isalnum())

    if len(matf) > 8 :
        matf=matf[0:7]
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    else:
        pass

    # if name in list change name by appending number
<<<<<<< HEAD:get_tag_values_matlib.sh
    if matf.upper() in flukamat_list:
        if number <= 9:
            matf = matf.rstrip(matf[-1])
            matf = matf + str(number)
        else:
            for i in range(2):
                matf = matf.rstrip(matf[-1])
                matf = matf + str(number)

        flukamat_list.append(matf.upper())
    # otherwise uppercase
    else:
        flukamat_list.append(matf.upper())

    matl.attrs['fluka_name'] = matf.upper()
=======
    if matf.upper() in flukamat_list :
        if number <= 9 :
            matf=matf.rstrip(matf[-1])
            matf=matf+str(number)
        else:
            for i in range(2) :
                matf=matf.rstrip(matf[-1])
                matf=matf+str(number)

        flukamat_list.append(matf.upper())   
    # otherwise uppercase
    else :    
        flukamat_list.append(matf.upper())

    matl.attrs['fluka_name']=matf.upper()   
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
    return

"""
Function write_mats, writes material objects to hdf5 file
-------
material_list: list of PyNE Material Objects
filename: filename to write the objects to
"""


def write_mats_h5m(materials_list, filename):
    for material in materials_list:
        material.write_hdf5(filename)

"""
function to parse the script, adding options:
defining 
-f  : the .h5m file path
-d  : nuc_data path
"""


def parsing(parsescript):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', action='store', dest='datafile', help='the path to the .h5m file')
    parser.add_argument('-d', action='store', dest='nuc_data',
                        help='the path to the PyNE materials library   nuc_data.h5')
    parser.add_argument('-c', action='store', dest='code',
                        help='the format of the output h5m file; mcnp or fluka')
    parser.add_argument(
        '-o', action='store', dest='output', help='the name of the output file')
    args = parser.parse_args()
    if args.datafile:
<<<<<<< HEAD:get_tag_values_matlib.sh
        # print args
        global datafile
        datafile = args.datafile
    else:
        print('h5m file not found!!. [-f] is not set')
        exit()
    if args.nuc_data:
        global nuc_data
        nuc_data = args.nuc_data
    else:
        print('nuc_data file not found!!. [-d] is not set')
        exit()
    if args.code:
        global code
        code = args.code
    else:
        print('output file format is not specified!!. [-c] is not set')
        exit()
    if args.output:
        output = args.output
    else:
        output = 'output.h5m'

    return output

# parse the script
output = parsing(1)
=======
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
       output=args.output
    else :
       output='output.h5m'

    return output

#parse the script
output = parsing(1)            

>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
# get list of tag values
tag_values = get_tag_values(datafile)
# now load material library
mat_lib = load_mat_lib(nuc_data)
# check that material tags exist in library
# material_list is list of pyne objects in problem
mat_list = check_matname(tag_values)
# create material objects from library
<<<<<<< HEAD:get_tag_values_matlib.sh
material_object_list = check_and_create_materials(mat_list, mat_lib)
#print material_object_list
# write materials to file
write_mats_h5m(material_object_list, output)
=======
material_object_list = check_and_create_materials(mat_list,mat_lib)
# write materials to file
write_mats_h5m(material_object_list,output)
>>>>>>> b3d06a6cf82c7ae1a7805a1e2d009c1051a35628:get_tag_values_matlib.py
exit()
 
