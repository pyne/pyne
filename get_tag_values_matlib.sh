#!/usr/bin/python

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
   # get all entities
    ents = mesh.getEntities(iBase.Type.all, iMesh.Topology.triangle)
    # get mesh set
    mesh_set = mesh.getEntSets()
    tag_values = []  # list of tag values
    found_all_tags = 0
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
                a = []
                # since we have a byte type tag loop over the 32 elements
                for part in tag:
                    # if the byte char code is non 0
                    if (part != 0):
                        # convert to ascii
                        a.append(str(unichr(part)))
                        # join to end string
                        test = ''.join(a)
                        # the the string we are testing for is not in the list of found
                        # tag values, add to the list of tag_values
                # if not already in list append to lilst
                if not any(test in s for s in tag_values):
                    tag_values.append(test)
                # last tag we are done
                if any('impl_complement' in s for s in tag_values):
                    found_all_tags = 1
    print('The matrials group names found in the h5m file are: ')
    print tag_values
    return tag_values

"""
function which loads pyne material library
"""


def load_mat_lib(filename):
    mat_lib = material.Material()
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(
        filename, datapath="/material_library/materials", nucpath="/material_library/nucid")
    return mat_lib

"""
function to check that materials exist in library
-------------------------------------------------
tag_values - list of tags from dagmc file
mat_lib - pyne material library instance
"""


def check_matname(tag_values):
    # loop over tags
    # a dictionary of material names as values and density as keys
    mat_dict = {}
    # loop over the tags in the file
    for tag in tag_values:
        # look for mat, this id's material in group name
        try:
            if "mat" in tag:
            # split on the basis of "/" being delimiter and split colons from
            # name
                if "/" in tag:
                    mat_name = tag.split("/")
                    # list of material name only
                    matname = mat_name[0].split(":")
                    matdensity = mat_name[1].split(":")
                    mat_dict[str(matdensity[1])] = matname[1]
                # otherwise we have only "mat:"
                else:
                    matname = tag.split(":")
                    mat_dict[""] = matname[1]
        except:
            print("Could not find group name in appropriate format"), tag
            exit()
    print mat_dict
    # error conditions, not tags found
    if len(mat_dict) == 0:
        print("no group names found")
        exit()
    # print mat_list
    return mat_dict

"""
----------------------------------
"""


def check_and_create_materials(material_dict, mat_lib):
    flukamaterials_list = []
    material_object_list = []
    d = 0  # counter of the materials to set mat_number for mcnp
    # loop over materials in geometry
    for dkey in mat_dict:
        material = material_dict[dkey]
        # loop over materials in library
        for key in mat_lib.iterkeys():
            if material == key:
                d = d + 1
                # get the material
                new_mat = mat_lib.get(key)[:]
                # copy attrs 'cos python is dumb
                copy_attrs(new_mat, mat_lib.get(key))
                flukamaterials_list.append(material)
                # set the mcnp material number or fluka material name
                set_attrs(new_mat, d, code, flukamaterials_list, dkey)
                material_object_list.append(new_mat)
                break
            if mat_lib.keys().index(key) == len(mat_lib.keys()) - 1:
                print(
                    'material {%s} doesn''t exist in pyne material lib' % material)
                print_near_match(material, mat_lib)

    # check that there are as many materials as there are groups
    if d != len(mat_dict):
        print "There are insuficient materials"
        exit()
    # return the list of material objects to write to file
    return material_object_list

"""
------------------------------
"""


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
    return

"""
function to set the attributes of the materials:
"""


def set_attrs(mat, number, code, flukamat_list, density):
    if code is 'mcnp' or 'both':
        mat.attrs['mat_number'] = str(number)
        mat.attrs['mat_density'] = density
    if code == 'fluka' or 'both':
        fluka_material_naming(mat, flukamat_list)
    return

"""
Function to prepare fluka material names:
"""


def fluka_material_naming(matl, flukamat_list):
    matf = matl.attrs['name']
    matf = ''.join(c for c in matf if c.isalnum())
    if len(matf) > 8:
        matf = matf[0:7]
    else:
        pass
    # if name in list change name by appending number
    if matf.upper() in flukamat_list:
        for a in range(len(flukamat_list)):
            a = a + 1
            if a <= 9:
                matf = matf.rstrip(matf[-1])
                matf = matf + str(a)
            else:
                for i in range(2):
                    matf = matf.rstrip(matf[-1])
                    matf = matf + str(a)
            if matf.upper() in flukamat_list:
                continue
            else:
                flukamat_list.append(matf.upper())
                break
    # otherwise uppercase
    else:
        flukamat_list.append(matf.upper())
    matl.attrs['fluka_name'] = matf.upper()
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
# get list of tag values
tag_values = get_tag_values(datafile)
# now load material library
mat_lib = load_mat_lib(nuc_data)
# check that material tags exist in library
# material_list is list of pyne objects in problem
mat_dict = check_matname(tag_values)
# create material objects from library
material_object_list = check_and_create_materials(mat_dict, mat_lib)
print material_object_list
# write materials to file
write_mats_h5m(material_object_list, output)
exit()
