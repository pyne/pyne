#!/usr/bin/python

import string
import argparse
try:
    from itaps import iMesh, iBase
except:
    raise ImportError("The PyTAPS dependency could not be imported")
try:
    from pyne import material
    from pyne.material import Material, MaterialLibrary
except:
    raise ImportError("The PyNE dependency could not be imported")

"""
function that gets all tags on dagmc geometry
------------------------------
filename : the dagmc filename
return vector of tag_values
"""


def get_tag_values(filename):
    mesh = iMesh.Mesh()
    mesh.load(filename)
    # get all entities
    ents = mesh.getEntities(iBase.Type.all, iMesh.Topology.triangle)
    # get mesh set
    mesh_sets = mesh.getEntSets()
    # tag_values = []  # list of tag values
    tag_values = []
    found_all_tags = 0
    for s in mesh_sets:
        if found_all_tags == 1:
            break
        # get all the tags
        tags = mesh.getAllTags(s)
        # loop over the tags checking for name
        for t in tags:
            # look for NAME tag
            if t.name == 'NAME':
                # the handle is the tag name
                t_handle = mesh.getTagHandle(t.name)
                # get the data for the tag, with taghandle in element i
                tag = t_handle[s]
                tag_to_script(tag, tag_values)
                # last tag we are done
                if any('impl_complement' in s for s in tag_values):
                    found_all_tags = 1
    print('The groups found in the h5m file are: ')
    print tag_values
    return tag_values

"""
function to transform the tags into strings
tag : string of the tag to add to tag_list
tag_list : vector of tags in the problem
returns tag_list
"""


def tag_to_script(tag, tag_list):
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
    # if not already in list append to list
    #    if not any(test in s for s in tag_list):
    # the original code was incorrectly missing groups when one of the same
    # name with/without rho was added
    if test not in tag_list:
        tag_list.append(test)
    return tag_list

"""
function which loads pyne material library
filename : string of the name of the material library
returns PyNE MaterialLibrary instance
"""


def load_mat_lib(filename):
    mat_lib = material.MaterialLibrary()
    mat_lib.from_hdf5(
        filename, datapath='/material_library/materials', nucpath='/material_library/nucid')
    return mat_lib

"""
function to check that material group names exist and creates
a list of names and densities if provided
-------------------------------------------------
tag_values - list of tags from dagmc file
mat_lib - pyne material library instance
returns mat_dens_list, a zipped pair of density and material name
"""


def check_matname(tag_values):
    # loop over tags
    g = 0
    mat_list_matname = []  # list of material names
    mat_list_density = []  # list of density if provided in the group names
    # loop over the tags in the file
    for tag in tag_values:
        if ('Graveyard' in tag) or ('graveyard' in tag):
            g = 1
            continue
    # look for mat, this id's material in group name
        if tag[0:3] == 'mat':
            # split on the basis of "/" being delimiter and split colons from
            # name
            if '/' in tag:
                splitted_group_name = mat_dens_split(tag)
            # otherwise we have only "mat:"
            elif ':' in tag:
                splitted_group_name = mat_split(tag)
            else:
                raise Exception(
                    "Couldn\'t find group name in appropriate format; ': is absent' in  %s" % tag)
            mat_list_matname.append(splitted_group_name['material'])
            mat_list_density.append(splitted_group_name['density'])
    if g == 0:
        raise Exception(
            "Graveyard group is missing! You must have a graveyard")
    mat_dens_list = zip(mat_list_matname, mat_list_density)
    # error conditions, no tags found
    if len(mat_dens_list) == 0:
        raise Exception(
            "No material group names found, you must have materials")

    return mat_dens_list


"""
function that splits group name on the basis of '/'
group name containing both material name and density
"""


def mat_dens_split(tag):
    splitted_group_name = {}
    mat_name = tag.split('/')
    if ':' not in mat_name[0]:
        raise Exception(
            "Couldn\'t find group name in appropriate format; ':' is absent in %s" % tag)
    # list of material name only
    matname = mat_name[0].split(':')
    if len(matname) > 2:
        raise Exception(
            "Wrong format for group names! %s. correct: mat:NAME/rho:VALUE or mat:NAME" % tag)
    if matname[1] == '':
        raise Exception(
            "Couldn\'t find group name in appropriate format; wrong material name in %s" % tag)
    splitted_group_name['material'] = matname[1]
    if mat_name[1] == '':
        raise Exception(
            "Couldn\'t find group name in appropriate format; extra \'/\' in %s" % tag)
    if ':' not in mat_name[1]:
        raise Exception(
            "Couldn\'t find group name in appropriate format; ':' is absent after the '/' in %s" % tag)
    matdensity = mat_name[1].split(':')
    try:
        matdensity_test = float(matdensity[1])
    except:
        raise Exception(
            "Couldn\'t find density in appropriate format!; density is not a float in %s" % tag)
    splitted_group_name['density'] = matdensity[1]
    return splitted_group_name

"""
function that splits group name on the basis of ':'
group name containing only material name
"""


def mat_split(tag):
    splitted_group_name = {}
    matname = tag.split(':')
    if len(matname) > 2:
        raise Exception(
            "Wrong format for group names! %s. correct: mat:NAME/rho:VALUE or mat:NAME" % tag)
    if matname[1] == '':
        raise Exception(
            "Couldn\'t find group name in appropriate format; wrong material name in %s" % tag)
    splitted_group_name['material'] = matname[1]
    splitted_group_name['density'] = ''
    return splitted_group_name

"""
function that checks the existence of material names on the PyNE library 
and creates a list of materials with attributes set
-------------------------------------------------------------
material_list : vector of material_name & density pairs
mat_lib : PyNE Material library object
"""


def check_and_create_materials(material_list, mat_lib):
    flukamaterials_list = []
    material_object_list = []
    d = 0
    # loop over materials in geometry
    for g in range(len(material_list)):
        material = material_list[g][0]
        # loop over materials in library
        for key in mat_lib.iterkeys():
            if material == key:
                d = d + 1
                # get the material
                new_mat = mat_lib.get(key)[:]
                flukamaterials_list.append(material)
                copy_metadata(new_mat, mat_lib.get(key))
                # set the mcnp material number and fluka material name
                set_metadata(new_mat, d, flukamaterials_list)

                # rename the material to match the group
                group_name = "mat:" + material_list[g][0]
                if material_list[g][1] is not '':
                    group_name += "/rho:" + material_list[g][1]
                print "grp2", group_name
                new_mat.metadata['name'] = group_name

                if material_list[g][1] != '':
                    new_mat.density = float(material_list[g][1])

                material_object_list.append(new_mat)
                break
            if mat_lib.keys().index(key) == len(mat_lib.keys()) - 1:
                print(
                    'Material {%s} doesn\'t exist in pyne material lib' % material)
                print_near_match(material, mat_lib)
                raise Exception(
                    'Couldn\'t find exact match in material library for : %s' % material)

    # check that there are as many materials as there are groups
    if d != len(material_list):
        raise Exception("There are insuficient materials")

    # return the list of material objects to write to file
    print material_object_list
    return material_object_list

"""
function to copy the metadata of materials from the PyNE material library
-------------------------------------
material : PyNE material object to copy data into 
material_from_lib : PyNE material objec to copy data from
"""


def copy_metadata(material, material_from_lib):
    # copy metadata from lib to material
    for key in list(material_from_lib.metadata.keys()):
        material.metadata[key] = material_from_lib.metadata[key]

    material.density = material_from_lib.density
    material.mass = material_from_lib.mass
    material.atoms_per_molecule = material_from_lib.atoms_per_molecule
    material.comp = material_from_lib.comp
    return material


"""
function to set the attributes of the materials:
----------------------------------------
mat : PyNE Material Object
number : mcnp material number
flukamat_list : 
returns : PyNE Material Object
"""


def set_metadata(mat, number, flukamat_list):
    mat.metadata['mat_number'] = str(number)
    mat.metadata['fluka_material_index'] = str(number + 25)
    fluka_material_naming(mat, flukamat_list)
    return mat


"""
Function to prepare fluka material names:
"""


def fluka_material_naming(material, flukamat_list):
    matf = material.metadata['name']
    matf = ''.join(c for c in matf if c.isalnum())
    if len(matf) > 8:
        matf = matf[0:8]
    else:
        pass
    # if name is in list, change name by appending number
    if matf.upper() in flukamat_list:
        for a in range(len(flukamat_list)):
            a = a + 1
            if a <= 9:
                matf = matf.rstrip(matf[-1])
                matf = matf + str(a)
            else:
                for i in range(len(a)):
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
    material.metadata['fluka_name'] = matf.upper()
    return material

""" 
function to print near matches to material name
"""


def print_near_match(material, material_library):
    list_of_matches = []
    for item in material_library.iterkeys():
        if (material.lower() in item.lower()) or (material.upper() in item.upper()):
            print("Near matches to %s are :" % material)
            print item
        list_of_matches.append(item)
    return list_of_matches

"""
Function that writes material objects to hdf5 file
-------
material_list: list of PyNE Material Objects
filename: filename to write the objects to
"""


def write_mats_h5m(materials_list, filename):
    new_matlib = MaterialLibrary()
    for material in materials_list:
        # using fluka name as index since this is unique
        new_matlib[material.metadata['name']] = material
    new_matlib.write_hdf5(filename)

"""
function to parse the script, adding options:
defining 
-f  : the .h5m file path
-d  : nuc_data path
-o  : name of the output h5m file "NAME.h5m"
"""


def parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', action='store', dest='datafile', help='The path to the .h5m file')
    parser.add_argument('-d', action='store', dest='nuc_data',
                        help='The path to the PyNE materials library nuc_data.h5')
    parser.add_argument(
        '-o', action='store', dest='output', help='The name of the output file ***.h5m')
    args = parser.parse_args()
    if not args.datafile:
        raise Exception('h5m file path not specified!!. [-f] is not set')
    if not args.nuc_data:
        raise Exception('nuc_data file path not specified!!. [-d] is not set')
    if not args.output:
        args.output = 'output.h5m'
    return args

"""
main
"""


def main():
    # parse the script
    args = parsing()
    # get list of tag values
    tag_values = get_tag_values(args.datafile)
    # now load material library
    mat_lib = load_mat_lib(args.nuc_data)
    # check that material tags exist in library # material_list is list of
    # pyne objects in problem
    mat_dens_list = check_matname(tag_values)
    # create material objects from library
    material_object_list = check_and_create_materials(
        mat_dens_list, mat_lib)
    # write materials to file
    write_mats_h5m(material_object_list, args.output)

if __name__ == "__main__":
    main()
