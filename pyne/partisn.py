#!/usr/bin/env python

""" Module for the production of PartiSn input decks. PartiSn is a discrete
ordinates code produced by Los Almos National Laboratory (LANL). Can be used
to produce neutron, photon, or coupled neutron photon prblems, adjoint or
forward or time dependent problems can be run.

Module is designed to work on 1D, 2D, or 3D Cartesian geometries.

If PyTaps not installed then this module will not work.
"""

from __future__ import print_function, division
import sys
import collections
import string
import struct
import math
import os
import linecache
import datetime
from warnings import warn
from pyne.utils import QAWarning
import itertools
from sets import Set

import numpy as np
import tables

from pyne import dagmc
from pyne.material import Material
from pyne.material import MultiMaterial
from pyne.material import MaterialLibrary

from pyne import nucname
from pyne.binaryreader import _BinaryReader, _FortranRecord

warn(__name__ + " is not yet QA compliant.", QAWarning)

# Mesh specific imports
try:
    from itaps import iMesh
    HAVE_PYTAPS = True
except ImportError:
    warn("the PyTAPS optional dependency could not be imported. "
                  "All aspects of the partisn module are not imported.",
                  QAWarning)
    HAVE_PYTAPS = False

if HAVE_PYTAPS:
    from pyne.mesh import Mesh, StatMesh, MeshError, IMeshTag


def write_partisn_input(mesh, hdf5, ngroup, nmq, **kwargs):
    """This definition reads all necessary attributes from a material-laden 
    geometry file, a pre-made PyNE mesh object, and the nuclear data cross 
    section library, and any optional inputs that are necessary for creating a 
    PARTISN input file. It then writes a PARTISN text input file for blocks 1-5.
    Note that comments appear in the created input file where more variables
    must be set. 
    
    Notes:
        This does not write out all necessary inputs for the solver and cross
        section library (block 3 and 5). There is no assumed cross section 
        library type.
    
    Parameters:
    -----------
        mesh : PyNE mesh object, a premade mesh object that conforms to the 
            geometry. Bounds of the mesh must correspond to the desired PartiSn
            fine mesh intervals. One fine mesh per coarse mesh will be created. 
            Can be 1-D, 2-D, or 3-D.
            Note: Only Cartesian meshes are currently supported.
        hdf5 : file path, a material-laden dagmc geometry file.
        ngroup : int, the number of energy groups in the cross section library
        nmq : int, the number of moments in a P_n expansion of the source
        
    Optional Parameters:
    --------------------
        data_hdf5path : str, the path in the heirarchy to the data table in an 
            HDF5 file. (for MaterialLibrary)
            default = material_library/materials
        nuc_hdf5path : str, the path in the heirarchy to the nuclide array in 
            an HDF5 file. (for MaterialLibrary)
            default = material_library/nucid
        names_dict : dict, pyne element/isotope names to bxslib name assignment,
            keys are pyne nucids (int) and values are bxslib names (str)
                Example: names_dict[250550000] ='mn55'
        input_file : str, desired name of generated PARTISN input file
            If no name is provided the name will default to 
            '<hdf5 file name>_partisn.inp'. Any file already existing by the 
            same name will be overwritten.
        num_rays : (for discretize_geom), int, optional, default = 10
            Structured mesh only. The number of rays to fire in each mesh row 
            for each direction.
        grid : (for discretize_geom), boolean, optional, default = False
            Structured mesh only. If false, rays starting points are chosen 
            randomly (on the boundary) for each mesh row. If true, a linearly 
            spaced grid of starting points is used, with dimension 
            sqrt(num_rays) x sqrt(num_rays). In this case, "num_rays" must be 
            a perfect square.
    
    Output:
    -------
        PARTISN Input file named by 'input_file' above or the default name.
            Note: read comments generated in file. Not all variables will be 
            assigned that are necessary.
    """
    
    # Load the geometry
    dagmc.load(hdf5)
    
    # Initialize dictionaries for each PARTISN block
    block01 = {}
    block02 = {}
    block03 = {}
    block04 = {}
    block05 = {}
    
    # Read optional inputs:
    
    # discretize_geom inputs
    if 'num_rays' in kwargs:
        num_rays = kwargs['num_rays']
    else:
        num_rays = 10
    
    if 'grid' in kwargs:
        grid = kwargs['grid']
    else:
        grid = False
    
    # hdf5 paths
    if 'data_hdf5path' in kwargs:
        data_hdf5path = kwargs['data_hdf5path']  
    else:
        data_hdf5path = '/material_library/materials'
    
    if 'nuc_hdf5path' in kwargs:
        nuc_hdf5path = kwargs['nuc_hdf5path']
    else:
        nuc_hdf5path = '/material_library/nucid'
    
    # input file name
    if 'input_file' in kwargs:
        input_file = kwargs['input_file']
        input_file_tf = True
    else:
        input_file_tf = False
    
    # Dictionary of hdf5 names and cross section library names
    # Assumes PyNE naming convention in the cross section library if no dict
    # provided.
    if 'names_dict' in kwargs:
        nuc_names = kwargs['names_dict']
        mat_lib = _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, nuc_names=nuc_names)
        mat_xs_names = _nucid_to_xs(mat_lib, nuc_names=nuc_names)
    else:
        mat_lib = _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path)
        mat_xs_names = _nucid_to_xs(mat_lib)
    
    # Set input variables
    
    block04['matls'] = mat_xs_names
    
    xs_names = _get_xs_names(mat_xs_names)
    block01['niso'] = len(xs_names)
    block03['names'] = xs_names

    block01['igeom'], bounds = _get_coord_sys(mesh)
    block01['ngroup'] = ngroup
    block01['mt'] = len(mat_lib)
    
    block02['zones'], zones = _get_zones(mesh, hdf5, bounds, num_rays, grid)
    block01['nzone'] = len(zones)
    block04['assign'] = zones
    
    for key in bounds.keys():
        if key == 'x':
            n = len(bounds[key]) - 1
            block01['im'] = n
            block01['it'] = block01['im']
            block02['xmesh'] = bounds[key]
            block05['sourcx'] = np.zeros(shape=(n, nmq), dtype=float)
            block05['sourcx'][:,0] = 1.0
        elif key == 'y':
            n = len(bounds[key]) - 1
            block01['jm'] = n
            block01['jt'] = block01['jm']
            block02['ymesh'] = bounds[key]
            block05['sourcy'] = np.zeros(shape=(n, nmq), dtype=float)
            block05['sourcy'][:,0] = 1.0
        elif key == 'z':
            n = len(bounds[key]) - 1
            block01['km'] = n
            block01['kt'] = block01['km']
            block02['zmesh'] = bounds[key]
            block05['sourcz'] = np.zeros(shape=(n, nmq), dtype=float)
            block05['sourcz'][:,0] = 1.0

    block05['source'] = np.zeros(shape=(ngroup, nmq), dtype=float)
    block05['source'][:,0] = 1.0
    
    # create title
    title = _title(hdf5)
    
    # call function to write to file
    if input_file_tf:
        _write_input(title, block01, block02, block03, block04, block05, name=input_file)
    else:
        _write_input(title, block01, block02, block03, block04, block05)
  

def _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, **kwargs):
    """Read material properties from the loaded dagmc geometry.
    """
    
    # If a set of nuc_names is provided, then collapse elements
    if 'nuc_names' in kwargs:
        nuc_names = kwargs['nuc_names']
        collapse = True
        # set of exception nuclides for collapse_elements
        mat_except = Set(nuc_names.keys())
    else:
        collapse = False
    
    # collapse isotopes into elements (if required)
    mats = MaterialLibrary(hdf5, datapath=data_hdf5path, nucpath=nuc_hdf5path)
    mats_collapsed = {}
    for mat_name in mats.keys():
        if collapse:
            mats_collapsed[mat_name] = mats[mat_name].collapse_elements(mat_except)
        else:
            mats_collapsed[mat_name] = mats[mat_name]

    # convert mass fraction to atom density in units [at/b-cm]
    mat_lib = {}
    for mat_name, comp in mats_collapsed.iteritems():
        atom_dens_dict = comp.to_atom_dens()
        comp_list = {}
        for nucid, dens in atom_dens_dict.iteritems():
            # convert from [at/cc] to [at/b-cm]
            comp_list[nucid] = dens*10.**-24
        mat_lib[mat_name] = comp_list

    return mat_lib


def _nucid_to_xs(mat_lib, **kwargs):
    """Replace nucids with xs library names.
    """
    
    if 'nuc_names' in kwargs:
        nuc_names = kwargs['nuc_names']
        names_tf = True
    else:
        names_tf = False
    
    mat_xs_names = {}
    for mat in mat_lib.keys():
        mat_xs_names[mat] = {}
        for nucid in mat_lib[mat].keys():
            
            if names_tf:
                if nucid in nuc_names.keys():
                    name = nuc_names[nucid]
                    mat_xs_names[mat][name] = mat_lib[mat][nucid]
                else:
                    warn("Nucid {0} does not exist in the provided nuc_names dictionary.".format(nucid))
                    mat_xs_names[mat]["{0}".format(nucid)] = mat_lib[mat][nucid]
            else:
                mat_xs_names[mat][nucname.name(nucid)] = mat_lib[mat][nucid]

    return mat_xs_names
    

def _get_xs_names(mat_xs_names):
    """Create list of names (strings) of the nuclides that appear in the cross
    section library from the list of nuc_names.
    """
    
    xs_names = []
    for mat, nuc_set in mat_xs_names.iteritems():
        for name in nuc_set.keys():
            if name not in xs_names:
                xs_names.append(name)
    
    return xs_names


def _get_coord_sys(mesh):
    """Determine coordinate system and get bounds
    """
    
    # get number of divisions
    nx = len(mesh.structured_get_divisions("x"))
    ny = len(mesh.structured_get_divisions("y"))
    nz = len(mesh.structured_get_divisions("z"))
    
    coord_sys = ""
    if nx > 2:
        coord_sys += "x"
    if ny > 2:
        coord_sys += "y"
    if nz > 2:
        coord_sys += "z"

    # collect values of mesh boundaries for each coordinate
    bounds = {}
    fine = {}
    for i in coord_sys:
        bounds[i] = mesh.structured_get_divisions(i)

    # Determine IGEOM
    # assumes a Cartesian system
    if len(coord_sys) == 1:
        igeom = 'SLAB'
    elif len(coord_sys) == 2:
        igeom = 'X-Y'
    elif len(coord_sys) == 3:
        igeom = 'X-Y-Z'
    
    return igeom, bounds


def _get_zones(mesh, hdf5, bounds, num_rays, grid):
    """Get the minimum zone definitions for the geometry.
    """
    
    # Descretize the geometry and get cell fractions
    dg = dagmc.discretize_geom(mesh, num_rays=num_rays, grid=grid)
    
    # Reorganize dictionary of each voxel's info with the key the voxel number 
    # and values of cell and volume fraction   
    voxel = {}
    for i in dg:
        idx = i[0]  # voxel number
        if idx not in voxel.keys():
            voxel[idx] = {}
            voxel[idx]['cell'] = []
            voxel[idx]['vol_frac'] = []
        voxel[idx]['cell'].append(i[1])
        voxel[idx]['vol_frac'].append(i[2])

    # get material to cell assignments
    mat_assigns = dagmc.materials_to_cells(hdf5)
    
    # Replace cell numbers with materials, eliminating duplicate materials
    # within single zone definition
    zones = {}
    for z in voxel.keys():
        zones[z] = {}
        zones[z]['vol_frac'] = []
        zones[z]['mat'] = []
        for i, cell in enumerate(voxel[z]['cell']):
            if mat_assigns[cell] not in zones[z]['mat']:
                # create new entry
                zones[z]['mat'].append(mat_assigns[cell])
                zones[z]['vol_frac'].append(voxel[z]['vol_frac'][i])
            else:
                # update value that already exists with new volume fraction
                for j, val in enumerate(zones[z]['mat']):
                    if mat_assigns[cell] == val:
                        vol_frac = zones[z]['vol_frac'][j] + voxel[z]['vol_frac'][i]
                        zones[z]['vol_frac'][j] = vol_frac
    
    # Eliminate duplicate zones and assign each voxel a zone number.
    # Assign zone = 0 if vacuum or graveyard and eliminate material definition.
    voxel_zone = {}
    zones_mats = {}
    z = 0
    match = False
    first = True    
    for i, vals in zones.iteritems():
        for zone, info in zones_mats.iteritems():
            if (vals['mat'] == info['mat']) and \
                    np.allclose(np.array(vals['vol_frac']), \
                                np.array(info['vol_frac']), rtol=1e-5):
                match = True
                y = zone
                break
            else:
                match = False
        if first or not match:
            if vals['mat'] in [['mat:Vacuum'], ['mat:vacuum'], 
                    ['mat:graveyard'], ['mat:Graveyard']]:
                voxel_zone[i] = 0
            else:
                z += 1
                zones_mats[z] = zones[i]
                voxel_zone[i] = z
                first = False
        else:
            if vals['mat'] in [['mat:Vacuum'], ['mat:vacuum'], 
                    ['mat:graveyard'], ['mat:Graveyard']]:
                voxel_zone[i] = 0
            else:
                voxel_zone[i] = y
    
    # Put zones into format for PARTISN input
    if 'x' in bounds.keys():
        im = len(bounds['x']) - 1
    else:
        im = 1
    
    if 'y' in bounds.keys():
        jm = len(bounds['y']) - 1
    else:
        jm = 1
    
    if 'z' in bounds.keys():
        km = len(bounds['z']) - 1
    else:
        km = 1

    n = 0
    zones_formatted = np.zeros(shape=(im, jm*km), dtype=int)
    for i in range(im):
        for jk in range(jm*km):
            zones_formatted[i,jk] = voxel_zone[n]
            n += 1
            
    return zones_formatted, zones_mats
    

def _title(hdf5):
    """Create a title for the input based on the geometry name.
    """
    if "/" in hdf5:
        name = hdf5.split("/")[len(hdf5.split("/"))-1].split(".")[0]
    else:
        name = hdf5.split(".")[0]
    
    return name


def _write_input(title, block01, block02, block03, block04, block05, **kwargs):
    """Write all variables and comments to a file.
    """
    
    # Create file to write to
    if 'name' in kwargs:
        f = open(kwargs['name'], 'w')
    else:
        file_name = str(title) + '_partisn.inp'
        f = open(file_name, 'w')
    
    # Write title
    f.write("     1     0     0\n")
    f.write(str(title))

    f.write("\n\ ")
    f.write("\n\ Notes: This input assumes a volumetric source calculation using")
    f.write("\n\ default PARTISN values in many cases. Please refer to the comments")
    f.write("\n\ throughout and the PARTISN input manual.")
    f.write("\n\ Variables that MUST be set in each block (other defaults and \n")
    f.write("\ optional variables may exist):")
    f.write("\n\     Block 1:  ISN")
    f.write("\n\     Block 3:  LIB, MAXORD, IHM, IHT")
    f.write("\n\     Block 6:  no input is provided for block 6")
    
    ###########################################
    #              Write Block 1              #
    ###########################################
    f.write("\n\ \n")
    f.write("\ ------------ Block 1 (Control and Dimensions) ------------")
    f.write("\n\ \n")
    f.write("igeom='{0}'".format(block01['igeom']))
    f.write("  ngroup={0}".format(block01['ngroup']))
    f.write("  niso={0}".format(block01['niso']))
    f.write("  mt={0}".format(block01['mt']))
    f.write("  nzone={0}\n".format(block01['nzone']))
    
    f.write("\ Please provide input for ISN variable:\n")
    f.write("\ isn=  \n")
    
    if 'im' in block01.keys():
        f.write("im={0}".format(block01['im']))
        f.write("  it={0}  ".format(block01['it']))
    if 'jm' in block01.keys():
        f.write("jm={0}".format(block01['jm']))
        f.write("  jt={0}  ".format(block01['jt']))
    if 'km' in block01.keys():
        f.write("km={0}".format(block01['km']))
        f.write("  kt={0}  ".format(block01['kt']))
    
    f.write("\n")
    f.write('t')
    
    ###########################################
    #              Write Block 2              #
    ###########################################
    f.write("\n\ \n")
    f.write("\ ------------ Block 2 (Geometry) ------------")
    f.write("\n\ \n")
    
    if 'xmesh' in block02.keys():
        f.write("xmesh= ")
        count = 0
        for i, val in enumerate(block02['xmesh']):
            count += 1
            f.write("{:.3f} ".format(val))
            if count == 8:
                if i != len(block02['xmesh'])-1:
                    f.write("\n       ")
                count = 0
        f.write("\nxints= ")
        f.write("{0}R {1}".format(len(block02['xmesh'])-1, 1))
        f.write("\n")
        
    if 'ymesh' in block02.keys():
        f.write("ymesh= ")
        count = 0
        for i, val in enumerate(block02['ymesh']):
            count += 1
            f.write("{:.3f} ".format(val))
            if count == 8:
                if i != len(block02['ymesh'])-1:
                    f.write("\n       ")
                count = 0
        f.write("\nyints= ")
        f.write("{0}R {1}".format(len(block02['ymesh'])-1, 1))
        f.write("\n")
        
    if 'zmesh' in block02.keys():
        f.write("zmesh= ")
        count = 0
        for i, val in enumerate(block02['zmesh']):
            count += 1
            f.write("{:.3f} ".format(val))
            if count == 8:
                if i != len(block02['zmesh'])-1:
                    f.write("\n       ")
                count = 0
        f.write("\nzints= ")
        f.write("{0}R {1}".format(len(block02['zmesh'])-1, 1))
        f.write("\n")
        
    f.write("zones= ")
    for i, row in enumerate(block02['zones']):
        string = format_repeated_vector(row)
        list_string = string.split()
        count = 0
        for num in list_string:
            f.write("{} ".format(num))
            count += 1
            if count == 20:
                f.write("\n       ")
                count = 0
        f.write(";")
        if i != len(block02['zones'])-1:
            f.write("\n       ")
        else:
            f.write("\n")

    f.write("t")
    
    ###########################################
    #              Write Block 3              #
    ###########################################
    f.write("\n\ \n")
    f.write("\ ------------ Block 3 (Nuclear Data) ------------")
    f.write("\n\ \n")
    f.write("\ Please provide input for the following variables:\n")
    f.write("\ lib=\n")
    f.write("\ maxord=\n")
    f.write("\ ihm=\n")
    f.write("\ iht=\n")
    
    f.write("\ Note: NAMES is not all inclusive. Only NAMES that are present in\n")
    f.write("\ meshed area are listed.\n")
    f.write("names= ")
    count = 0
    for i, name in enumerate(block03['names']):
        count += 1
        f.write("{0} ".format(name))
        if count == 10:
            if i != len(block03['names'])-1:
                f.write("\n       ")
            count = 0
    
    #f.write("\n\ \n")
    f.write("\nt")
    
    ###########################################
    #              Write Block 4              #
    ###########################################
    f.write("\n\ \n")
    f.write("\ ------------ Block 4 (Cross-Section Mixing) ------------")
    f.write("\n\ \n")
    
    f.write("matls= ")
    for i, mat in enumerate(block04['matls']):
        f.write("{0} ".format(mat))
        count = 0
        j = 0
        for iso, dens in block04['matls'][mat].iteritems():
            count += 1
            j += 1
            if j != len(block04['matls'][mat]):
                f.write("{} {:.4e}, ".format(iso, dens))
                if count == 3:
                    if j != len(block04['matls'][mat]):
                        f.write("\n       ")
                    count = 0
            else:
                if i == len(block04['matls']) - 1:
                    f.write("{} {:.4e};\n".format(iso, dens))
                else:
                    f.write("{} {:.4e};\n       ".format(iso, dens))
            
    
    f.write("assign= ")
    for i, z in enumerate(block04['assign']):
        f.write("{0} ".format(z))
        count = 0
        for j, mat in enumerate(block04['assign'][z]['mat']):
            count += 1
            if j != len(block04['assign'][z]['mat'])-1:
                f.write("{} {:.4e}, ".format(mat, block04['assign'][z]['vol_frac'][j]))
                if count == 3:
                    if i != len(block04['assign'][z]['mat'])-1:
                        f.write("\n          ")
                    count = 0
            else:
                if i == len(block04['assign']) - 1:
                    f.write("{} {:.4e};\n".format(mat, block04['assign'][z]['vol_frac'][j]))
                else:
                    f.write("{} {:.4e};\n        ".format(mat, block04['assign'][z]['vol_frac'][j]))
    
    f.write("t")
    
    ###########################################
    #              Write Block 5              #
    ###########################################
    f.write("\n\ \n")
    f.write("\ ------------ Block 5 (Solver Inputs) ------------")
    f.write("\n\ \n")
    f.write("\ This input assumes a volumetric source calculation with vacuum boundary conditions.\n")
    f.write("\ Change inputs below if otherwise.\n")
    f.write("ievt=0      \ source calculation\n")
    f.write("\ isct=     \ Legendre order of scattering (default=0)\n")
    f.write("\ ith=      \ 0/1/2= direct/adjoint/POI calculation (default=0)\n")
    f.write("\ ibl=      \ left BC (default=0, vacuum)\n")
    f.write("\ ibr=      \ right BC (default=0, vacuum)\n")
    f.write("\ ibt=      \ top BC (default=0, vacuum)\n")
    f.write("\ ibb=      \ bottom BC (default=0, vacuum)\n")
    f.write("\ ibfrnt=   \ front BC (default=0, vacuum)\n")
    f.write("\ ibback=   \ back BC (default=0, vacuum)\n")
    f.write("\ \n")
    
    f.write("\ Source is in format of option 3 according to PARTISN input manual.\n")
    f.write("\ Default is an evenly distributed volume source.\n")
    f.write("source= ")
    count = 0
    tot = 0
    for row in block05['source']:
        formatted_string = format_repeated_vector(row)
        f.write(formatted_string)
        tot += 1
        count += 1
        if count == 4:
            if tot != len(block05['source']):
                f.write(";\n        ")
            else:
                f.write(";")
            count = 0
        else:
            f.write("; ")
    f.write("\n")
        
    if 'sourcx' in block05.keys():
        f.write("sourcx= ")
        count = 0
        tot = 0
        for row in block05['sourcx']:
            formatted_string = format_repeated_vector(row)
            f.write(formatted_string)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcx']):
                    f.write(";\n        ")
                else:
                    f.write(";")
                count = 0
            else:
                f.write("; ")
        f.write("\n")
                
    if 'sourcy' in block05.keys():
        f.write("sourcy= ")
        count = 0
        tot = 0
        for row in block05['sourcy']:
            formatted_string = format_repeated_vector(row)
            f.write(formatted_string)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcy']):
                    f.write(";\n        ")
                else:
                    f.write(";")
                count = 0
            else:
                f.write("; ")
        f.write("\n")
            
    if 'sourcz' in block05.keys():
        f.write("sourcz= ")
        count = 0
        tot = 0
        for row in block05['sourcz']:
            formatted_string = format_repeated_vector(row)
            f.write(formatted_string)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcz']):
                    f.write(";\n        ")
                else:
                    f.write(";")
                count = 0
            else:
                f.write("; ")
        f.write("\n")


def format_repeated_vector(vector):
    """Creates string out of a vector with the PARTISN format for repeated
    numbers.
    
    Parameters:
    -----------
        vector: list of numbers, desired list to be formatted
    
    Returns:
    --------
        string: str, formatted string representation of the vector
    
    Example:
        vector = [1, 2, 0, 0, 0, 7, 8, 3, 3]
        string = "1 2 3R 0 7 8 2R 3"
    """
    
    # put vector into a list of lists formatted as 
    # [[number , R], [number, R], ...]
    # where 'R' is the number of times that 'number' is repeated
    tot = 0
    repeats = []
    for i, val in enumerate(vector):
        if tot == 0:
            repeats.append([val, 1])
            tot += 1
        else:
            if val == repeats[tot-1][0]:
                repeats[tot - 1][1] += 1
            else:
                repeats.append([val, 1])
                tot += 1
    
    # make into a string of characters
    string = ""
    n = 0
    for pair in repeats:
        if pair[1] == 1:
            string += "{} ".format(pair[0])
            n =+ 1
        else:
            string += "{0}R {1} ".format(pair[1], pair[0])
            n += 2

    return string
