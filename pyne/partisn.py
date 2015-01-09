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


def write_partisn_input(mesh, hdf5, ngroup, isn, nmq, **kwargs):
    """This class reads all necessary attributes from a material-laden 
    geometry file, a pre-made PyNE mesh object, and the nuclear data cross 
    section library, and any optional inputs that are necessary for creating a 
    PARTISN input file. It then writes a PARTISN text input file.
    
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
        isn : int, S_n order to be used MIGHT NOT KEEP AS AN INPUT
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
    
    Output:
    -------
        PARTISN Input file named by 'input_file' above or the default name.
            Note: read comments generated in file. Not all variables will be 
            assigned that are necessary.
    """
    
    # Load the geometry
    dagmc.load(hdf5)
    
    # Read optional inputs:
    
    # hdf5 paths
    if 'data_hdf5path' in kwargs:
        data_hdf5path = kwargs['data_hdf5path']  
    else:
        data_hdf5path = '/material_library/materials'
    
    if 'nuc_hdf5path' in kwargs:
        nuc_hdf5path = kwargs['nuc_hdf5path']
    else:
        nuc_hdf5path = '/material_library/nucid'
    
    # Dictionary of hdf5 names and cross section library names
    # Assumes PyNE naming convention in the cross section library if no dict
    # provided.
    if 'names_dict' in kwargs:
        nuc_names = kwargs['names_dict']
    else:
        # read a function
        pass
    
    if 'input_file' in kwargs:
        input_file = kwargs['input_file']
        input_file_tf = True
    else:
        input_file_tf = False
    
    # Initialize dictionaries for each PARTISN block
    block01 = {}
    block02 = {}
    block03 = {}
    block04 = {}
    block05 = {}
    
    # Set input variables
    
    block01['igeom'], bounds = _get_coord_sys(mesh)
    block01['ngroup'] = ngroup
    block01['isn'] = isn
    
    xs_names = _get_xs_names(nuc_names)
    block01['niso'] = len(xs_names)
    
    mat_lib = _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, nuc_names)
    block01['mt'] = len(mat_lib)
    
    block02['zones'], zones = _get_zones(mesh, hdf5, bounds)
    block01['nzone'] = len(zones)
    
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
    
    block03['names'] = xs_names
    
    mat_xs_names = _nucid_to_xs(mat_lib, xs_names, nuc_names)
    block04['matls'] = mat_xs_names
    
    block04['assign'] = zones
    
    block05['ievt'] = 0 # default? 0 = source
    block05['source'] = np.zeros(shape=(ngroup, nmq), dtype=float)
    block05['source'][:,0] = 1.0
    
    title = _title(hdf5)
    
    if input_file_tf:
        _write_input(title, block01, block02, block03, block04, block05, name=input_file)
    else:
        _write_input(title, block01, block02, block03, block04, block05)
    
    
def _get_xs_names(nuc_names):
    xs_names = []
    for nucid, name in nuc_names.iteritems():
        xs_names.append(name)
    
    return xs_names


def _get_coord_sys(mesh):
    
    # Determine coordinate system and get bounds
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
        IGEOM = 'SLAB'
    elif len(coord_sys) == 2:
        IGEOM = 'X-Y'
    elif len(coord_sys) == 3:
        IGEOM = 'X-Y-Z'
    
    return IGEOM, bounds


def _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, nuc_names):
    # reads material properties from the loaded dagmc_geometry
    
    # set of exception nuclides for collapse_elements
    mat_except = Set(nuc_names.keys())
    
    # collapse isotopes into elements
    mats = MaterialLibrary(hdf5, datapath=data_hdf5path, nucpath=nuc_hdf5path)
    mats_collapsed = {}
    for mat_name in mats.keys():
        mats_collapsed[mat_name] = mats[mat_name].collapse_elements(mat_except)      
    
    # convert mass fraction to atom density in units [at/b-cm]
    mat_lib = {}
    comp_list = {}
    for mat_name, comp in mats_collapsed.iteritems():
        atom_dens_dict = comp.to_atom_dens()
        for nucid, dens in atom_dens_dict.iteritems():
            # convert from [at/cc] to [at/b-cm]
            comp_list[nucid] = dens*10.**-24
        mat_lib[mat_name] = comp_list

    return mat_lib


def _get_zones(mesh, hdf5, bounds):
    
    # Descretize the geometry and get cell fractions
    dg = dagmc.discretize_geom(mesh)
    
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
            if mat_assigns[cell].strip(":")[1] not in zones[z]['mat']:
                # create new entry
                zones[z]['mat'].append(mat_assigns[cell].strip(":")[1])
                zones[z]['vol_frac'].append(voxel[z]['vol_frac'][i])
            else:
                # update value that already exists with new volume fraction
                for j, val in enumerate(zones[z]['mat']):
                    if mat_assigns[cell].strip(":")[1] == val:
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
                                np.array(info['vol_frac']), rtol=1e-8):
                match = True
                y = zone
                break
            else:
                match = False
        if first or not match:
            if vals['mat'] in [['Vacuum'], ['vacuum'], 
                    ['graveyard'], ['Graveyard']]:
                voxel_zone[i] = 0
            else:
                z += 1
                zones_mats[z] = zones[i]
                voxel_zone[i] = z
                first = False
        else:
            if vals['mat'] in [['Vacuum'], ['vacuum'], 
                    ['graveyard'], ['Graveyard']]:
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
    

def _nucid_to_xs(mat_lib, xs_names, nuc_names):

    # replace nucids with xs library names
    mat_xs_names = {}
    for mat in mat_lib.keys():
        mat_xs_names[mat] = {}
        for nucid in mat_lib[mat].keys():
            if nucid in nuc_names.keys():
                name = nuc_names[nucid]
                mat_xs_names[mat][name] = mat_lib[mat][nucid]
            else:
                warn("Nucid {0} does not exist in the hdf5 geometry.".format(nucid))
                mat_xs_names[mat]["{0}".format(nucid)] = mat_lib[mat][nucid]

    return mat_xs_names


def _title(hdf5):
    
    if "/" in hdf5:
        name = hdf5.split("/")[len(hdf5.split("/"))-1].split(".")[0]
    else:
        name = hdf5.split(".")[0]
    
    dt = datetime.datetime.now()
    
    title = [name, dt]
    
    return title


def _write_input(title, block01, block02, block03, block04, block05, **kwargs):
    
    # Create file to write to
    if 'name' in kwargs:
        f = open(kwargs['name'], 'w')
    else:
        file_name = str(title[0]) + '_partisn.inp'
        f = open(file_name, 'w')
    
    # Write title
    f.write("     1     0     0\n")
    f.write(str(title[0])+'  ')
    f.write(str(title[1]))
    
    # Write Block 1
    f.write("\n\ \n")
    f.write("\ ------------ Block 1 (Control and Dimensions) ------------")
    f.write("\n\ \n")
    f.write("igeom='{0}'".format(block01['igeom']))
    f.write("  ngroup={0}".format(block01['ngroup']))
    #f.write("  isn={0}".format(block01['isn']))
    f.write("  niso={0}".format(block01['niso']))
    f.write("  mt={0}".format(block01['mt']))
    f.write("  nzone={0}\n".format(block01['nzone']))
    
    
    
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
    
    f.write("\ Must provide a value for variable ISN\n")
    f.write("\ isn=\n")
    
    f.write("\ \n")
    f.write('t')
    
    # Write Block 2
    f.write("\n\ \n")
    f.write("\ ------------ Block 2 (Geometry) ------------")
    f.write("\n\ \n")
    
    if 'xmesh' in block02.keys():
        f.write("xmesh= ")
        count = 0
        for i in block02['xmesh']:
            count += 1
            f.write("{:.3f} ".format(i))
            if count == 8:
                f.write("\n       ")
                count = 0
        f.write("\nxints= ")
        f.write("{0}R {1}".format(len(block02['xmesh'])-1, 1))
        f.write("\n")
        
    if 'ymesh' in block02.keys():
        f.write("ymesh= ")
        count = 0
        for i in block02['ymesh']:
            count += 1
            f.write("{:.3f} ".format(i))
            if count == 8:
                f.write("\n       ")
                count = 0
        f.write("\nyints= ")
        f.write("{0}R {1}".format(len(block02['ymesh'])-1, 1))
        f.write("\n")
        
    if 'zmesh' in block02.keys():
        f.write("zmesh= ")
        count = 0
        for i in block02['zmesh']:
            count += 1
            f.write("{:.3f} ".format(i))
            if count == 8:
                f.write("\n       ")
                count = 0
        f.write("\nzints= ")
        f.write("{0}R {1}".format(len(block02['xmesh'])-1, 1))
        f.write("\n")
    
    f.write("zones= ")
    for row in block02['zones']:
        count = 0
        for i in row:
            count += 1
            f.write("{:3d} ".format(i))
            if count == 15:
                f.write("\n       ")
                count = 0
        f.write(";\n       ")
    
    f.write("\ \n")
    f.write("t")
    
    # Write Block 3
    f.write("\n\ \n")
    f.write("\ ------------ Block 3 (Nuclear Data) ------------")
    f.write("\n\ \n")
    f.write("\ The following variables must be set (optional variables are not set):")
    f.write("\ lib=\n")
    f.write("\ maxord=\n")
    f.write("\ ihm=\n")
    f.write("\ iht=\n")
    
    f.write("names= ")
    count = 0
    for name in block03['names']:
        count += 1
        f.write("{0} ".format(name))
        if count == 10:
                f.write("\n       ")
                count = 0
    
    f.write("\n\ \n")
    f.write("t")
    
    # Write Block 4
    f.write("\n\ \n")
    f.write("\ ------------ Block 4 (Cross-Section Mixing) ------------")
    f.write("\n\ \n")
    
    f.write("matls= ")
    for mat in block04['matls']:
        f.write("{0} ".format(mat))
        count = 0
        for iso, dens in block04['matls'][mat].iteritems():
            count += 1
            f.write("{} {:.4e}, ".format(iso, dens))
            if count == 3:
                f.write("\n       ")
                count = 0
        f.write(";\n       ")
    
    f.write("\n\ ")
    f.write("assign= ")
    for z in block04['assign']:
        f.write("{0} ".format(z))
        count = 0
        for i, mat in enumerate(block04['assign'][z]['mat']):
            count += 1
            f.write("{} {:.4e}, ".format(mat, block04['assign'][z]['vol_frac'][i]))
            if count == 3:
                f.write("\n       ")
                count = 0
        f.write(";\n        ")
