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
                  "All aspects of the PartiSn module are not imported.",
                  VnVWarning)
    HAVE_PYTAPS = False

if HAVE_PYTAPS:
    from pyne.mesh import Mesh, StatMesh, MeshError, IMeshTag


#class PartisnRead(object):
""" This class reads all necessary attributes from a material-laden 
geometry file, a pre-made PyNE mesh object, and the nuclear data cross 
section library, and any optional inputs that are necessary for creating a 
PARTISN input file. Supported are 1D, 2D, and 3D geometries.

Parameters
----------
    mesh : mesh object, a premade mesh object that conforms to the geometry. 
        Bounds of the mesh must correspond to the desired PartiSn fine mesh. 
        One fine mesh per coarse mesh will be created. Can be 1-D, 2-D, or 3-D.
        Only Cartesian based geometries are currently supported.
    hdf5 : file, a material-laden dagmc geometry file.
    nucdata : file, nuclear data cross section library.
                note: only BXSLIB format is currently supported.
    nuc_names : dict, pyne element/isotope names to bxslib name assignment,
                keys are pyne nucids (int) and values are bxslib names (str)
    datapath : str, optional, The path in the heirarchy to the data table 
            in an HDF5 file. (for MaterialLibrary)
                default = material_library/materials
    nucpath : str, optional, The path in the heirarchy to the 
            nuclide array in an HDF5 file. (for MaterialLibrary)
                default = material_library/nucid
"""

def read_hdf5_mesh(mesh, hdf5, nucdata, nuc_names, **kwargs):
    dagmc.load(hdf5)
    # optional inputs
    datapath = kwargs['datapath'] if 'datapath' in kwargs else '/material_library/materials'
    nucpath = kwargs['nucpath'] if 'nucpath' in kwargs else '/material_library/nucid'

    # get coordinate system and mesh bounds from mesh       
    # not necessary
    coord_sys, bounds = _read_mesh(mesh)
    
    # Read the materials from the hdf5 and convert to correct naming convention
    mat_lib = _get_materials(hdf5, datapath, nucpath, nuc_names)
    
    # Assign materials to cells   
    mat_assigns = _materials_to_cells(hdf5)
    
    # determine the zones
    zones, voxel_zone = _define_zones(mesh, mat_assigns)
    
    # read nucdata
    xs_names = _read_bxslib(nucdata)

    return coord_sys, bounds, mat_lib, zones, voxel_zone, xs_names
 
    
def _read_mesh(mesh):
    # determines the system geometry (1-D, 2-D, or 3-D Cartesian)
    # currently cartesian is only supported
    nx = len(mesh.structured_get_divisions("x"))
    ny = len(mesh.structured_get_divisions("y"))
    nz = len(mesh.structured_get_divisions("z"))
    
    # Check for dimensions with >1 voxel (>2 bounds)
    # This determines 1-D, 2-D, or 3-D
    
    # !!! change coord_sys to string "xyz" etc
    
    dim = 0
    i = False
    j = False
    k = False
    if nx > 2:
        dim += 1
        i = "x"
    if ny > 2:
        dim += 1
        if not i:
            i = "y"
        else:
            j = "y"
    if nz > 2:
        dim += 1
        if not i:
            i = "z"
        elif not j:
            j = "z"
        else:
            k = "z"
    
    # coordinate system data
    if dim == 1:
        coord_sys = [i]
    elif dim == 2:
        coord_sys = [i, j]
    elif dim == 3:
        coord_sys = [i, j, k]
        
    # collect values of mesh boundaries for each coordinate
    bounds = {}
    fine = {}
    for i in coord_sys:
        bounds[i] = mesh.structured_get_divisions(i)
        #fine[i] = [1]*(len(bounds[i]) - 1)
    
    return coord_sys, bounds

 
def _read_bxslib(nucdata):
    # read entire file
    binary_file = _BinaryReader(nucdata, mode='rb')
    record = binary_file.get_fortran_record()
    #print(binary_file)
    ##print(record.get_double())
    #print(record.get_string(28))   
        
    bxslib = open(nucdata, 'rb')
    string = ""
    edits = ""
    xs_names=[]
    # 181st byte is the start of xsnames
    bxslib.seek(180)
    done = False
    while not done:
        for i in range(0,8):
            bytes = bxslib.read(1)
            pad1=struct.unpack('s',bytes)[0]
            if '\x00' in pad1:
                done = True
                return xs_names
            string += pad1
        xs_names.append(string.strip(" "))
        string=""
    

def _get_materials(hdf5, datapath, nucpath, nuc_names):
    # reads material properties from the loaded dagmc_geometry
    
    # set of exception nuclides for collapse_elements
    mat_except = Set(nuc_names.keys())
    
    # collapse isotopes into elements
    mats = MaterialLibrary(hdf5,datapath=datapath,nucpath=nucpath)
    mats_collapsed = {}
    for mat_name in mats.keys():
        mats_collapsed[mat_name] = mats[mat_name].collapse_elements(mat_except)
    
    # Check that the materials are valid:
    #   1) non zero and non-negative densities (density = True)
    #   2) set of nuclides is not empty (else it is vacuum) (empty = False)
    #   3) nucids appear in nuc_names
    # might put 2 and 3 later      
    
    # convert mass fraction to atom fraction and then to [at/b-cm]
    Na = 6.022*(10.**23) # Avagadro's number [at/mol]
    barn_conv = 10.**-24 # [cm^2/b]
    mat_lib = {}
    for mat_name, comp in mats_collapsed.iteritems():
        #print(comp)
        comp_atom_frac = comp.to_atom_frac() # atom fractions
        density = comp.mass_density() # [g/cc]
        
        if density < 0.0:
            warn("Material {0} has an invalid negative density.".format(mat_name))
        
        mol_mass = comp.molecular_mass() # [g/mol]
        comp_list = {}
        
        for nucid, frac in comp_atom_frac.iteritems():
            comp_list[nucid] = frac*density*Na*barn_conv/mol_mass # [at/b-cm]
        
        mat_lib[mat_name] = comp_list

    return mat_lib


def _materials_to_cells(hdf5):
    """Takes the material-laden geometry and matches cells to materials
    """
    # Load the geometry
    dag_geom = iMesh.Mesh()
    dag_geom.load(hdf5)
    dag_geom.getEntities()
    mesh_sets = dag_geom.getEntSets()

    # Get tag handle
    cat_tag = dag_geom.getTagHandle('CATEGORY')
    id_tag = dag_geom.getTagHandle('GLOBAL_ID')
    name_tag = dag_geom.getTagHandle('NAME')

    # Get list of materials and list of cells
    mat_assigns={}
    
    # loop over all mesh_sets in model
    for mesh_set in mesh_sets:
        tags = dag_geom.getAllTags(mesh_set)
        
        # check for mesh_sets that are groups
        if name_tag in tags and cat_tag in tags \
                and _tag_to_string(cat_tag[mesh_set]) == 'Group':
            child_sets = mesh_set.getEntSets()
            name = _tag_to_string(name_tag[mesh_set])
            
            # if mesh_set is a group with a material name_tag, loop over child
            # mesh_sets and assign name to cell
            if 'mat:' in name:
                for child_set in child_sets:
                    child_tags = dag_geom.getAllTags(child_set)
                    if id_tag in child_tags:
                        cell = id_tag[child_set]
                        mat_assigns[cell] = name
                        
    return mat_assigns


def _tag_to_string(tag):
    a = []
    # since we have a byte type tag loop over the 32 elements
    for part in tag:
        # if the byte char code is non 0
        if (part != 0):
            # convert to ascii
            a.append(str(unichr(part)))
            # join to end string
            string = ''.join(a)
    return string


def _define_zones(mesh, mat_assigns):
    """This function takes results of discretize_geom and finds unique voxels
    """
    
    dg = dagmc.discretize_geom(mesh)
    
    # Create dictionary of each voxel's info    
    voxel = {}
    for i in dg:
        idx = i[0]
        if idx not in voxel.keys():
            voxel[idx] = {}
            voxel[idx]['cell'] = []
            voxel[idx]['vol_frac'] = []
            #voxel[idx]['rel_error'] = []
            
        voxel[idx]['cell'].append(i[1])
        voxel[idx]['vol_frac'].append(i[2])
        #voxel[idx]['rel_error'].append(i[3])

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
                                np.array(info['vol_frac']), rtol=1e-8):
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
    
    return zones_mats, voxel_zone


#class PartisnWrite(object):

def write_partisn_input(coord_sys, bounds, mat_lib, zones, voxel_zone, xs_names, nuc_names, ngroup, isn, nmq, hdf5, input_file):
    """This function writes out the necessary information to a text partisn 
    input file.
    
    Parameters
    ----------
        coord_sys : list of str, indicator of either 1-D, 2-D, or 3-D Cartesian
            geometry. 
                1-D: [i]
                2-D: [i ,j]
                3-D: [i, j, k]
                where i, j, and k are either "x", "y", or "z".
        bounds : dict of list of floats, coarse mesh bounds for each dimension. 
            Dictionary keys are the dimension "x", "y" or "z" for Cartesian. 
            Must correspond to a 1:1 fine mesh to coarse mesh interval.
        mat_lib : dict of dicts, keys are names of PyNE materials whose keys 
            are bxslib names and their value is atomic density in units 
            [at/b-cm].
        zones : dict of dict of lists, first dict key is PartiSn zone number 
            (int). Inner dict keys are "cell" with a list of cell numbers (int) 
            as values and "vol_frac" with a list of corresponding cell volume 
            fractions (float).
        xs_names : list of str, names of isotope/elements from the bxslib
    
    """
    title = _title(hdf5)
    
    block01 = _block01(coord_sys, xs_names, mat_lib, zones, bounds, ngroup, isn)
    #print(block01)
    
    block02 = _block02(bounds, voxel_zone)
    #print(block02)
    
    block03 = _block03(xs_names)
    #print(block03)
    
    block04 = _block04(mat_lib, xs_names, nuc_names, zones)
    #print(block04)
    
    block05 = _block05(ngroup, bounds, nmq)
    #print(block05)
    
    _write(title, block01, block02, block03, block04, block05)


def _title(hdf5):
    
    if "/" in hdf5:
        name = hdf5.split("/")[len(hdf5.split("/"))-1].split(".")[0]
    else:
        name = hdf5.split(".")[0]
    
    dt = datetime.datetime.now()
    
    title = [name, dt]
    
    return title

        
def _block01(coord_sys, xs_names, mat_lib, zones, bounds, ngroup, isn):
    block01 = {}
    
    # !!!! do get_coord_sys here and pull bounds here (not in read step)
    
    # Determine IGEOM
    if len(coord_sys) == 1:
        block01['IGEOM'] = 'SLAB'
    elif len(coord_sys) == 2:
        block01['IGEOM'] = 'X-Y' # assuming cartesian
    elif len(coord_sys) == 3:
        block01['IGEOM'] = 'X-Y-Z' # assuming cartesian
    
    block01['NGROUP'] = ngroup
    block01['ISN'] = isn
    
    # !!! ISN - have to read from bxslib still
    
    block01['NISO'] = len(xs_names)
    block01['MT'] = len(mat_lib)
    block01['NZONE'] = len(zones)
    
    # Number of Fine and Coarse Meshes
    # one fine mesh per coarse by default
    for key in bounds.keys():
        if key == 'x':
            block01['IM'] = len(bounds[key]) - 1
            block01['IT'] = block01['IM']
        elif key == 'y':
            block01['JM'] = len(bounds[key]) - 1
            block01['JT'] = block01['JM']
        elif key == 'z':
            block01['KM'] = len(bounds[key]) - 1
            block01['KT'] = block01['KM']
    
    # Optional Input IQUAD
    block01['IQUAD'] = 1 # default
    
    return block01


def _block02(bounds, voxel_zone):
    block02 = {}
    
    # fine intervals are 1 by default
    for key in bounds.keys():
        if key == 'x':
            block02['XMESH'] = bounds[key]
            block02['XINTS'] = 1
        elif key == 'y':
            block02['YMESH'] = bounds[key]
            block02['YINTS'] = 1
        elif key == 'z':
            block02['XZMESH'] = bounds[key]
            block02['ZINTS'] = 1  
    
    # zones list 
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
    block02['ZONES'] = np.zeros(shape=(im, jm*km), dtype=int)
    for i in range(im):
        for jk in range(jm*km):
            block02['ZONES'][i,jk] = voxel_zone[n]
            n += 1
            
    return block02
    

def _block03(xs_names):
    block03 = {}
    
    block03['LIB'] = 'BXSLIB' # default
    block03['NAMES'] = xs_names
    #block03['MAXORD'] = 
    
    # Figure out the rest of this block later
    
    return block03


def _block04(mat_lib, xs_names, nuc_names, zones):
    block04 = {}
    
    # replace nucids with bsxlib names
    mat_bxslib = {}
    for mat in mat_lib.keys():
        mat_bxslib[mat] = {}
        for nucid in mat_lib[mat].keys():
            if nucid in nuc_names.keys():
                name = nuc_names[nucid]
                mat_bxslib[mat][name] = mat_lib[mat][nucid]
            else:
                warn("Nucid {0} does not exist in nuc_names dictionary.".format(nucid))
                mat_bxslib[mat]["{0}".format(nucid)] = mat_lib[mat][nucid]
    
    block04['MATLS'] = mat_bxslib
    block04['ASSIGN'] = zones
    
    
    
    return block04


def _block05(ngroup, bounds, nmq):
    block05 = {}
    # need volumetric source def here
    # Calculation/Solver inputs
    
    block05['IEVT'] = 0 # default? 0 = source

    # create source definition
    # Volumetric source option 3
    for i in bounds.keys():
        if i == 'x':
            it = len(bounds[i]) - 1
            block05['SOURCX'] = np.zeros(shape=(it, nmq), dtype=float)
            block05['SOURCX'][:,0] = 1.0
        elif i == 'y':
            jt = len(bounds[i]) - 1
            block05['SOURCY'] = np.zeros(shape=(jt, nmq), dtype=float)
            block05['SOURCY'][:,0] = 1.0
        elif i == 'z':
            kt = len(bounds[i]) - 1
            block05['SOURCZ'] = np.zeros(shape=(kt, nmq), dtype=float)
            block05['SOURCZ'][:,0] = 1.0
    
    block05['SOURCE'] = np.zeros(shape=(ngroup, nmq), dtype=float)
    block05['SOURCE'][:,0] = 1.0
    
    return block05
    

def _write(title, block01, block02, block03, block04, block05):
    pass
    
    
