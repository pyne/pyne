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
    coord_sys, bounds = _read_mesh(mesh)
    
    # Read the materials from the hdf5 and convert to correct naming convention
    mat_lib = _get_materials(hdf5, datapath, nucpath, nuc_names)
    #print(mat_lib)
    #mat_lib = Material.collapse_elements(mat_lib_expanded)
    
    # determine the zones
    zones = _define_zones(mesh)

    # read nucdata
    bxslib = open(nucdata, 'rb')
    xs_names = _read_bxslib(bxslib)

    return coord_sys, bounds, mat_lib, zones, xs_names
 
    
def _read_mesh(mesh):
    # determines the system geometry (1-D, 2-D, or 3-D Cartesian)
    # currently cartesian is only supported
    nx = len(mesh.structured_get_divisions("x"))
    ny = len(mesh.structured_get_divisions("y"))
    nz = len(mesh.structured_get_divisions("z"))
    
    # Check for dimensions with >1 voxel (>2 bounds)
    # This determines 1-D, 2-D, or 3-D
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

 
def _read_bxslib(bxslib):
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
    # cell # -> material name & vol fract -> isotope name & dens
    
    # set of exception nuclides for collapse_elements
    mat_exceptions = Set(nuc_names.keys())
    
    # collapse isotopes into elements
    mats = MaterialLibrary(hdf5,datapath=datapath,nucpath=nucpath)
    mats_collapsed = {}
    for mat_name in mats.keys():
        mats_collapsed[mat_name] = mats[mat_name].collapse_elements(mat_exceptions)
    
    # Check that the materials are valid:
    #   1) non zero and non-negative densities (density = True)
    #   2) set of nuclides is not empty (else it is vacuum) (empty = False)
    #   3) nucids appear in nuc_names
    # might put 2 and 3 later      
    
    # convert mass fraction to atom fraction and then to [at/b-cm]
    Na = 6.022*(10.**23) # Avagadro's number [at/mol]
    barn_conv = 10.**-24 # [cm^2/b]
    mat_lib = {}
    for name, comp in mats_collapsed.iteritems():
        print(comp)
        mat_name = '{0}'.format(name.split(':')[1])
        comp_atom_frac = comp.to_atom_frac() # atom fractions
        #print(comp_atom_frac)
        
        density = comp.mass_density() # [g/cc]
        if density < 0.0:
            warn("Material {0} has an invalid negative density.".format(mat_name))
        
        mol_mass = comp.molecular_mass() # [g/mol]
        comp_list = {}
        total = 0
        for nucid, frac in comp_atom_frac.iteritems():
            comp_list[nucid] = frac*density*Na*barn_conv/mol_mass # [at/b-cm]
            total += frac*density*Na*barn_conv/mol_mass # [at/b-cm]
        mat_lib[mat_name] = comp_list
        #print(total)    
    
    #print(mat_lib)
        
    return mat_lib


def _define_zones(mesh):
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
    
    # determine which voxels are identical and remove
    z = 0
    zones = {} # defined by cell number
    match = False
    first = True    
    for idx, vals in voxel.iteritems():
        for zone, info in zones.iteritems():
            if vals == info:
                match = True
                break
            else:
                match = False
        if first or not match:
            z += 1
            zones[z] = voxel[idx]
            first = False
            
    
    ## Add section here which replaces cell number with materials
            
    return zones


#class PartisnWrite(object):

def write_partisn_input(coord_sys, bounds, mat_lib, zones, xs_names, input_file):
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
    _block01(coord_sys, xs_names, mat_lib, zones, bounds)

def _title():
    # figure out what to make the title
    pass
        
def _block01(coord_sys, xs_names, mat_lib, zones, bounds):
    # Determine IGEOM
    if len(coord_sys) == 1:
        IGEOM = 'SLAB'
    elif len(coord_sys) == 2:
        IGEOM = 'X-Y' # assuming cartesian
    elif len(coord_sys) == 3:
        IGEOM = 'X-Y-Z' # assuming cartesian
    
    # NGROUP
    
    # ISN
    
    NISO = len(xs_names)
    MT = len(mat_lib)
    NZONE = len(zones)
    
    # Number of Fine and Coarse Meshes
    # one fine mesh per coarse by default
    for key in bounds.keys():
        if key == 'x':
            IM = len(bounds[key]) - 1
            IT = IM
        elif key == 'y':
            JM = len(bounds[key]) - 1
            JT = IM
        elif key == 'z':
            KM = len(bounds[key]) - 1
            KT = IM
    
    # Optional Input IQUAD
    IQUAD = 1 # default

def _block02():
    pass

def _block03():
    pass

def _block04():
    pass

def _write():
    pass
    
    
