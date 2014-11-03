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
from pyne.utils import VnVWarning
import itertools

import numpy as np
import tables

from pyne import dagmc
from pyne.material import Material
from pyne.material import MultiMaterial
from pyne.material import MaterialLibrary

from pyne import nucname
from pyne.binaryreader import _BinaryReader, _FortranRecord

warn(__name__ + " is not yet V&V compliant.", VnVWarning)

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
    mesh :: mesh object, a premade mesh object that conforms to the 
            geometry.
                Can be 1-D, 2-D, or 3-D.
                note: only Cartesian based geometries are supported.
    hdf5 :: file, a material-laden dagmc geometry file.
    nucdata :: file, nuclear data cross section library.
                note: only BXSLIB format is currently supported.
    nuc_names :: dict, pyne element/isotope names to bxslib name assignment,
                keys are pyne nucids (int) and values are bxslib names (str)
    datapath :: str, optional, The path in the heirarchy to the data table 
            in an HDF5 file. (for MaterialLibrary)
                default = material_library/materials
    nucpath :: str, optional, The path in the heirarchy to the 
            nuclide array in an HDF5 file. (for MaterialLibrary)
                default = material_library/nucid
    ****fine :: dict of lists, optional, number of fine mesh intervals per coarse
            mesh interval. Fine mesh is used for solver.
                keys must be 'x', 'y', or 'z'. Value is list of ints. List 
                can be of length 1 so that the value is applied to all coarse
                meshes in that direction. Or list can be length of coarse
                meshes in that direction.
                    example: fine = {'x':[5], 'y':[3,5,6,8,2]}
                default: 10 in all directions.
                NOTE: make user supply fine mesh only, make it default one fine
                mesh per coarse mesh
Attributes
----------
    dim :: list of str, specifies the dimensions in problem. Currently
            only Cartesian is supported so can be any combination of one or
            more of 'x', 'y', or 'z'.
    bounds :: dict of lists of floats, values for the coarse mesh bounds
            in each dimension present
    matlib :: dict of dict, keys are names of pyne materials whose keys are
            *pyne* element/isotope names and their value is *density*
    xs_names :: list of strings, names of isotope/elements from the bxslib
    fine :: dict of lists, number of fine mesh intervals per coarse mesh in 
            bounds, keys are the dimensions (x, y, or z), values is list of
            length of bounds in each direction.        
"""

def read_hdf5_mesh(mesh, hdf5, nucdata, nuc_names, **kwargs):
    dagmc.load(hdf5)
    # optional inputs
    datapath = kwargs['datapath'] if 'datapath' in kwargs else 'material_library/materials'
    nucpath = kwargs['nucpath'] if 'nucpath' in kwargs else 'material_library/nucid'

    # get coordinate system and mesh bounds from mesh       
    coord_sys, bounds = _read_mesh(mesh)

    # Read the materials from the hdf5 and convert to correct naming convention
    mat_lib = _get_materials(hdf5, datapath, nucpath, nuc_names)       
    
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
    for i in coord_sys:
        bounds[i] = mesh.structured_get_divisions(i)
    
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
    
    mats = MaterialLibrary(hdf5,datapath=datapath,nucpath=nucpath)
    mat_lib = {}
    
    for key in mats.keys():
        mat_name = '{0}'.format(key.split(':')[1])
        mat_lib[mat_name] = {}
        for nuc in mats[key]:
            if nuc in nuc_names.keys():
                nuc_name = nuc_names[nuc]
            else:
                nuc_name = nuc
                warn("PyNE nuclide {0} not listed in nuc_names".format(nuc))
            # convert mass fraction to atom density [at/b*cm]
            # rho/mw*Na    
            rho_frac = mats[key][nuc_name]
            mw = mats[key].molecular_mass()
            Na = 60.22
            mat_lib[mat_name][nuc_name] = rho_frac/mw*Na
            
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
    zones = {}
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
            
    return zones


#class PartisnWrite(object):
"""This class writes out the information stored by PartisnRead to
a text partisn input file.
"""
def write_partisn(coord_sys, bounds, mat_lib, zones, xs_names, input_file):
    block01(coord_sys, xs_names, mat_lib, zones, bounds)

def title():
    # figure out what to make the title
    pass
        
def block01(coord_sys, xs_names, mat_lib, zones, bounds):
    
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

def block02():
    pass

def block03():
    pass

def block04():
    pass
    
    
