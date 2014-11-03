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
    coord_sys :: int, optional, defines the coordinate system used.
                1 = Cartesian (default)
                2 = cylindrical
                3 = spherical
    datapath :: str, optional, The path in the heirarchy to the data table 
            in an HDF5 file. (for MaterialLibrary)
                default = material_library/materials
    nucpath :: str, optional, The path in the heirarchy to the 
            nuclide array in an HDF5 file. (for MaterialLibrary)
                default = material_library/nucid
        
    fine :: dict of lists, optional, number of fine mesh intervals per coarse
            mesh interval. Fine mesh is used for solver.
                keys must be 'x', 'y', or 'z'. Value is list of ints. List 
                can be of length 1 so that the value is applied to all coarse
                meshes in that direction. Or list can be length of coarse
                meshes in that direction.
                    example: fine = {'x':[5], 'y':[3,5,6,8,2]}
                default: 10 in all directions.
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
    
    coord_sys = kwargs['coord_sys'] if 'coord_sys' in kwargs else 1
    if coord_sys != 1:
        warn("Only Cartesian geometries are currently supported")
    
    datapath = kwargs['datapath'] if 'datapath' in kwargs else 'material_library/materials'
    nucpath = kwargs['nucpath'] if 'nucpath' in kwargs else 'material_library/nucid'

    dagmc_geom = dagmc.load(hdf5)
    dg = dagmc.discretize_geom(mesh)
           
    coord_sys, bounds = _read_mesh(mesh)

    #xs_names = _read_nucdata(nucdata)

    # Read the materials from the hdf5
    matlib = _get_materials(hdf5, datapath, nucpath, nuc_names)       
    
        
    
    # set fine mesh
    fine = {}
    for i in coord_sys:
        if 'fine' in kwargs:
            if len(kwargs['fine'][i]) == 1:
                fine[i] = kwargs['fine'][i]*len(bounds[i])
            else:
                if len(kwargs['fine'][i]) == len(bounds[i]):
                    fine[i] = kwargs['fine'][i]
                else:
                    warn("Number of fine mesh entries must be same length as coarse mesh entries. Using default value of 10.")
                    fine[i] = [10]*len(bounds[i])
        else:
            fine[i] = [10]*len(bounds[i])
    
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

def _read_nucdata(nucdata_file):
    # figure out how to generically determine the form of lib
    nucdata = open(nucdata_file,'rb')
    libname = 'bxslib'
    if libname == 'bxslib':
        return _read_bxslib(nucdata)
    # make it so that more lib types can be included
    
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
        # mat_name will be used for the 
        mat_name = '{0}'.format(key.split(':')[1])
        mat_lib[mat_name] = {}
        for nuc in mats[key]:
            # insert conversion of element name to bxslib name here and then
            # use that new name as the dict key instead of 'element'
            
            if nuc in nuc_names.keys():
                nuc_name = nuc_names[nuc]
            else:
                nuc_name = nuc
                warn("Pyne nuclide {0} not listed in nuc_names".format(nuc))
            # convert mass fraction to atom density [at/b*cm]
            # rho/mw*Na    
            rho_frac = mats[key][nuc_name]
            mw = mats[key].molecular_mass()
            Na = 60.22
            mat_lib[mat_name][nuc_name] = rho_frac/mw*Na
            
    return mat_lib
        
def _define_zones(dg):
    ### !!! NOT FINSIHED !!! ###
    # defines the "zones" based on unique discretize_geom results
    # dg = discretize_geom record array
    
    voxel = {}
    
    # define a single voxel
    for i in dg:
        idx = i[0]
        #print(i[1])
        #print(idx)
        if idx not in voxel.keys():
            voxel[idx] = {}
            voxel[idx]['cell'] = []
            voxel[idx]['vol_frac'] = []
            #voxel[idx]['rel_error'] = []
            
        voxel[idx]['cell'].append(i[1])
        voxel[idx]['vol_frac'].append(i[2])
        #voxel[idx]['rel_error'].append(i[3])
    
    # !!!! Come back to this later !!!!
    # determine which voxels are identical
    z = 1    # start zone counter
    zones = {}
    for idx in voxel.keys():
        if z not in zones.keys():
            zones[z] = {}
            zones[z]['cell'] = voxel[idx]['cell']
            zones[z]['vol_frac'] = voxel[idx]['vol_frac']
        else:
            for zz in zones.keys():
                c_tf = False
                vf_tf = False
                if zones[zz]['cell'] == voxel[idx]['cell']:
                    c_tf = True
                if zones[zz]['vol_frac'] == voxel[idx]['vol_frac']:
                    vf_tf = True
                if c_tf and vf_tf:
                    break
    print(zones)
                    
    
    #for item in voxel.iteritems():
    #    print(item)
        
    names = dg.dtype.names
    print(names)


#class PartisnWrite(object):
"""This class writes out the information stored by PartisnRead to
a text partisn input file.
"""
def write_partisn(self):
    pass
    
