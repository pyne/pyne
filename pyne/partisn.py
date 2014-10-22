#!/usr/bin/env python

""" Module for the production of PartiSn input decks. PartiSn is a discrete
ordinates code produced by Los Almos National Laboratory (LANL). Can be used
to produce neutron, photon, or coupled neutron photon prblems, adjoint or
forward or time dependent problems can be run.

The module is designed to operate on either 2D or 3D meshes, and produce the
appropriate input. It would be lovely if we eventually manage to get it working
with 1D as well as this appears to be a common mode of operation for PartiSn.

The Input class is the on being worked on currently and should need the least work
to improve. Fundamental inputs to the PartiSn class are:
    cell_fracs, a list of cell fractions with the number of materials
                as produced by DG
    mesh, a PyNE mesh instance including materials
    bxslib, the filename of the PartiSn cross section file

Next should be the Output class to read the output file and rtflux file

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


class PartisnRead(object):
    """ This class reads all necessary attributes from a material-laden
    geometry file, a pre-made PyNE mesh object, and the nuclear data
    cross section library.
    
    Parameters
    ----------
        mesh :: mesh object, a premade mesh object that conforms to the 
            geometry.
            Can be 1-D, 2-D, or 3-D.
            note: only Cartesian based geometries are supported.
        hfm :: file, a material-laden dagmc geometry file.
        nucdata :: file, nuclear data cross section library.
            note: only BXSLIB format is currently supported.
        coord_sys :: int, (optional), defines the coordinate system used.
            1 = Cartesian (default)
            2 = cylindrical
            3 = spherical
        fine :: int, (optional), number of fine mesh intervals per
            coarse mesh.
            default = 10
        
    Attributes
    ----------
        bounds :: dict of lists of floats, values for the mesh bounds
            in each dimension present
        fine :: int, (optional), number of fine mesh intervals
        

        
    """
    
    def __init__(self, mesh, h5m, nucdata, **kwargs):
        
        # read optional inputs
        coord_sys = kwargs['coord_sys'] if 'coord_sys' in kwargs else 1
        if coord_sys != 1:
            warn("Only Cartesian geometries are currently supported")
            
        #self.fine = kwargs['fine'] if 'fine' in kwargs else 10
        
        dagmc_geom = dagmc.load(h5m)
        dg = dagmc.discretize_geom(mesh)
               
        # determine if 1D, 2D, or 3D
        dim = self.get_dimensions(mesh)

        
        # collect the bounds data
        self.bounds = {}
        for i in dim:
            self.bounds[i] = mesh.structured_get_divisions(i)
        
        self._define_zones(dg)
        
        
    def get_dimensions(self, mesh):
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
        
        # Return dimension data
        if dim == 1:
            return [i]
        elif dim == 2:
            return [i, j]
        elif dim == 3:
            return [i, j, k]
            
    
    def _read_materials(self, dagmc_geom):
        # reads material properties from the loaded dagmc_geometry
        # cell # -> material name & vol fract -> isotope name & dens
        pass
        
    def _define_zones(self, dg):
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
        
       

class PartisnWrite(object):
    """This class writes out the information stored by PartisnRead to
    a text partisn input file.
    """
    def __init__(self):
        pass
        
