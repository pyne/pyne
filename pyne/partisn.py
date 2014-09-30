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

class PartisnReadMesh(object):
    """This class will create the pyne mesh based a provided h5m geometry
    and user-specified values for the mesh. These user specified values
    are contained within a supplementary file.
    """
    def __init__(self, h5m_f, sup_f):
        """Opens the two provided files and creates a mesh
        h5m_f :: path to h5m file
        sup_f :: path to supplementary file containing bound data
        """

        fs = open(sup_f) # supplementary file containing mesh bounds
        fg = open(h5m_f) # h5m file containing

        self._create_bounds(fs)
        
    def _create_bounds(self, fs):
        # Reads the file and stores the bound data in a list
        # Note: This will store cartesian data ONLY
        """Data format:
            - Lines can be in any order
            - Must provide at least 2 different dimensions
            - "ic" refers to the to the list of coarse mesh floating
                point values for the ith direction
            - "if" refers to the list integer number of fine mesh intervals
                in each coarse mesh for the ith direction
            - the number of if entries must be one less than the number of
                ic entries for the ith direction
            - the jth entry in the fine mesh list refers to the number of 
                intervals between value j and j+1 in the coarse mesh list

                    xc = [a,b,c]
                    yc = [a,b,c]
                    zc = [a,b,c]
                    xf = [d,e]
                    yf = [d,e]
                    zf = [d,e]
        """
        x_tf = False
        y_tf = False
        z_tf = False
        line = fs.readline()
        coarse_list = {}
        fine_list = {}
        while line:
            if not line.isspace():
                if (line.strip()[0][0] != "#"):

                    coord = line.split("=")[0].strip()[0]
                    if (coord == 'x') or (coord == 'X'):
                        x_tf = True
                    elif (coord == 'y') or (coord == 'Y'):
                        y_tf = True
                    elif (coord == 'z') or (coord == 'Z'):
                        z_tf = True
                    else:
                        print("Coordinate system not supported")
    
                    if line.split("=")[0].strip()[1] == 'c':
                        # coarse mesh values
                        if coord not in coarse_list.keys():
                            coarse_list[coord] = []
                        coarse_list[coord] += [float(x) for x in line.split('[')[1].split(']')[0].split(',')]

                    elif line.split("=")[0].strip()[1] == 'f':
                        # fine mesh values
                        if coord not in fine_list.keys():
                            fine_list[coord] = []
                        fine_list[coord] += [int(x) for x in line.split('[')[1].split(']')[0].split(',')]

            line = fs.readline()
        print(coarse_list, fine_list)

h5m_f = "/userspace/k/kiesling/Documents/CNERG/partisn_geoms/fngn_checked_zip.h5m"
sup_f = "/userspace/k/kiesling/Documents/CNERG/partisn_geoms/test_bounds"
my_file = PartisnReadMesh(h5m_f, sup_f)
