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
    """ This class read in a material-laden mesh object and an h5m 
    geometry file and stores each attribute. It then manipulates/
    rearranges and stores the information into the necessary 
    Partisn blocks.
    """
    def __init__(self, mesh_obj, h5m_file):
        
        h5mf = open(h5m_file)
        discretized_mesh = dagmc.discretize_geom(mesh_obj)
        

class PartisnWrite(object):
    """This class writes out the information stored by PartisnRead to
    a text partisn input file.
    """
    def __init__(self):
        pass
        
