#!/usr/bin/env python
""" Module for the production of PartiSn input decks. PartiSn is a discrete
ordinates code produced by Los Almos National Laboratory (LANL). Can be used
to produce neutron, photon, or coupled neutron photon prblems, adjoint or
forward or time dependent problems can be run.

The Input class is the on being worked on currently and should need the least work

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

class PartisnInput():

    def _block1(mesh):
        """This function reads a structured  mesh object and returns the 1st data
        block of input, the 1st data block contains the definition of the 
        problem initialisation, specifically, the number of energy groups, the
        Sn order 
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block1_str : str
        A string containing the full block 1 definition for a PartiSn input
        """
        
        block1_str = "/A# block 1 \n"
        block1_str += "t \n"
        return block1_str


    def _block2(mesh):
        """This function reads a structured  mesh object and returns the 2nd data
        block of input, the 2nd data block contains the full definition for the 
        problem geometry and material assignments
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block2_str : str
        A string containing the full block 2 definition for a PartiSn input
        """
        
        block2_str = "/A# block 2 \n"
        block2_str += PartisnInput._partisn_geom(mesh)
        block2_str += PartisnInput._partisn_material(mesh)
        block2_str += "t \n"

        return block2_str

    def _block3(mesh):
        """This function reads a structured  mesh object and returns the 3rd data
        block of input, the 3rd data block specifies the cross section library
        and edits that you would like. This will be a wall of text that is valid
        PartiSn input
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block3_str : str
        A string containing the full block 3 definition for a PartiSn input
        """
        
        block3_str = "/A# block 3 \n"
        block3_str += "t \n"
        return block3_str

    def _block4(mesh):
        """This function reads a structured  mesh object and returns the 4th data
        block of input, the 4th data block specifies the material definitions
        from the problem. The "pure" materials which make up the problem are defined
        and then the mixtures that make up the zones are required, i.e. each unique
        mixture needs to be represented
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block4_str : str
        A string containing the full block 4 definition for a PartiSn input
        """
        
        block4_str = "/A# block 4 \n"
        block4_str += "t \n"
        return block4_str


    def _block5(mesh):
        """This function reads a structured  mesh object and returns the 5th data
        block of input, the 5th data block specifies the source behaviour and 
        normalisation
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        block5_str : str
        A string containing the full block 5 definition for a PartiSn input
        """
        
        block5_str = "/A# block 5 \n"
        block5_str += "t \n"
        return block5_str
    
    def _partisn_geom(mesh):
        """This function reads a structured  mesh object and returns the mesh 
        portion of a PartiSn input deck
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        geom : str
        A string containing the PartiSn mesh boundaries and geometry that 
        can be written directly as valid syntax 
        """
        geom = "* temporary placeholder for geom \n"
        return geom

    def _partisn_material(mesh):
        """This function reads a structured  mesh object and returns the material
        assignnent portion of a PartiSn input deck
        
        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        material : str
        A string containing the PartiSn material assignments and is valid
        PartiSn syntax
        """
        material = "* temporary placeholder for material \n"
        
        return material

    def write_partisn(mesh, filename):
        """This function reads a structured mesh object and returns the 
        the complete PartiSn input deck and writes it out to file

        Parameters
        ----------
        mesh : PyNE Mesh object
        
        Returns
        -------
        filename : str
        A file to write the output to, note this will always overwrite the
        file if it exists
        """
        
        output_data = ("* PartiSn input deck produced automatially by PyNE \n"
                       "* pyne.partisn module \n")
        output_data += PartisnInput._block1(mesh)
        output_data += PartisnInput._block2(mesh)
        output_data += PartisnInput._block3(mesh)
        output_data += PartisnInput._block4(mesh)
        output_data += PartisnInput._block5(mesh)
        output_data += ("/ ********************* \n"
                        "* You must produce your own edits\n")

        f = open(filename,'w')
        f.write(output_data) # python will convert \n to os.linesep
        f.close()
