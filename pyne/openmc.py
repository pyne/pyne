"""Utilities for handling OpenMC.
"""
from __future__ import print_function
import os
import io
import sys
from warnings import warn
from collections import namedtuple
import numpy as np
import tables as tb

if sys.version_info[0] == 2:
    from HTMLParser import HTMLParser
else:
    from html.parser import HTMLParser

from pyne import nucname
from pyne.utils import QAWarning
from pyne.mesh import HAVE_PYMOAB
try:
    from pymoab import core as mb_core, hcoord, scd, types
    from pymoab.rng import subtract
    from pymoab.tag import Tag
    from pymoab.types import _eh_py_type, _TAG_TYPE_STRS
    HAVE_PYMOAB = True
except ImportError:
    HAVE_PYMOAB = False
    warn("The PyMOAB optional dependency could not be imported. "
         "Some aspects of the openmc module may be incomplete.", QAWarning)


from pyne.mesh import Mesh

warn(__name__ + " is not yet QA compliant.", QAWarning)

if sys.version_info[0] > 2:
    basestring = str


class AceTable(namedtuple('_AceTable', ['alias', 'awr', 'location', 'metastable',
                                        'name', 'path', 'temperature', 'zaid'])):
    """A simple data structure reprsenenting an <ace_table /> tag in a
    cross_sections.xml file.
    """
    def __new__(cls, alias=None, awr=None, location=None, metastable=None,
                name=None, path=None, temperature=None, zaid=None,
                cross_sections_path=None):
        return super(AceTable, cls).__new__(cls, alias=alias, awr=awr,
                                            location=location,
                                            metastable=metastable, name=name,
                                            path=path, temperature=temperature,
                                            zaid=zaid)

    def __init__(self, alias=None, awr=None, location=None, metastable=None,
                 name=None, path=None, temperature=None, zaid=None,
                 cross_sections_path=None):
        """Parameters
        ----------
        alias : str, optional
            ace_table attribute.
        awr : str, optional
            ace_table attribute.
        location : str, optional
            ace_table attribute.
        metastable : str, optional
            ace_table attribute.
        name : str, optional
            ace_table attribute.
        path : str, optional
            ace_table attribute.
        temperature : str, optional
            ace_table attribute.
        zaid : str, optional
            ace_table attribute. If set or non-zero then the nucid attribute
            will be set.
        cross_sections_path : str, optional
            If this and path are both present then the abspath attribute will be
            set.
        """
        super(AceTable, self).__init__()
        nuc = None
        if zaid is not None or zaid != '0':
            meta = "0" if metastable is None else metastable
            nuc = nucname.zzaaam_to_id(zaid + meta)
            if nuc == 0:
                pass
            elif not nucname.isnuclide(nuc):  # then it's in MCNP metastable form
                nuc = nucname.mcnp_to_id(zaid)
        self.nucid = nuc
        abspath = None
        if path is not None and cross_sections_path is not None:
            if os.path.isdir(cross_sections_path):
                d = cross_sections_path
            else:
                d = os.path.dirname(cross_sections_path)
            abspath = os.path.abspath(os.path.join(d, path))
        self.abspath = abspath

    def xml(self):
        """Creates an XML representation of the ACE Table.
        """
        s = '<ace_table '
        s += " ".join(['{0}="{1}"'.format(f, getattr(self, f)) for f in self._fields
                       if getattr(self, f) is not None])
        s += '/>'
        return s


class CrossSections(HTMLParser):
    """This class represents an OpenMC cross_sections.xml file.
    """

    def __init__(self, f=None):
        """Parameters
        ----------
        f : str, file-like, or None, optional
            This is a path to the cross_sections.xml file, a file handle, or
            None indicating an empty container.
        """
        # HTMLParser is only a new-style class in python 3
        if sys.version_info[0] > 2:
            super(CrossSections, self).__init__()
        else:
            HTMLParser.__init__(self)
        self.reset()
        self._tag = None
        self.path = None
        self.filetype = 'ascii'
        self.ace_tables = []
        if f is None:
            return
        opened_here = False
        if isinstance(f, str):
            opened_here = True
            self.path = f
            f = io.open(f, 'r')
        raw = f.read()
        self.feed(raw)
        if opened_here:
            f.close()

    def handle_starttag(self, tag, attrs):
        self._tag = tag
        if tag == 'ace_table':
            self.handle_ace_table(attrs)

    def handle_endtag(self, tag):
        self._tag = None

    def handle_startendtag(self, tag, attrs):
        if tag == 'ace_table':
            self.handle_ace_table(attrs)

    def handle_data(self, data):
        if self._tag == 'filetype':
            self.filetype = data
        elif self._tag == 'directory':
            self.path = data.strip()

    def handle_ace_table(self, attrs):
        ace_table = AceTable(cross_sections_path=self.path, **dict(attrs))
        self.ace_tables.append(ace_table)

    def xml(self):
        """Returns an XML representation of the cross sections file.
        """
        template = ('<?xml version="1.0" ?>\n'
                    '<cross_sections>\n'
                    '  <filetype>{filetype}</filetype>\n'
                    '  {ace_tables}\n'
                    '</cross_sections>\n')
        ace_tables = "\n  ".join([a.xml() for a in self.ace_tables])
        s = template.format(filetype=self.filetype, ace_tables=ace_tables)
        return s

def get_tally_results_from_openmc_sp(filename, tally_num):
    """
    This function reads a OpenMC state point file to get the results data for
    specific tally number.
    """
    # check tally_num exist
    tally_name = ''.join(["tally ", str(tally_num)])
    with tb.open_file(filename) as h5f:
        try:
            tally_results = h5f.root.tallies._f_get_child(
                    tally_name)._f_get_child('results')[:]
        except:
            raise ValueError("Tally {0} not found in file: {1}".format(
                str(tally_num), filename))
    return tally_results

def mesh_from_openmc_statepoint(filename, tally_num, particle_type='n'):
    """
    This function creates a Mesh instance from OpenMC statepoint file.

    Parameters:
    -----------
    filename : str
        Filename of the OpenMC statepoint file. It ends with ".h5",
        eg: "statepoint.10.h5".
    tally_num : int
        Tally number of specific mesh tally.
    particle_type : str
        Type of the tallied particle.

    Returns:
    --------
    mesh : Mesh
        PyNE Mesh instance.
    """
    # check tally_num exist
    tally_name = ''.join(["tally ", str(tally_num)])
    with tb.open_file(filename) as h5f:
        try:
            tally_results = get_tally_results_from_openmc_sp(filename,
                    tally_num)
            meshes = h5f.root.tallies._f_get_child('meshes')
            if meshes._v_nchildren != 1:
                raise ValueError("Only one mesh is support for each Tally now")
            mesh_str = meshes._v_groups.__str__()
            mesh_name = _get_mesh_name(mesh_str)
            mesh = meshes._f_get_child(mesh_name)
            structured_coords = calc_structured_coords(
                    mesh.lower_left[:],
                    mesh.upper_right[:],
                    mesh.dimension[:])
        except:
            raise ValueError("Tally {0} not found in {1}".format(str(tally_num), filename))

    # parameters to create mesh
    mesh = Mesh(mesh=None, structured=True, structured_coords=structured_coords)
    mesh.tag_flux_error_from_openmc_tally_results(tally_results)
    return mesh

def calc_structured_coords(lower_left, upper_right, dimension):
    """
    This function calculate the structured mesh coordinations from OpenMC mesh
    parameters.
    x_bounds, y_bounds, z_bounds = [], [], []

    Parameters:
    -----------
    lower_left : numpy array of float
        The lower left coordinate of the mesh. A numpy array of lenght 3.
        Format: [x_min, y_min, z_min].
    upper_right : numpy array of float
        The upper right coordinate of the mesh. A numpy array of lenght 3.
        Format: [x_max, y_max, z_max].
    dimension : numpy array of int
        Number of mesh intervals in each dimension. A numpy array of length 3.
        Format: [x_ints, y_ints, z_ints].

    Returns:
    --------
    structured_coords : numpy array
        A nested numpy array definning the boundaries of the mesh element in
        each dimension. Format: [[x_bounds1, x_bounds2. ...],
                                 [y_bounds1, y_bounds2, ...],
                                 [z_bounds1, z_bounds2, ...]]
    """
    # check the length of parameters
    if len(lower_left) != 3 or len(upper_right) != 3 or len(dimension) != 3:
        raise ValueError("Only 3D OpenMC mesh is supported!")
    structured_coords = []
    for dim in range(3):
        bounds = []
        step = (upper_right[dim] - lower_left[dim]) / dimension[dim]
        for i in range(dimension[dim] + 1):
            bounds.append(lower_left[dim] + i * step)
        structured_coords.append(bounds)
    return structured_coords


def _get_mesh_name(mesh_str):
    """
    This function is used to get mesh name from a string contain it.
    A mesh string contains the content such as:
    "{'mesh 14': /tallies/meshes/mesh 14 (Group)"

    Parameters:
    -----------
    mesh_str : str
        A mesh string contains mesh name.
    
    Returns:
    --------
    mesh_name : str
        The mesh name, Eg: "mesh 14"
    """
    ls = mesh_str.strip().split(':')
    mesh_name = ls[0].split("'")[1]
    return mesh_name
