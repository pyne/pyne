"""Utilities for handling OpenMC.
"""
from __future__ import print_function
import os
import io
import sys
from warnings import warn

try:
    from collections.abc import namedtuple
except ImportError:
    from collections import namedtuple
import numpy as np
import tables as tb

if sys.version_info[0] == 2:
    from HTMLParser import HTMLParser
else:
    from html.parser import HTMLParser

from pyne.mesh import MeshTally, HAVE_PYMOAB

if HAVE_PYMOAB:
    from pyne.mesh import NativeMeshTag
else:
    warn(
        "The PyMOAB optional dependency could not be imported. "
        "Some aspects of the openmc module may be incomplete.",
        ImportWarning,
    )
from pyne import nucname
from pyne.utils import QA_warn

QA_warn(__name__)

try:
    import openmc
except:
    warn(
        "The openmc (OpenMC Python API) could not be imported. "
        "Some aspects of the openmc module may be incomplete."
    )

if sys.version_info[0] > 2:
    basestring = str


class AceTable(
    namedtuple(
        "_AceTable",
        [
            "alias",
            "awr",
            "location",
            "metastable",
            "name",
            "path",
            "temperature",
            "zaid",
        ],
    )
):
    """A simple data structure reprsenenting an <ace_table /> tag in a
    cross_sections.xml file.
    """

    def __new__(
        cls,
        alias=None,
        awr=None,
        location=None,
        metastable=None,
        name=None,
        path=None,
        temperature=None,
        zaid=None,
        cross_sections_path=None,
    ):
        return super(AceTable, cls).__new__(
            cls,
            alias=alias,
            awr=awr,
            location=location,
            metastable=metastable,
            name=name,
            path=path,
            temperature=temperature,
            zaid=zaid,
        )

    def __init__(
        self,
        alias=None,
        awr=None,
        location=None,
        metastable=None,
        name=None,
        path=None,
        temperature=None,
        zaid=None,
        cross_sections_path=None,
    ):
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
        if zaid is not None or zaid != "0":
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
        """Creates an XML representation of the ACE Table."""
        s = "<ace_table "
        s += " ".join(
            [
                '{0}="{1}"'.format(f, getattr(self, f))
                for f in self._fields
                if getattr(self, f) is not None
            ]
        )
        s += "/>"
        return s


class CrossSections(HTMLParser):
    """This class represents an OpenMC cross_sections.xml file."""

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
        self.filetype = "ascii"
        self.ace_tables = []
        if f is None:
            return
        opened_here = False
        if isinstance(f, str):
            opened_here = True
            self.path = f
            f = io.open(f, "r")
        raw = f.read()
        self.feed(raw)
        if opened_here:
            f.close()

    def handle_starttag(self, tag, attrs):
        self._tag = tag
        if tag == "ace_table":
            self.handle_ace_table(attrs)

    def handle_endtag(self, tag):
        self._tag = None

    def handle_startendtag(self, tag, attrs):
        if tag == "ace_table":
            self.handle_ace_table(attrs)

    def handle_data(self, data):
        if self._tag == "filetype":
            self.filetype = data
        elif self._tag == "directory":
            self.path = data.strip()

    def handle_ace_table(self, attrs):
        ace_table = AceTable(cross_sections_path=self.path, **dict(attrs))
        self.ace_tables.append(ace_table)

    def xml(self):
        """Returns an XML representation of the cross sections file."""
        template = (
            '<?xml version="1.0" ?>\n'
            "<cross_sections>\n"
            "  <filetype>{filetype}</filetype>\n"
            "  {ace_tables}\n"
            "</cross_sections>\n"
        )
        ace_tables = "\n  ".join([a.xml() for a in self.ace_tables])
        s = template.format(filetype=self.filetype, ace_tables=ace_tables)
        return s


def get_e_bounds_from_openmc_sp(filename, tally_id):
    """
    This function reads OpenMC state point file to get the energy boundaries
    for a specific tally number.

    Parameters:
    -----------
    filename : str
        The OpenMC state point file name.
    tally_id : int
        Tally id to read.

    Returns:
    --------
    e_bounds : numpy array
        Energy boundries with size of (num_e_groups + 1).
    """

    sp = openmc.StatePoint(filename)
    tally = sp.get_tally(id=tally_id)
    # check filter to find EnergyFilter
    for flt in tally.filters:
        if isinstance(flt, openmc.filter.EnergyFilter):
            energy_filter = flt
    e_bounds = energy_filter.values
    return e_bounds


def get_structured_coords_from_openmc_sp(filename, tally_id):
    """
    This function reads the OpenMC state point file and gets the structured
    coordinates of the mesh.

    Parameters:
    -----------
    filename : str
        OpenMC state point filename.
    tally_id : int
        Tally id.

    Returns:
    --------
    structured_coords : numpy array
        A nested numpy array definning the boundaries of the mesh element in
        each dimension. Format: [[x_bounds1, x_bounds2. ...],
                                 [y_bounds1, y_bounds2, ...],
                                 [z_bounds1, z_bounds2, ...]]
    """

    sp = openmc.StatePoint(filename)
    tally = sp.get_tally(id=tally_id)
    # check filters to find MeshFilter
    for flt in tally.filters:
        if isinstance(flt, openmc.filter.MeshFilter):
            if isinstance(flt.mesh, openmc.RegularMesh):
                mesh_filter = flt
            else:
                raise NotImplementedError(
                    "The mesh on this tally: {0} is not "
                    "yet supported. Please use an openmc.RegularMesh".format(tally_id)
                )
    structured_coords = calc_structured_coords(
        mesh_filter.mesh.lower_left[:],
        mesh_filter.mesh.upper_right[:],
        mesh_filter.mesh.dimension[:],
    )
    return structured_coords


def calc_structured_coords(lower_left, upper_right, dimension):
    """
    This function calculates the structured mesh coordinates from an OpenMC mesh
    parameters.

    Parameters:
    -----------
    lower_left : numpy array of float
        The lower left coordinate of the mesh. A numpy array of length 3.
        Format: [x_min, y_min, z_min].
    upper_right : numpy array of float
        The upper right coordinate of the mesh. A numpy array of length 3.
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
        raise NotImplementedError("Only 3D OpenMC mesh is supported!")
    structured_coords = []
    for dim in range(3):
        bounds = []
        step = (upper_right[dim] - lower_left[dim]) / dimension[dim]
        for i in range(dimension[dim] + 1):
            bounds.append(lower_left[dim] + i * step)
        structured_coords.append(bounds)
    return structured_coords


def get_result_error_from_openmc_sp(filename, m):
    """
    Convert the openmc flux into result, rel_err, res_tot, rel_err_tot.

    Parameters:
    -----------
    filename : str
        The OpenMC state point file name.
    m : MeshTally
        MeshTally for the tally.

    Returns:
    --------
    result : numpy array
        This numpy array contains the flux data read from the OpenMC state point
        file. The shape of this numpy array is (num_ves*num_e_groups), with z
        changes fastest.
    rel_error: numpy array
        This numpy array contains the relative error data.
    res_tot : numpy array
        The total results.
    rel_err_tot : numpy array
        Relative error of total results.
    """

    sp = openmc.StatePoint(filename)
    tally = sp.get_tally(scores=["flux"], id=m.tally_number)
    flux = tally.get_slice(scores=["flux"])

    num_ves = len(m)
    # assumption: openmc mesh voxels volumes are uniform for a regular mesh
    ve_vol = m.structured_hex_volume(0, 0, 0)
    num_e_groups = len(flux.mean.flatten()) // num_ves

    # get result and res_tot
    result = flux.mean.flatten()
    result = np.reshape(result, newshape=(num_e_groups, num_ves))
    result = result.transpose()
    res_tot = np.sum(result, axis=1)

    # calculate rel_err
    rel_err = np.zeros_like(flux.std_dev)
    nonzero = flux.mean > 0
    rel_err[nonzero] = flux.std_dev[nonzero] / flux.mean[nonzero]
    rel_err = np.reshape(rel_err.flatten(), newshape=(num_e_groups, num_ves))
    rel_err = rel_err.transpose()
    # calculate rel_err_tot, assuming independent in each energy bin (not true)
    # due to the lack of covariance, this could be smaller than the true value
    rel_err_tot = np.zeros_like(res_tot)
    std_dev = np.reshape(flux.std_dev.flatten(), newshape=(num_e_groups, num_ves))
    std_dev = std_dev.transpose()
    var_tot = np.sum(np.square(std_dev), axis=1)
    nonzero = res_tot > 0
    rel_err_tot[nonzero] = np.sqrt(var_tot[nonzero]) / res_tot[nonzero]

    # get volume averged results
    result = np.divide(result, ve_vol)
    res_tot = np.divide(res_tot, ve_vol)
    # In the result of OpenMC, x changes fastest
    # change the order to z changes fastest
    result = result_changes_order(result, m.dims[3:6])
    rel_err = result_changes_order(rel_err, m.dims[3:6])
    res_tot = result_changes_order(res_tot, m.dims[3:6])
    rel_err_tot = result_changes_order(rel_err_tot, m.dims[3:6])
    return result, rel_err, res_tot, rel_err_tot


def result_changes_order(result, dims):
    """
    This function changes the order of openmc flux data. The default order of
    OpenMC flux is x changes the fastest. This function reorganizes the data
    such that z changes fastest (the same order as the MCNP meshtal file).
    which z changes the fasest (the same order as MCNP meshtal).

    Parameters:
    -----------
    result : numpy array
        This numpy array contains the flux data read from OpenMC state point
        file. The shape of this numpy array is (num_ves, num_e_groups) or
        (num_ves, ), with x changing the fastest.
    dims : list
        The dimension of x, y, z.

    Returns:
    --------
    result_ : numpy array
        This numpy array contains the flux data read from OpenMC state point
        file. The shape of this numpy array is (num_ves, num_e_groups) or
        (num_ves, ), with z chaning the fastest.
    """

    x_dim, y_dim, z_dim = dims[0], dims[1], dims[2]
    num_e_groups = result.size // (x_dim * y_dim * z_dim)
    # the original result, energy changes the fastest, x next, y next, then z
    # reshape result to 4 dimension (z, y, x, e, with e changes fastest) array
    result_ = np.reshape(result, newshape=(z_dim, y_dim, x_dim, num_e_groups))

    # the target sequence, energy changes fastest, z next, y next, then x
    # move the axis from [0, 1, 2, 3] (z, y, x, e) to
    #                    [2, 1, 0, 3] (x, y, z, e)
    result_ = np.moveaxis(result_, [0, 1, 2, 3], [2, 1, 0, 3])

    # reshape the 4 dimension (x, y, z, e, with e changes fastest) array to
    #             2 dimension (x_dim*y_dim*z_dim, num_e_groups) array
    result_ = np.reshape(result_, newshape=(x_dim * y_dim * z_dim, num_e_groups))

    return result_


def create_meshtally(
    filename, tally_id, particle=None, tag_names=None, mesh_has_mats=False
):
    """
    This function creates a MeshTally instance from OpenMC state point file.

    Parameters:
    -----------
    filename : str
        Filename of the OpenMC statepoint file. It ends with ".h5",
        eg: "statepoint.10.h5".
    tally_id : int
        Tally number.
    particle : str
        The particle type, 'neutron' or 'photon'.
    tag_names : iterable, optional
        Four strs that specify the tag names for the results, relative
        errors, total results and relative errors of the total results.
    mesh_has_mats: bool
        If false, Meshtally objects will be created without PyNE material
        objects.

    Returns:
    --------
    m : MeshTally object
        The MeshTally object created from OpenMC state point file with tally
        number of tally_id.
    """

    m = MeshTally()
    # assign tally_number
    m.tally_number = tally_id
    # assign particle
    if particle != None:
        m.particle = particle
    # assign tag_names
    if tag_names is None:
        m.set_default_tag_names()
    else:
        m.tag_names = tag_names
    # parameters to create mesh
    structured_coords = get_structured_coords_from_openmc_sp(
        filename, tally_id=m.tally_number
    )
    m.x_bounds = tuple(structured_coords[0])
    m.y_bounds = tuple(structured_coords[1])
    m.z_bounds = tuple(structured_coords[2])
    m.dims = [len(m.x_bounds) - 1, len(m.y_bounds) - 1, len(m.z_bounds) - 1]
    m.num_ves = np.prod(m.dims)
    m.e_bounds = tuple(get_e_bounds_from_openmc_sp(filename, m.tally_number))
    m.num_e_groups = len(m.e_bounds) - 1
    if mesh_has_mats:
        raise NotImplementedError(
            "mesh_has_mats is currently not supported " "for OpenMC"
        )
    else:
        mats = None
    super(MeshTally, m).__init__(
        structured_coords=structured_coords, structured=True, mats=mats
    )
    result, rel_err, res_tot, rel_err_tot = get_result_error_from_openmc_sp(filename, m)
    m.tag_flux_error_from_tally_results(result, rel_err, res_tot, rel_err_tot)
    return m
