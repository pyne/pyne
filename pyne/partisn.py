#!/usr/bin/env python

""" Module for the production of PartiSn input decks. PartiSn is a discrete
ordinates code produced by Los Almos National Laboratory (LANL). Can be used
to produce neutron, photon, or coupled neutron photon prblems, adjoint or
forward or time dependent problems can be run.

Module is designed to work on 1D, 2D, or 3D Cartesian geometries.

If PyMOAB not installed then this module will not work.
"""

from __future__ import print_function, division
from pyne.mesh import HAVE_PYMOAB
import sys
import string
import struct
import math
import os
import linecache
import datetime
from textwrap import wrap
from warnings import warn
from pyne.utils import QA_warn
import itertools

import numpy as np
import tables

from pyne.material import Material
from pyne.material import MultiMaterial
from pyne.material_library import MaterialLibrary

from pyne import nucname
from pyne.binaryreader import _BinaryReader, _FortranRecord

QA_warn(__name__)

# Mesh specific imports

if HAVE_PYMOAB:
    from pyne.mesh import Mesh, StatMesh, MeshError, NativeMeshTag
else:
    warn(
        "The PyMOAB optional dependency could not be imported. "
        "All aspects of the partisn module are not imported.",
        ImportWarning,
    )

try:
    from pyne import dagmc

    HAVE_DAGMC = True
except:
    HAVE_DAGMC = False


def write_partisn_input(mesh, hdf5, ngroup, **kwargs):
    """This function reads a material-laden geometry file and a pre-made PyNE
    mesh object and writes a PARTISN text input file for blocks 1-5. The
    following cards are included:

    block 1: igeom, ngroup, niso,  mt, nzone, im, it, jm, jt, km, kt,
    block 2: xmesh, xints, ymesh, yints, zmesh, zints, zones,
    block 3: names,
    block 4: matls, assign,
    block 5: source.

    The 'source' card that appears by default is uniform in space and energy
    and isotropic in direction. In addition, "suggested" cards are printed,
    commented-out. These suggested cards are the additional cards required for a
    minimum working PARTISN input file:

    block 1: isn, maxscm, maxlcm,
    block 2: lib, lng, maxord, ihm, iht, ihs, ifido, ititl.

    Using the 'cards' parameter, any one of these cards (including 'source') as
    well as any PARTISN input card not specified here can be supplied and
    printed to the input file. Supplied cards will be ommitted from suggested
    cards.

    Parameters
    ----------
    mesh : PyNE mesh
        A premade mesh object that conforms to the geometry. Bounds of the mesh
        must correspond to the desired PARTISN coarse mesh intervals. By default
        one fine mesh inverval per coarse mesh will be used. This can be changed
        with the fine_per_coarse parameter. The sum of all fine mesh intervals
        in the problem must be greater than or equal to 7. Mesh can be 1-D
        (Nx1x1 mesh), 2-D (NxMx1 mesh), or 3-D (NxMxP mesh). Only Cartesian
        meshes are currently supported.
    hdf5 : string
        File path to a material-laden dagmc geometry file.
    ngroup : int
        The number of energy groups in the cross section library.
    input_file : string, optional, default = '<hdf5 file name>_partisn.inp'
        Desired path of generated PARTISN input file. Any file already existing
        by the same name will be overwritten.
    cards : dict, optional, default = {}
        This is a dictionary with the following keys: 'block1', 'block2',
        'block3', 'block4', 'block5'. The values are each dicts in the format:
         <partisn_card_name>:<partisn_card_value>. These cards will be printed
        out in the input file produced by this function. When specifying a
        source via this method, the key source be 'source' and the value should
        the entire source card, including the card name (e.g. source, sourcx,
        sourcef, etc) and '='.
    names_dict : dict, optional, default = None
        PyNE element/isotope names to bxslib name assignment. Keys are PyNE
        nucids (int) and values are bxslib names (str)
        Example: names_dict[250550000] ='mn55'
    num_rays : int, optional, default = 10
        For discretize_geom. Structured mesh only. The number of rays to fire
        in each mesh row for each direction.
    grid : boolean, optional, default = False
        For discretize_geom. Structured mesh only. If false, rays starting
        points are chosen randomly (on the boundary) for each mesh row. If
        true, a linearly spaced grid of starting points is used, with dimension
        sqrt(num_rays) x sqrt(num_rays). In this case, "num_rays" must be a
        perfect square.
    dg : record array, optional, default = None
        The output of pyne.dagmc.discretize_geom(). Use this input option if
        discretize_geom() has already been run, to avoid duplicating this
        expensive step. If HAVE_DAGMC=False, then this must be supplied.
    mat_assigns : dict, optional, default = None
        The output from pyne.cell_material_assignments().
        Dictionary of the cell to material assignments. Keys are cell
        numbers and values are material names. If HAVE_DAGMC=False, then
        this must be supplied.
    fine_per_coarse : int, optional, default = 1
        The number of fine mesh intervals per coarse mesh interval.
    data_hdf5path : string, optional, default = /materials
        the path in the heirarchy to the data table in an HDF5 file.
    nuc_hdf5path : string, optional, default = /nucid
        the path in the heirarchy to the nuclide array in an HDF5 file.
    """

    # Initialize dictionaries for each PARTISN block
    block01 = {}
    block02 = {}
    block03 = {}
    block04 = {}
    block05 = {}

    # Read optional inputs:
    cards = kwargs.get("cards", {})
    dg = kwargs.get("dg", None)
    mat_assigns = kwargs.get("mat_assigns", None)
    num_rays = kwargs.get("num_rays", 10)
    grid = kwargs.get("grid", False)
    if dg is not None and ("num_rays" in kwargs or "grid" in kwargs):
        warn("discretize_geom() options not used due to 'dg' argument")
    fine_per_coarse = kwargs.get("fine_per_coarse", 1)
    data_hdf5path = kwargs.get("data_hdf5path", "/materials")
    nuc_hdf5path = kwargs.get("nuc_hdf5path", "/nucid")
    title = hdf5.split("/")[-1].split(".")[0]
    input_file = kwargs.get("input_file", "{}_partisn.inp".format(title))

    # Dictionary of hdf5 names and cross section library names
    # Assumes PyNE naming convention in the cross section library if no dict
    # provided.
    if "names_dict" in kwargs:
        nuc_names = kwargs["names_dict"]
        mat_lib, unique_names = _get_material_lib(
            hdf5, data_hdf5path, nuc_names=nuc_names
        )
        mat_xs_names = _nucid_to_xs(mat_lib, nuc_names=nuc_names)
    else:
        mat_lib, unique_names = _get_material_lib(hdf5, data_hdf5path)
        mat_xs_names = _nucid_to_xs(mat_lib)

    # Set input variables
    block04["matls"] = mat_xs_names

    xs_names = _get_xs_names(mat_xs_names)
    block01["niso"] = len(xs_names)
    block03["names"] = xs_names

    block01["igeom"], bounds = _get_coord_sys(mesh)
    block01["ngroup"] = ngroup
    block01["mt"] = len(mat_lib)

    block02["zones"], block04["assign"] = _get_zones(
        mesh, hdf5, bounds, num_rays, grid, dg, mat_assigns, unique_names
    )
    block01["nzone"] = len(block04["assign"])
    block02["fine_per_coarse"] = fine_per_coarse

    for dim in bounds:
        if dim == "x":
            n = len(bounds[dim]) - 1
            block01["im"] = n
            block01["it"] = block01["im"] * fine_per_coarse
            block02["xmesh"] = bounds[dim]
        elif dim == "y":
            n = len(bounds[dim]) - 1
            block01["jm"] = n
            block01["jt"] = block01["jm"] * fine_per_coarse
            block02["ymesh"] = bounds[dim]
        elif dim == "z":
            n = len(bounds[dim]) - 1
            block01["km"] = n
            block01["kt"] = block01["km"] * fine_per_coarse
            block02["zmesh"] = bounds[dim]

    _check_fine_mesh_total(block01)

    # call function to write to file
    _write_input(title, block01, block02, block03, block04, block05, cards, input_file)


def _get_material_lib(hdf5, data_hdf5path, **kwargs):
    """Read material properties from the loaded dagmc geometry."""

    # If a set of nuc_names is provided, then collapse elements
    if "nuc_names" in kwargs:
        nuc_names = kwargs["nuc_names"]
        collapse = True
        # set of exception nuclides for collapse_elements
        mat_except = set(nuc_names.keys())
    else:
        collapse = False

    # collapse isotopes into elements (if required)
    mats = MaterialLibrary(hdf5, datapath=data_hdf5path)

    mats_collapsed = {}
    unique_names = {}

    for mat_name in mats:
        mat_name = mat_name.decode("utf-8")
        fluka_name = mats[mat_name].metadata["fluka_name"]
        if sys.version_info[0] > 2:
            unique_names[mat_name] = str(fluka_name.encode(), "utf-8")
        else:
            unique_names[mat_name] = fluka_name.decode("utf-8")

        if collapse:
            mats_collapsed[fluka_name] = mats[mat_name].collapse_elements(mat_except)
        else:
            mats_collapsed[fluka_name] = mats[mat_name]

    # convert mass fraction to atom density in units [at/b-cm]
    mat_lib = {}
    for mat_name in mats_collapsed:
        comp = mats_collapsed[mat_name]
        atom_dens_dict = comp.to_atom_dens()
        comp_list = {}
        for nucid, dens in atom_dens_dict.items():
            # convert from [at/cc] to [at/b-cm]
            comp_list[nucid] = dens * 10.0**-24
        mat_lib[mat_name] = comp_list

    return mat_lib, unique_names


def _nucid_to_xs(mat_lib, **kwargs):
    """Replace nucids with xs library names."""
    if "nuc_names" in kwargs:
        nuc_names = kwargs["nuc_names"]
        names_tf = True
    else:
        names_tf = False

    mat_xs_names = {}
    for mat in mat_lib:
        mat_xs_names[mat] = {}
        for nucid in mat_lib[mat]:
            if names_tf:
                if nucid in nuc_names:
                    name = nuc_names[nucid]
                    mat_xs_names[mat][name] = mat_lib[mat][nucid]
                else:
                    warn(
                        "Nucid {0} does not exist in the provided nuc_names dictionary.".format(
                            nucid
                        )
                    )
                    mat_xs_names[mat]["{0}".format(nucid)] = mat_lib[mat][nucid]
            else:
                mat_xs_names[mat][nucname.name(nucid)] = mat_lib[mat][nucid]

    return mat_xs_names


def _get_xs_names(mat_xs_names):
    """Create list of names (strings) of the nuclides that appear in the cross
    section library from the list of nuc_names.
    """

    xs_names = set()
    list(map(xs_names.update, mat_xs_names.values()))
    return list(xs_names)


def _get_coord_sys(mesh):
    """Determine coordinate system and get bounds"""

    # get number of divisions
    nx = len(mesh.structured_get_divisions("x"))
    ny = len(mesh.structured_get_divisions("y"))
    nz = len(mesh.structured_get_divisions("z"))

    coord_sys = ""
    if nx > 2:
        coord_sys += "x"
    if ny > 2:
        coord_sys += "y"
    if nz > 2:
        coord_sys += "z"

    # collect values of mesh boundaries for each coordinate
    bounds = {}
    fine = {}
    for i in coord_sys:
        bounds[i] = mesh.structured_get_divisions(i)

    # Determine IGEOM
    # assumes a Cartesian system
    if len(coord_sys) == 1:
        igeom = "slab"
    elif len(coord_sys) == 2:
        igeom = "x-y"
    elif len(coord_sys) == 3:
        igeom = "x-y-z"

    return igeom, bounds


def _get_zones(mesh, hdf5, bounds, num_rays, grid, dg, mat_assigns, unique_names):
    """Get the minimum zone definitions for the geometry."""

    # Discretize the geometry and get cell fractions
    if dg is None:
        if not HAVE_DAGMC:
            raise RuntimeError(
                "DAGMC is not available." "Unable to discretize the geometry."
            )
        else:
            dagmc.load(hdf5)
            dg = dagmc.discretize_geom(mesh, num_rays=num_rays, grid=grid)

    # Reorganize dictionary of each voxel's info with the key the voxel number
    # and values of cell and volume fraction
    voxel = {}
    for i in dg:
        idx = i[0]  # voxel number
        if idx not in voxel:
            voxel[idx] = {}
            voxel[idx]["cell"] = []
            voxel[idx]["vol_frac"] = []
        voxel[idx]["cell"].append(i[1])
        voxel[idx]["vol_frac"].append(i[2])

    # get material to cell assignments
    if mat_assigns is None:
        if not HAVE_DAGMC:
            raise RuntimeError(
                "DAGMC is not available." "Unable to get cell material assignments."
            )
        else:
            mat_assigns = dagmc.cell_material_assignments(hdf5)

    # Replace the names in the material assignments with unique names
    temp = {}
    for i, name in mat_assigns.items():
        if "vacuum" in name.lower() or "graveyard" in name.lower():
            temp[i] = name
        else:
            temp[i] = unique_names[name]
    mat_assigns = temp

    # Replace cell numbers with materials, eliminating duplicate materials
    # within single zone definition
    zones = {}
    for z in voxel:
        zones[z] = {}
        zones[z]["vol_frac"] = []
        zones[z]["mat"] = []
        for i, cell in enumerate(voxel[z]["cell"]):
            if mat_assigns[cell] not in zones[z]["mat"]:
                # create new entry
                zones[z]["mat"].append(mat_assigns[cell])
                zones[z]["vol_frac"].append(voxel[z]["vol_frac"][i])
            else:
                # update value that already exists with new volume fraction
                for j, val in enumerate(zones[z]["mat"]):
                    if mat_assigns[cell] == val:
                        vol_frac = zones[z]["vol_frac"][j] + voxel[z]["vol_frac"][i]
                        zones[z]["vol_frac"][j] = vol_frac

    # Remove vacuum or graveyard from material definition if not vol_frac of 1.0
    skip_array = [["mat:Vacuum"], ["mat:vacuum"], ["mat:Graveyard"], ["mat:graveyard"]]
    skip_list = ["mat:Vacuum", "mat:vacuum", "mat:Graveyard", "mat:graveyard"]
    zones_compressed = {}
    for z, info in zones.items():
        # check first if the definition is 100% void, keep same if is
        if zones[z]["mat"] in skip_array and zones[z]["vol_frac"] == [1.0]:
            zones_compressed[z] = info
        else:
            # check for partial void
            zones_compressed[z] = {"mat": [], "vol_frac": []}
            for i, mat in enumerate(zones[z]["mat"]):
                if mat not in skip_list:
                    zones_compressed[z]["mat"].append(mat)
                    zones_compressed[z]["vol_frac"].append(zones[z]["vol_frac"][i])

    # Eliminate duplicate zones and assign each voxel a zone number.
    # Assign zone = 0 if vacuum or graveyard and eliminate material definition.
    voxel_zone = {}
    zones_mats = {}
    z = 0
    match = False
    first = True
    for i, vals in zones_compressed.items():
        # Find if the zone already exists
        for zone, info in zones_mats.items():
            # Iterate through both sets to disregard order
            match_all = np.empty(len(vals["mat"]), dtype=bool)
            match_all.fill(False)
            for ii, mat in enumerate(vals["mat"]):
                for jj, mat_info in enumerate(info["mat"]):
                    if mat == mat_info and np.allclose(
                        np.array(vals["vol_frac"][ii]),
                        np.array(info["vol_frac"][jj]),
                        rtol=1e-5,
                    ):
                        match_all[ii] = True
                        break
            if match_all.all() == True:
                match = True
                y = zone
                break
            else:
                match = False
        # Create a new zone if first zone or does not match other zones
        if first or not match:
            # Check that the material is not 100% void (assign zone 0 otherwise)
            if vals["mat"] in skip_array:
                voxel_zone[i] = 0
            else:
                z += 1
                zones_mats[z] = zones_compressed[i]
                voxel_zone[i] = z
                first = False
        else:
            if vals["mat"] in skip_array:
                voxel_zone[i] = 0
            else:
                voxel_zone[i] = y

    # Remove any instances of graveyard or vacuum in zone definitions
    zones_novoid = {}
    for z in zones_mats:
        zones_novoid[z] = {"mat": [], "vol_frac": []}
        for i, mat in enumerate(zones_mats[z]["mat"]):
            if mat not in skip_list:
                zones_novoid[z]["mat"].append(mat)
                zones_novoid[z]["vol_frac"].append(zones_mats[z]["vol_frac"][i])

    # Put zones into format for PARTISN input
    if "x" in bounds:
        im = len(bounds["x"]) - 1
    else:
        im = 1

    if "y" in bounds:
        jm = len(bounds["y"]) - 1
    else:
        jm = 1

    if "z" in bounds:
        km = len(bounds["z"]) - 1
    else:
        km = 1

    n = 0
    zones_formatted = np.zeros(shape=(jm * km, im), dtype=int)
    for i in range(im):
        temp = np.zeros(shape=(jm * km), dtype=int)
        for jk in range(jm * km):
            temp[jk] = voxel_zone[n]
            n += 1
        temp = np.reshape(temp, (jm, km))
        temp = np.transpose(temp)
        temp = np.reshape(temp, jm * km)
        zones_formatted[:, i] = temp

    return zones_formatted, zones_novoid


def _check_fine_mesh_total(block01):
    """Check that the fine mesh total is greater than or equal to 7."""
    total = 0
    for key in block01:
        if key in ["it", "jt", "kt"]:
            total += block01[key]

    if total < 7:
        warn(
            "Please supply a larger mesh. Number of fine mesh intervals is less than 7."
        )


def _write_input(title, block01, block02, block03, block04, block05, cards, file_name):
    """Write all variables and comments to a file."""

    # Create file to write to
    f = open(file_name, "w")
    partisn = ""

    # NOTE: header is prepended at the end of this function.

    ###########################################
    #              Write Block 1              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 1 (Control and Dimensions) ------------"
    partisn += "\n/ \n"
    partisn += "igeom={0}".format(block01["igeom"])
    partisn += "  ngroup={0}".format(block01["ngroup"])
    partisn += "  niso={0}".format(block01["niso"])
    partisn += "  mt={0}".format(block01["mt"])
    partisn += "  nzone={0}\n".format(block01["nzone"])

    if "im" in block01:
        partisn += "im={0}".format(block01["im"])
        partisn += "  it={0}  ".format(block01["it"])
    if "jm" in block01:
        partisn += "jm={0}".format(block01["jm"])
        partisn += "  jt={0}  ".format(block01["jt"])
    if "km" in block01:
        partisn += "km={0}".format(block01["km"])
        partisn += "  kt={0}  ".format(block01["kt"])

    partisn += "\n"

    block1_cards = []
    if "block1" in cards:
        for card, value in sorted(cards["block1"].items()):
            partisn += "{}={}\n".format(card, value)
            block1_cards.append(card)

    missing_1 = set(["isn", "maxscm", "maxlcm"]) - set(block1_cards)
    if len(missing_1) > 0:
        partisn += "/ Please provide input for the following variables:\n"
        for mis in sorted(missing_1):
            partisn += "/{}=\n".format(mis)
    partisn += "t"

    ###########################################
    #              Write Block 2              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 2 (Geometry) ------------"
    partisn += "\n/ \n"

    if "xmesh" in block02:
        partisn += "xmesh= "
        count = 0
        for i, val in enumerate(block02["xmesh"]):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02["xmesh"]) - 1:
                    partisn += "\n       "
                count = 0
        partisn += "\nxints= "
        partisn += "{0}R {1}".format(
            len(block02["xmesh"]) - 1, block02["fine_per_coarse"]
        )
        partisn += "\n"

    if "ymesh" in block02:
        partisn += "ymesh= "
        count = 0
        for i, val in enumerate(block02["ymesh"]):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02["ymesh"]) - 1:
                    partisn += "\n       "
                count = 0
        partisn += "\nyints= "
        partisn += "{0}R {1}".format(
            len(block02["ymesh"]) - 1, block02["fine_per_coarse"]
        )
        partisn += "\n"

    if "zmesh" in block02:
        partisn += "zmesh= "
        count = 0
        for i, val in enumerate(block02["zmesh"]):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02["zmesh"]) - 1:
                    partisn += "\n       "
                count = 0
        partisn += "\nzints= "
        partisn += "{0}R {1}".format(
            len(block02["zmesh"]) - 1, block02["fine_per_coarse"]
        )
        partisn += "\n"

    partisn += "zones= "
    for i, row in enumerate(block02["zones"]):
        count = 0
        for num in row:
            partisn += "{} ".format(num)
            count += 1
            if count == 10:
                partisn += "\n       "
                count = 0
        partisn += ";"
        if i != len(block02["zones"]) - 1:
            partisn += "\n       "
        else:
            partisn += "\n"

    if "block2" in cards:
        for card, value in sorted(cards["block2"].items()):
            partisn += "{}={}\n".format(card, value)

    partisn += "t"

    ###########################################
    #              Write Block 3              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 3 (Nuclear Data) ------------"
    partisn += "\n/ \n"

    partisn += "/ Note: NAMES is not all inclusive. Only NAMES that are present in\n"
    partisn += "/ meshed area are listed.\n"
    partisn += "names= "
    count = 0
    for i, name in enumerate(sorted(block03["names"])):
        count += 1
        partisn += "{0} ".format(name)
        if count == 10:
            if i != len(block03["names"]) - 1:
                partisn += "\n       "
            count = 0

    partisn += "\n"

    block3_cards = []
    if "block3" in cards:
        for card, value in sorted(cards["block3"].items()):
            partisn += "{}={}\n".format(card, value)
            block3_cards.append(card)

    missing_3 = set(
        ["lib", "lng", "maxord", "ihm", "iht", "ihs", "ifido", "ititl"]
    ) - set(block3_cards)
    if len(missing_3) > 0:
        partisn += "/ Please provide input for the following variables:\n"
        for mis in sorted(missing_3):
            partisn += "/{}=\n".format(mis)
    partisn += "t"

    ###########################################
    #              Write Block 4              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 4 (Cross-Section Mixing) ------------"
    partisn += "\n/ \n"

    partisn += "matls= "
    for i, mat in enumerate(sorted(block04["matls"])):
        partisn += "{0} ".format(mat)
        count = 0
        j = 0
        for iso, dens in sorted(block04["matls"][mat].items()):
            count += 1
            j += 1
            if j != len(block04["matls"][mat]):
                partisn += "{} {:.4e}, ".format(iso, dens)
                if count == 3:
                    if j != len(block04["matls"][mat]):
                        partisn += "\n       "
                    count = 0
            else:
                if i == len(block04["matls"]) - 1:
                    partisn += "{} {:.4e};\n".format(iso, dens)
                else:
                    partisn += "{} {:.4e};\n       ".format(iso, dens)

    partisn += "assign= "
    for i, z in enumerate(block04["assign"]):
        partisn += "{0} ".format(z)
        count = 0
        for j, mat in enumerate(block04["assign"][z]["mat"]):
            if j != len(block04["assign"][z]["mat"]) - 1:
                count += 1
                partisn += "{} {:.4e}, ".format(
                    mat, block04["assign"][z]["vol_frac"][j]
                )
                if count == 3:
                    if i != len(block04["assign"][z]["mat"]) - 1:
                        partisn += "\n          "
                    count = 0
            else:
                if i == len(block04["assign"]) - 1:
                    partisn += "{} {:.4e};\n".format(
                        mat, block04["assign"][z]["vol_frac"][j]
                    )
                else:
                    partisn += "{} {:.4e};\n        ".format(
                        mat, block04["assign"][z]["vol_frac"][j]
                    )

    if "block4" in cards:
        for card, value in cards["block4"].items():
            partisn += "{}={}\n".format(card, value)

    partisn += "t"

    ###########################################
    #              Write Block 5              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 5 (Solver Inputs) ------------"
    partisn += "\n/ \n"
    if "block5" in cards and "source" in cards["block5"]:
        partisn += cards["block5"]["source"]
        if partisn[-1] != "\n":
            partisn += "\n"
        default_source = False
    else:
        # default source
        partisn += "source={}R 1\n".format(block01["ngroup"])
        default_source = True

    if "block5" in cards:
        for card, value in cards["block5"].items():
            if card != "source":
                partisn += "{}={}\n".format(card, value)
    partisn += "t\n"

    ###########################################
    #              Write Header               #
    ###########################################
    header = "     1     0     0\n"
    header += "{}\n".format(title)
    header += "/\n"
    if default_source:
        header += "/ NOTE: This input includes a default source that is isotropic\n"
        header += "/       in direction and uniform in space and energy.\n"
    if len(missing_1) > 0 or len(missing_3) > 0:
        header += "/ NOTE: The follow commented out cards must be filled in for\n"
        header += "/       a complete PARTISN input file:\n"
        if len(missing_1) > 0:
            header += "/       Block 1:"
            for mis in sorted(missing_1):
                header += " {},".format(mis)
            header += "\n"
        if len(missing_3) > 0:
            header += "/       Block 3:"
            for mis in sorted(missing_3):
                header += " {},".format(mis)
            header += "\n"
    header += "/"
    # Prepend header to begining of file
    partisn = header + partisn

    # Write to the file
    f.write(partisn)


def format_repeated_vector(vector):
    """Creates string out of a vector with the PARTISN format for repeated
    numbers.

    Parameters:
    -----------
    vector: list
        Desired list to be formatted

    Returns:
    --------
    string: string
        Formatted string representation of the vector

    Example:
        vector = [1, 2, 0, 0, 0, 7, 8, 3, 3]
        string = "1 2 3R 0 7 8 2R 3"
    """

    # put vector into a list of lists formatted as
    # [[number , R], [number, R], ...]
    # where 'R' is the number of times that 'number' is repeated
    tot = 0
    repeats = []
    for i, val in enumerate(vector):
        if tot == 0:
            repeats.append([val, 1])
            tot += 1
        else:
            if val == repeats[tot - 1][0]:
                repeats[tot - 1][1] += 1
            else:
                repeats.append([val, 1])
                tot += 1

    # make into a string of characters
    string = ""
    n = 0
    for pair in repeats:
        if pair[1] == 1:
            string += "{} ".format(pair[0])
            n = +1
        else:
            string += "{0}R {1} ".format(pair[1], pair[0])
            n += 2

    return string


def mesh_to_isotropic_source(m, tag):
    """This function reads an isotropic source definition from a supplied mesh
    and creates a corresponding PARTISN SOURCF input card. The input card wraps
    to 80 characters and utilizes the "R" repeation notation for zero values
    (e.g. 4R 0 = 0 0 0 0).

    Parameters:
    -----------
    m : PyNE Mesh
        The mesh tagged with an energy-dependent source distribution.
    tag : str
        The tag on the mesh with the source information.

    Returns:
    --------
    s : str
        SOURF input card wrapped to 80 characters.
    """

    # get data
    temp = m.structured_ordering
    m.structured_ordering = "zyx"
    m.src = NativeMeshTag(name=tag)
    data = m.src[:].transpose()[::-1]
    m.structured_ordering = temp
    ninti = len(m.structured_coords[0]) - 1

    # format output
    s = "sourcf="
    count = 1
    zero_count = 0
    for e_row in data:
        for src in e_row:
            if src == 0.0:
                zero_count += 1
            else:
                if zero_count != 0:
                    if zero_count == 1:
                        s += " 0"
                    else:
                        s += " {}R 0".format(zero_count)
                    zero_count = 0
                s += " {0:9.5E}".format(src)
            if count % ninti == 0:
                if zero_count != 0:
                    if zero_count == 1:
                        s += " 0"
                    else:
                        s += " {}R 0".format(zero_count)
                    zero_count = 0
                s += ";"
            count += 1

    # wrap to 80 characters
    s = "\n".join(wrap(s, 80))
    return s


def isotropic_vol_source(geom, mesh, cells, spectra, intensities, **kwargs):
    """This function creates an isotropic volumetric source within each
    requested geometry cell of a DAGMC CAD geometry. This is done by Monte
    Carlo ray tracing to determine volume fractions of each source cell within
    each mesh volume element. The PARTISN SOURCF card is returned, as well as
    the ray-traced volume fractions (dagmc.discretize_geom output), so that ray
    tracing does not have to be done twice when using this function with
    write_partisn_input. The supplied mesh is also tagged with the calculated
    source distribution using this function.

    Parameters:
    -----------
    geom : str
        The DAGMC geometry (.h5m) file containing the geometry of interest.
    mesh : PyNE Mesh
        The superimposed Cartesian mesh that will be used to define the source
        for PARTISN transport.
    cells : list of ints
        The cell numbers of DAGMC geometry cells which have non-zero source
        intensity.
    spectra : list of list of floats
        The normalized energy spectrum for each of the cells. If spectra are
        not normalized, they will be normalized in this function.
    intensities : list of floats
        The volumetric intensity (i.e. s^-1 cm^-3) for each geometry cell.
    tag_name : str, optional, default = 'src'
        The name of the tag for which source data will be tagged on the mesh.
    num_rays : int, optional, default = 10
        For discretize_geom. The number of rays to fire in each mesh row for
        each direction.
    grid : boolean, optional, default = False
        For discretize_geom. If false, rays starting points are chosen randomly
        (on the boundary) for each mesh row. If true, a linearly spaced grid of
        starting points is used, with dimension sqrt(num_rays) x sqrt(num_rays).
        In this case, "num_rays" must be a perfect square.

    Returns:
    --------
    output : str
        PARTISN SOURF card representing the requested source
    dg : record array
        The output of dagmc.discretize_geom; stored in a one dimensional array,
        each entry containing the following
        fields:
        :idx: int
            The volume element index.
        :cell: int
            The geometry cell number.
        :vol_frac: float
            The volume fraction of the cell withing the mesh ve.
        :rel_error: float
            The relative error associated with the volume fraction.
        This array is returned in sorted order with respect to idx and cell, with
        cell changing fastest.
    """
    # discretize_geom inputs
    tag_name = kwargs.get("tag_name", "src")
    num_rays = kwargs.get("num_rays", 10)
    grid = kwargs.get("grid", False)

    # Check lengths of input
    if len(cells) != len(spectra) or len(cells) != len(intensities):
        raise ValueError("Cells, spectra, intensities must be the same length")
    lengths = [len(x) for x in spectra]
    if not all(lengths[0] == length for length in lengths):
        raise ValueError("Spectra must all be the same length")

    # Normalize spectra
    norm_spectra = []
    for spec in spectra:
        total = np.sum(spec)
        norm_spectra.append(np.array(spec) / total)

    norm_spectra = {cell: spec for cell, spec in zip(cells, norm_spectra)}
    intensities = {cell: inten for cell, inten in zip(cells, intensities)}

    # ray trace
    if not HAVE_DAGMC:
        raise RuntimeError(
            "DAGMC is not available." "Cannot run isotropic_vol_source()."
        )
    else:
        dagmc.load(geom)
        dg = dagmc.discretize_geom(mesh, num_rays=num_rays, grid=grid)

    # determine  source intensities
    data = np.zeros(shape=(len(mesh), len(spectra[0])))
    for row in dg:
        if row[1] in cells:
            data[row[0], :] += np.multiply(
                row[2] * intensities[row[1]], norm_spectra[row[1]]
            )

    mesh.tag = NativeMeshTag(len(spectra[0]), float, name=tag_name)
    mesh.tag[:] = data

    output = mesh_to_isotropic_source(mesh, tag_name)
    return output, dg
