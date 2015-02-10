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
#from sets import Set

import numpy as np
import tables

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
                  "All aspects of the partisn module are not imported.",
                  QAWarning)
    HAVE_PYTAPS = False

if HAVE_PYTAPS:
    from pyne.mesh import Mesh, StatMesh, MeshError, IMeshTag


def write_partisn_input(mesh, hdf5, ngroup, pn, **kwargs):
    """This definition reads all necessary attributes from a material-laden 
    geometry file, a pre-made PyNE mesh object, and the nuclear data cross 
    section library, and any optional inputs that are necessary for creating a 
    PARTISN input file. It then writes a PARTISN text input file for blocks 1-5.
    Note that comments appear in the created input file where more variables
    must be set. 
    
    Notes:
        This does not write out all necessary inputs for the solver and cross
        section library (block 3 and 5). There is no assumed cross section 
        library type.
    
    Parameters:
    -----------
    mesh : PyNE mesh 
        A premade mesh object that conforms to the geometry. Bounds of the mesh
        must correspond to the desired PARTISN fine mesh intervals. Two fine 
        mesh intervals per coarse mesh interval will be created. The sum of all
        fine mesh intervals in the problem must be greater than or equal to 7.
        Mesh can be 1-D (Nx1x1 mesh), 2-D (NxMx1 mesh), or 3-D (NxMxP mesh).
        Note: Only Cartesian meshes are currently supported.
    hdf5 : string
        File path to a material-laden dagmc geometry file.
    ngroup : int
        The number of energy groups in the cross section library.
    pn : int
        The number of moments in a P_n expansion of the source.
    data_hdf5path : string, optional, default = material_library/materials
        the path in the heirarchy to the data table in an HDF5 file.
    nuc_hdf5path : string, optional, default = material_library/nucid
        the path in the heirarchy to the nuclide array in an HDF5 file.
    names_dict : dict, optional
        PyNE element/isotope names to bxslib name assignment. Keys are PyNE
        nucids (int) and values are bxslib names (str)
        Example: names_dict[250550000] ='mn55'
    input_file : string, optional, default = '<hdf5 file name>_partisn.inp'
        Desired path of generated PARTISN input file. Any file already existing
        by the same name will be overwritten.
    num_rays : int, optional, default = 10
        For discretize_geom. Structured mesh only. The number of rays to fire 
        in each mesh row for each direction.
    grid : boolean, optional, default = False
        For discretize_geom. Structured mesh only. If false, rays starting 
        points are chosen randomly (on the boundary) for each mesh row. If 
        true, a linearly spaced grid of starting points is used, with dimension 
        sqrt(num_rays) x sqrt(num_rays). In this case, "num_rays" must be a 
        perfect square.
    
    Output:
    -------
    PARTISN Input file named by 'input_file' above or the default name.
        Note: read comments generated in file. Not all variables will be 
        assigned that are necessary.
    """
    
    # Initialize dictionaries for each PARTISN block
    block01 = {}
    block02 = {}
    block03 = {}
    block04 = {}
    block05 = {}
    
    # Read optional inputs:
    
    # discretize_geom inputs
    if 'num_rays' in kwargs:
        num_rays = kwargs['num_rays']
    else:
        num_rays = 10
    
    if 'grid' in kwargs:
        grid = kwargs['grid']
    else:
        grid = False
    
    # hdf5 paths
    if 'data_hdf5path' in kwargs:
        data_hdf5path = kwargs['data_hdf5path']  
    else:
        data_hdf5path = '/material_library/materials'
    
    if 'nuc_hdf5path' in kwargs:
        nuc_hdf5path = kwargs['nuc_hdf5path']
    else:
        nuc_hdf5path = '/material_library/nucid'
    
    # input file name
    if 'input_file' in kwargs:
        input_file = kwargs['input_file']
        input_file_tf = True
    else:
        input_file_tf = False
    
    # Dictionary of hdf5 names and cross section library names
    # Assumes PyNE naming convention in the cross section library if no dict
    # provided.
    if 'names_dict' in kwargs:
        nuc_names = kwargs['names_dict']
        mat_lib = _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, nuc_names=nuc_names)
        mat_xs_names = _nucid_to_xs(mat_lib, nuc_names=nuc_names)
    else:
        mat_lib = _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path)
        mat_xs_names = _nucid_to_xs(mat_lib)
    
    # Set input variables
    
    block04['matls'] = mat_xs_names
    
    xs_names = _get_xs_names(mat_xs_names)
    block01['niso'] = len(xs_names)
    block03['names'] = xs_names

    block01['igeom'], bounds, nmq = _get_coord_sys(mesh, pn)
    block01['ngroup'] = ngroup
    block01['mt'] = len(mat_lib)
    
    block02['zones'], block04['assign'] = _get_zones(mesh, hdf5, bounds, num_rays, grid)
    block01['nzone'] = len(block04['assign'])
    
    for dim in bounds:
        if dim == 'x':
            n = len(bounds[dim]) - 1
            block01['im'] = n
            block01['it'] = block01['im']*2
            block02['xmesh'] = bounds[dim]
            block05['sourcx'] = np.zeros(shape=(nmq, n), dtype=float)
            block05['sourcx'][:,0] = 1.0
        elif dim == 'y':
            n = len(bounds[dim]) - 1
            block01['jm'] = n
            block01['jt'] = block01['jm']*2
            block02['ymesh'] = bounds[dim]
            block05['sourcy'] = np.zeros(shape=(nmq, n), dtype=float)
            block05['sourcy'][:,0] = 1.0
        elif dim == 'z':
            n = len(bounds[dim]) - 1
            block01['km'] = n
            block01['kt'] = block01['km']*2
            block02['zmesh'] = bounds[dim]
            block05['sourcz'] = np.zeros(shape=(nmq, n), dtype=float)
            block05['sourcz'][:,0] = 1.0
    
    warn_fm = _check_fine_mesh_total(block01)
    if warn_fm:
        warn("Please supply a larger mesh. Number of fine mesh intervals is less than 7.")

    block05['source'] = np.zeros(shape=(nmq, ngroup), dtype=float)
    block05['source'][:,0] = 1.0
    
    # create title
    if "/" in hdf5:
        title = hdf5.split("/")[len(hdf5.split("/"))-1].split(".")[0]
    else:
        title = hdf5.split(".")[0]
    
    # call function to write to file
    input_file = input_file if input_file_tf else None
    _write_input(title, block01, block02, block03, block04, block05, name=input_file)


def _get_material_lib(hdf5, data_hdf5path, nuc_hdf5path, **kwargs):
    """Read material properties from the loaded dagmc geometry.
    """
    
    # If a set of nuc_names is provided, then collapse elements
    if 'nuc_names' in kwargs:
        nuc_names = kwargs['nuc_names']
        collapse = True
        # set of exception nuclides for collapse_elements
        mat_except = set(nuc_names.keys())
    else:
        collapse = False
    
    # collapse isotopes into elements (if required)
    mats = MaterialLibrary(hdf5, datapath=data_hdf5path, nucpath=nuc_hdf5path)
    mats_collapsed = {}
    for mat_name in mats:
        if collapse:
            mats_collapsed[mat_name] = mats[mat_name].collapse_elements(mat_except)
        else:
            mats_collapsed[mat_name] = mats[mat_name]

    # convert mass fraction to atom density in units [at/b-cm]
    mat_lib = {}
    for mat_name in mats_collapsed:
        comp = mats_collapsed[mat_name]
        atom_dens_dict = comp.to_atom_dens()
        comp_list = {}
        for nucid, dens in atom_dens_dict.iteritems():
            # convert from [at/cc] to [at/b-cm]
            comp_list[nucid] = dens*10.**-24
        mat_lib[mat_name] = comp_list

    return mat_lib


def _nucid_to_xs(mat_lib, **kwargs):
    """Replace nucids with xs library names.
    """
    
    if 'nuc_names' in kwargs:
        nuc_names = kwargs['nuc_names']
        names_tf = True
    else:
        names_tf = False
    
    mat_xs_names = {}
    for mat in mat_lib:
        mat_name = strip_mat_name(mat)
        mat_xs_names[mat_name] = {}
        for nucid in mat_lib[mat]:
            if names_tf:
                if nucid in nuc_names:
                    name = nuc_names[nucid]
                    mat_xs_names[mat_name][name] = mat_lib[mat][nucid]
                else:
                    warn("Nucid {0} does not exist in the provided nuc_names dictionary.".format(nucid))
                    mat_xs_names[mat_name]["{0}".format(nucid)] = mat_lib[mat][nucid]
            else:
                mat_xs_names[mat_name][nucname.name(nucid)] = mat_lib[mat][nucid]

    return mat_xs_names
    

def _get_xs_names(mat_xs_names):
    """Create list of names (strings) of the nuclides that appear in the cross
    section library from the list of nuc_names.
    """
    
    xs_names = set()
    list(map(xs_names.update, mat_xs_names.values()))
    return list(xs_names)


def _get_coord_sys(mesh, pn):
    """Determine coordinate system and get bounds
    """
    
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
        igeom = 'slab'
        nmq = pn + 1
    elif len(coord_sys) == 2:
        igeom = 'x-y'
        nmq = (pn + 1)*(pn + 2)/2
    elif len(coord_sys) == 3:
        igeom = 'x-y-z'
        nmq = (pn + 1)**2
    
    return igeom, bounds, nmq


def _get_zones(mesh, hdf5, bounds, num_rays, grid):
    """Get the minimum zone definitions for the geometry.
    """
    
    # load the geometry
    from pyne import dagmc
    dagmc.load(hdf5)
    
    # Descretize the geometry and get cell fractions
    dg = dagmc.discretize_geom(mesh, num_rays=num_rays, grid=grid)

    # Reorganize dictionary of each voxel's info with the key the voxel number 
    # and values of cell and volume fraction   
    voxel = {}
    for i in dg:
        idx = i[0]  # voxel number
        if idx not in voxel:
            voxel[idx] = {}
            voxel[idx]['cell'] = []
            voxel[idx]['vol_frac'] = []
        voxel[idx]['cell'].append(i[1])
        voxel[idx]['vol_frac'].append(i[2])

    # get material to cell assignments
    mat_assigns = dagmc.cell_material_assignments(hdf5)

    # Replace cell numbers with materials, eliminating duplicate materials
    # within single zone definition
    zones = {}
    for z in voxel:
        zones[z] = {}
        zones[z]['vol_frac'] = []
        zones[z]['mat'] = []
        for i, cell in enumerate(voxel[z]['cell']):
            if mat_assigns[cell] not in zones[z]['mat']:
                # create new entry
                zones[z]['mat'].append(mat_assigns[cell])
                zones[z]['vol_frac'].append(voxel[z]['vol_frac'][i])
            else:
                # update value that already exists with new volume fraction
                for j, val in enumerate(zones[z]['mat']):
                    if mat_assigns[cell] == val:
                        vol_frac = zones[z]['vol_frac'][j] + voxel[z]['vol_frac'][i]
                        zones[z]['vol_frac'][j] = vol_frac
    
    # Remove vacuum or graveyard from material definition if not vol_frac of 1.0
    skip_array = [['mat:Vacuum'], ['mat:vacuum'], ['mat:Graveyard'], ['mat:graveyard']]
    skip_list = ['mat:Vacuum', 'mat:vacuum', 'mat:Graveyard', 'mat:graveyard']
    zones_compressed = {}
    for z, info in zones.iteritems():
        # check first if the definition is 100% void, keep same if is
        if zones[z]['mat'] in skip_array and zones[z]['vol_frac'] == [1.0]:
            zones_compressed[z] = info
        else:
            # check for partial void
            zones_compressed[z] = {'mat':[], 'vol_frac':[]}
            for i, mat in enumerate(zones[z]['mat']):
                if mat not in skip_list:
                    zones_compressed[z]['mat'].append(mat)
                    zones_compressed[z]['vol_frac'].append(zones[z]['vol_frac'][i])
    
    # Eliminate duplicate zones and assign each voxel a zone number.
    # Assign zone = 0 if vacuum or graveyard and eliminate material definition.
    voxel_zone = {}
    zones_mats = {}
    z = 0
    match = False
    first = True    
    for i, vals in zones_compressed.iteritems():
        # Find if the zone already exists
        for zone, info in zones_mats.iteritems():
            # Iterate through both sets to disregard order
            match_all = np.empty(len(vals['mat']), dtype=bool)
            match_all.fill(False)
            for ii, mat in enumerate(vals['mat']):
                for jj, mat_info in enumerate(info['mat']):
                    if mat == mat_info and np.allclose(np.array(vals['vol_frac'][ii]), \
                                np.array(info['vol_frac'][jj]), rtol=1e-5):
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
            if vals['mat'] in skip_array:
                voxel_zone[i] = 0
            else:
                z += 1
                zones_mats[z] = zones_compressed[i]
                voxel_zone[i] = z
                first = False
        else:
            if vals['mat'] in skip_array:
                voxel_zone[i] = 0
            else:
                voxel_zone[i] = y
    
    # Remove any instances of graveyard or vacuum in zone definitions
    zones_novoid = {}
    for z in zones_mats:
        zones_novoid[z] = {'mat':[], 'vol_frac':[]}
        for i, mat in enumerate(zones_mats[z]['mat']):
            if mat not in skip_list:
                name = strip_mat_name(mat)
                zones_novoid[z]['mat'].append(name)
                zones_novoid[z]['vol_frac'].append(zones_mats[z]['vol_frac'][i])
    
    # Put zones into format for PARTISN input
    if 'x' in bounds:
        im = len(bounds['x']) - 1
    else:
        im = 1
    
    if 'y' in bounds:
        jm = len(bounds['y']) - 1
    else:
        jm = 1
    
    if 'z' in bounds:
        km = len(bounds['z']) - 1
    else:
        km = 1

    n = 0
    zones_formatted = np.zeros(shape=(jm*km, im), dtype=int)
    for i in range(im):
        for jk in range(jm*km):
            zones_formatted[jk, i] = voxel_zone[n]
            n += 1
    
    return zones_formatted, zones_novoid
    

def _check_fine_mesh_total(block01):
    """Check that the fine mesh total is greater than or equal to 7.
    """
    total = 0
    for key in block01:
        if key in ['it', 'jt', 'kt']:
            total += block01[key]
    
    if total >= 7:
        # no warning necessary
        return False
    else:
        # warn the user
        return True


def _write_input(title, block01, block02, block03, block04, block05, name=None):
    """Write all variables and comments to a file.
    """
    
    # Create file to write to
    file_name = str(title) + '_partisn.inp' if name is None else name
    f = open(file_name, 'w')
    partisn = ''
    
    # Write title
    partisn += "     1     0     0\n"
    partisn += str(title)

    partisn += "\n/ "
    partisn += "\n/ Notes: This input assumes a volumetric source calculation using"
    partisn += "\n/ default PARTISN values in many cases. Please refer to the comments"
    partisn += "\n/ throughout and the PARTISN input manual."
    partisn += "\n/ Variables that MUST be set in each block (other defaults and \n"
    partisn += "/ optional variables may exist):"
    partisn += "\n/     Block 1:  ISN"
    partisn += "\n/     Block 3:  LIB, MAXORD, IHM, IHT"
    partisn += "\n/     Block 6:  no input is provided for block 6"
    
    ###########################################
    #              Write Block 1              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 1 (Control and Dimensions) ------------"
    partisn += "\n/ \n"
    partisn += "igeom={0}".format(block01['igeom'])
    partisn += "  ngroup={0}".format(block01['ngroup'])
    partisn += "  niso={0}".format(block01['niso'])
    partisn += "  mt={0}".format(block01['mt'])
    partisn += "  nzone={0}\n".format(block01['nzone'])
    
    partisn += "/ Please provide input for ISN variable:\n"
    partisn += "/ isn=  \n"
    
    if 'im' in block01:
        partisn += "im={0}".format(block01['im'])
        partisn += "  it={0}  ".format(block01['it'])
    if 'jm' in block01:
        partisn += "jm={0}".format(block01['jm'])
        partisn += "  jt={0}  ".format(block01['jt'])
    if 'km' in block01:
        partisn += "km={0}".format(block01['km'])
        partisn += "  kt={0}  ".format(block01['kt'])

    partisn += "\n"
    partisn += 't'
    
    ###########################################
    #              Write Block 2              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 2 (Geometry) ------------"
    partisn += "\n/ \n"
    
    if 'xmesh' in block02:
        partisn += "xmesh= "
        count = 0
        for i, val in enumerate(block02['xmesh']):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02['xmesh'])-1:
                    partisn += "\n       "
                count = 0
        partisn += "\nxints= "
        partisn += "{0}R {1}".format(len(block02['xmesh'])-1, 2)
        partisn += "\n"
        
    if 'ymesh' in block02:
        partisn += "ymesh= "
        count = 0
        for i, val in enumerate(block02['ymesh']):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02['ymesh'])-1:
                    partisn += "\n       "
                count = 0
        partisn += "\nyints= "
        partisn += "{0}R {1}".format(len(block02['ymesh'])-1, 2)
        partisn += "\n"
        
    if 'zmesh' in block02:
        partisn += "zmesh= "
        count = 0
        for i, val in enumerate(block02['zmesh']):
            count += 1
            partisn += "{:.3f} ".format(val)
            if count == 8:
                if i != len(block02['zmesh'])-1:
                    partisn += "\n       "
                count = 0
        partisn += "\nzints= "
        partisn += "{0}R {1}".format(len(block02['zmesh'])-1, 2)
        partisn += "\n"
        
    partisn += "zones= "
    for i, row in enumerate(block02['zones']):
        count = 0
        for num in row:
            partisn += "{} ".format(num)
            count += 1
            if count == 20:
                partisn += "\n       "
                count = 0
        partisn += ";"
        if i != len(block02['zones'])-1:
            partisn += "\n       "
        else:
            partisn += "\n"

    partisn += "t"
    
    ###########################################
    #              Write Block 3              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 3 (Nuclear Data) ------------"
    partisn += "\n/ \n"
    partisn += "/ Please provide input for the following variables:\n"
    partisn += "/ lib=\n"
    partisn += "/ maxord=\n"
    partisn += "/ ihm=\n"
    partisn += "/ iht=\n"
    
    partisn += "/ Note: NAMES is not all inclusive. Only NAMES that are present in\n"
    partisn += "/ meshed area are listed.\n"
    partisn += "names= "
    count = 0
    for i, name in enumerate(block03['names']):
        count += 1
        partisn += "{0} ".format(name)
        if count == 10:
            if i != len(block03['names'])-1:
                partisn += "\n       "
            count = 0
    
    partisn += "\nt"
    
    ###########################################
    #              Write Block 4              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 4 (Cross-Section Mixing) ------------"
    partisn += "\n/ \n"
    
    partisn += "matls= "
    for i, mat in enumerate(block04['matls']):
        partisn += "{0} ".format(mat)
        count = 0
        j = 0
        for iso, dens in block04['matls'][mat].iteritems():
            count += 1
            j += 1
            if j != len(block04['matls'][mat]):
                partisn += "{} {:.4e}, ".format(iso, dens)
                if count == 3:
                    if j != len(block04['matls'][mat]):
                        partisn += "\n       "
                    count = 0
            else:
                if i == len(block04['matls']) - 1:
                    partisn += "{} {:.4e};\n".format(iso, dens)
                else:
                    partisn += "{} {:.4e};\n       ".format(iso, dens)
            
    partisn += "assign= "
    for i, z in enumerate(block04['assign']):
        partisn += "{0} ".format(z)
        count = 0
        for j, mat in enumerate(block04['assign'][z]['mat']):
            if j != len(block04['assign'][z]['mat'])-1:
                count += 1
                partisn += "{} {:.4e}, ".format(mat, block04['assign'][z]['vol_frac'][j])
                if count == 3:
                    if i != len(block04['assign'][z]['mat'])-1:
                        partisn += "\n          "
                    count = 0
            else:
                if i == len(block04['assign']) - 1:
                    partisn += "{} {:.4e};\n".format(mat, block04['assign'][z]['vol_frac'][j])
                else:
                    partisn += "{} {:.4e};\n        ".format(mat, block04['assign'][z]['vol_frac'][j])
    
    partisn += "t"
    
    ###########################################
    #              Write Block 5              #
    ###########################################
    partisn += "\n/ \n"
    partisn += "/ ------------ Block 5 (Solver Inputs) ------------"
    partisn += "\n/ \n"
    partisn += "/ This input assumes a volumetric source calculation with vacuum\n"
    partisn += "/ boundary conditions. Change inputs below if otherwise.\n"
    partisn += "ievt=0      / source calculation\n"
    partisn += "/ isct=     / Legendre order of scattering (default=0)\n"
    partisn += "/ ith=      / 0/1/2= direct/adjoint/POI calculation (default=0)\n"
    partisn += "/ ibl=      / left BC (default=0, vacuum)\n"
    partisn += "/ ibr=      / right BC (default=0, vacuum)\n"
    partisn += "/ ibt=      / top BC (default=0, vacuum)\n"
    partisn += "/ ibb=      / bottom BC (default=0, vacuum)\n"
    partisn += "/ ibfrnt=   / front BC (default=0, vacuum)\n"
    partisn += "/ ibback=   / back BC (default=0, vacuum)\n"
    partisn += "/ \n"
    
    partisn += "/ Source is in format of option 3 according to PARTISN input manual.\n"
    partisn += "/ Default is an evenly distributed volume source.\n"
    partisn += "source= "
    count = 0
    tot = 0
    for row in block05['source']:
        partisn += format_repeated_vector(row)
        tot += 1
        count += 1
        if count == 4:
            if tot != len(block05['source']):
                partisn += ";\n        "
            else:
                partisn += ";"
            count = 0
        else:
            partisn += "; "
    partisn += "\n"
        
    if 'sourcx' in block05:
        partisn += "sourcx= "
        count = 0
        tot = 0
        for row in block05['sourcx']:
            partisn += format_repeated_vector(row)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcx']):
                    partisn += ";\n        "
                else:
                    partisn += ";"
                count = 0
            else:
                partisn += "; "
        partisn += "\n"
                
    if 'sourcy' in block05:
        partisn += "sourcy= "
        count = 0
        tot = 0
        for row in block05['sourcy']:
            partisn += format_repeated_vector(row)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcy']):
                    partisn += ";\n        "
                else:
                    partisn += ";"
                count = 0
            else:
                partisn += "; "
        partisn += "\n"
            
    if 'sourcz' in block05:
        partisn += "sourcz= "
        count = 0
        tot = 0
        for row in block05['sourcz']:
            partisn += format_repeated_vector(row)
            tot += 1
            count += 1
            if count == 4:
                if tot != len(block05['sourcz']):
                    partisn += ";\n        "
                else:
                    partisn += ";"
                count = 0
            else:
                partisn += "; "
        partisn += "\n"
    
    partisn += "t\n"
    
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
            if val == repeats[tot-1][0]:
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
            n =+ 1
        else:
            string += "{0}R {1} ".format(pair[1], pair[0])
            n += 2

    return string


def strip_mat_name(mat_name):
    """Provide a material name (string) and receive a compacted name without 
    'mat:' or special characters.
    Assumes PyNE naming convention (must start with 'mat:').
    """
    
    # Remove 'mat:'
    tmp1 = mat_name.split(':')[1]
    
    # Remove other special characters
    special_char = [':', ',', ' ', ';', "'", '"']
    for char in special_char:
        tmp2 = tmp1.split(char)
        tmp1 = ''.join(tmp2)
    
    return tmp1
