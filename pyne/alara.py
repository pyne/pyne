"""This module contains functions relevant to the ALARA activation code.
"""
from __future__ import print_function
import os
import collections
import numpy as np
import tables as tb
import warnings

try:
    from itaps import iMesh, iBase, iMeshExtensions
except ImportError:
    warnings.warn("the PyTAPS optional dependency could not be imported. "
                  "Some aspects of the alara module may be incomplete.",
                  ImportWarning)

from pyne.mesh import Mesh, MeshError
from pyne.material import Material, from_atom_frac
from pyne.nucname import serpent, alara, znum, anum
from pyne.data import N_A

def mesh_to_fluxin(flux_mesh, flux_tag, fluxin="fluxin.out",
                   reverse=False):
    """This function creates an ALARA fluxin file from fluxes tagged on a PyNE
    Mesh object. Structured meshes are printed in xyz order (z changes fastest)
    and unstructured meshes are printed in the imesh.iterate() order.

    Parameters
    ----------
    flux_mesh : PyNE Mesh object
        Contains the mesh with fluxes tagged on each volume element.
    flux_tag : string
        The name of the tag of the flux mesh. Flux values for different energy
        groups are assumed to be represented as vector tags.
    fluxin : string
        The name of the ALARA fluxin file to be output.
    reverse : bool
        If true, fluxes will be printed in the reverse order as they appear in
        the flux vector tagged on the mesh.
    """
    tag_flux = flux_mesh.mesh.getTagHandle(flux_tag)

    if flux_mesh.structured:
        ves = flux_mesh.structured_iterate_hex("xyz")
    else:
        ves = flux.mesh.mesh.iterate(iBase.Type.region, iMesh.Toplogy.all)

    # find number of e_groups
    e_groups = flux_mesh.mesh.getTagHandle(flux_tag)[list(
        flux_mesh.mesh.iterate(iBase.Type.region,
                               iMesh.Topology.all))[0]]
    e_groups = np.atleast_1d(e_groups)
    num_e_groups = len(e_groups)

    # Establish for loop bounds based on if forward or backward printing
    # is requested
    if not reverse:
        start = 0
        stop = num_e_groups
        direction = 1
    else:
        start = num_e_groups - 1
        stop = -1
        direction = -1

    output = ""
    for ve in ves:
        # print flux data to file
        count = 0
        flux_data = np.atleast_1d(tag_flux[ve])
        for i in range(start, stop, direction):
            output += "{:.6E} ".format(flux_data[i])
            # fluxin formatting: create a new line after every 6th entry
            count += 1
            if count % 6 == 0:
                output += "\n"

        output += "\n\n"

    with open(fluxin, "w") as f:
        f.write(output)


def photon_source_to_hdf5(filename, chunkshape=(10000,)):
    """Converts a plaintext photon source file to an HDF5 version for
    quick later use.

    Parameters
    ----------
    filename : str
        The path to the file
    chunkshape : tuple of int
        A 1D tuple of the HDF5 chunkshape.

    Output
    ------
    A single HDF5 file named <filename>.h5 containing the table headings:

    idx : int
        The volume element index assuming the volume elements appear in xyz
        order (z changing fastest) within the photon source file in the case of
        a structured mesh or imesh.iterate() order for an unstructured mesh.
    nuc : str
        The nuclide name as it appears in the photon source file.
    time : str
        The decay time as it appears in the photon source file.
    phtn_src : 1D array of floats
        Contains the photon source density for each energy group.
    """
    f = open(filename, 'r')
    header = f.readline().strip().split('\t')
    f.seek(0)
    G = len(header) - 2

    dt = np.dtype([
        ('idx', np.int64),
        ('nuc', 'S6'),
        ('time', 'S20'),
        ('phtn_src', np.float64, G),
        ])

    filters = tb.Filters(complevel=1, complib='zlib')
    h5f = tb.openFile(filename + '.h5', 'w', filters=filters)
    tab = h5f.createTable('/', 'data', dt, chunkshape=chunkshape)

    chunksize = chunkshape[0]
    rows = np.empty(chunksize, dtype=dt)
    idx = 0
    old = ""
    for i, line in enumerate(f, 1):
        ls = line.strip().split('\t')

        # Keep track of the idx by delimiting by the last TOTAL line in a
        # volume element.
        if ls[0] != 'TOTAL' and old == 'TOTAL':
            idx += 1

        j = (i-1) % chunksize
        rows[j] = (idx, ls[0].strip(), ls[1].strip(),
                   np.array(ls[2:], dtype=np.float64))
        # Save the nuclide in order to keep track of idx
        old = ls[0]

        if i % chunksize == 0:
            tab.append(rows)
            rows = np.empty(chunksize, dtype=dt)

    if i % chunksize != 0:
        tab.append(rows[:j+1])

    h5f.close()
    f.close()


def photon_source_hdf5_to_mesh(mesh, filename, tags):
    """This function reads in an hdf5 file produced by photon_source_to_hdf5
    and tags the requested data to the mesh of a PyNE Mesh object. Any
    combinations of nuclides and decay times are allowed. The photon source
    file is assumed to be in xyz order (z changes fastest) if a stuctured mesh
    is supplied and imesh.iterate() order if an unstructured mesh is supplied.

    Parameters
    ----------
    mesh : PyNE Mesh
       The object containing the imesh instance to be tagged.
    filename : str
        The path of the hdf5 version of the photon source file.
    tags: dict
        A dictionary were the keys are tuples with two values. The first is a
        string denoting an nuclide in any form that is understood by
        pyne.nucname (e.g. '1001', 'U-235', '242Am') or 'TOTAL' for all
        nuclides. The second is a string denoting the decay time as it appears
        in the file (e.g. 'shutdown', '1 h' '3 d'). The values of the
        dictionary are the requested tag names for the combination of nuclide
        and decay time. For example if one wanted tags for the photon source
        densities from U235 at shutdown and from all nuclides at 1 hour, the
        dictionary could be:

        tags = {('U-235', 'shutdown') : 'tag1', ('TOTAL', '1 h') : 'tag2'}
    """
    # find number of energy groups
    with tb.openFile(filename) as h5f:
        num_e_groups = len(h5f.root.data[0][3])

    # create a dict of tag handles for all keys of the tags dict
    tag_handles = {}
    for tag_name in tags.values():
        tag_handles[tag_name] = \
            mesh.mesh.createTag(tag_name, num_e_groups, float)

    # iterate through each requested nuclide/dectay time
    for cond in tags.keys():
        with tb.openFile(filename) as h5f:
            # Convert nuclide to the form found in the ALARA phtn_src
            # file, which is similar to the Serpent form. Note this form is
            # different from the ALARA input nuclide form found in nucname.
            if cond[0] != "TOTAL":
                nuc = serpent(cond[0]).lower()
            else:
                nuc = "TOTAL"
            # create of array of rows that match the nuclide/decay criteria
            matched_data = h5f.root.data.readWhere(
                "(nuc == '{0}') & (time == '{1}')".format(nuc, cond[1]))

        if mesh.structured:
            ves = mesh.structured_iterate_hex("xyz")
        else:
            ves = mesh.imesh.iterate(iBase.Type.region, iMesh.Topology.all)

        idx = 0
        for i, ve in enumerate(list(ves)):
            if matched_data[idx][0] == i:
                tag_handles[tags[cond]][ve] = matched_data[idx][3]
                idx += 1
            else:
                tag_handles[tags[cond]][ve] = [0] * num_e_groups


def mesh_to_geom(mesh, geom_file, matlib_file):
    """This function reads the materials of a PyNE mesh object and prints the
    geometry and materials portion of an ALARA input file, as well as a
    corresponding matlib file. If the mesh is structured, xyz ordering is used
    (z changing fastest). If the mesh is unstructured iMesh.iterate order is
    used. 

    Parameters
    ----------
    mesh : PyNE Mesh object
        The Mesh object containing the materials to be printed.
    geom_file : str
        The name of the file to print the geometry and material blocks.
    matlib_file : str
        The name of the file to print the matlib.
    """
    # Create geometry information header. Note that the shape of the geometry
    # (rectangular) is actually inconsequential to the ALARA calculation so
    # unstructured meshes are not adversely affected. 
    geometry = "geometry rectangular\n\n"

    # Create three strings in order to create all ALARA input blocks in a
    # single mesh iteration.
    volume = "volume\n" # volume input block
    mat_loading = "mat_loading\n" # material loading input block
    mixture = "" # mixture blocks
    matlib = "" # ALARA material library string

    if mesh.structured:
        ves = mesh.structured_iterate_hex('xyz')
    else:
        ves = mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)

    for i, ve in enumerate(ves):
        volume += "    {0: 1.6E}    zone_{1}\n".format(mesh.elem_volume(ve), i)
        mat_loading += "    zone_{0}    mix_{0}\n".format(i)
        matlib += "mat_{0}    {1: 1.6E}    {2}\n".format(i, mesh.density[i], 
                                                         len(mesh.comp[i]))
        mixture += ("mixture mix_{0}\n"
                    "    material mat_{0} 1 1\nend\n\n".format(i))

        for nuc, comp in mesh.comp[i].iteritems():
            matlib += "{0}    {1: 1.6E}    {2}\n".format(alara(nuc), comp, 
                                                         znum(nuc))
        matlib += "\n"

    volume += "end\n\n"
    mat_loading += "end\n\n"

    with open(geom_file, 'w') as f:
        f.write(geometry + volume + mat_loading + mixture)
    
    with open(matlib_file, 'w') as f:
        f.write(matlib)

def num_density_to_mesh(lines, time, m):
    """This function reads ALARA output containing number density information
    and creates material objects which are then added to a supplied PyNE Mesh
    object. The volumes within ALARA are assummed to appear in the same order as
    the idx on the Mesh object. 

    Parameters
    ----------
    lines : list or str
        ALARA output from ALARA run with 'number_density' in the 'output' block
        of the input file. Lines can either be a filename or the equivalent to
        calling readlines() on an ALARA output file. If reading in ALARA output
        from stdout, call split('\n') before passing it in as the lines
        parameter.
    time : str
        The decay time for which number densities are requested (e.g. '1 h',
        'shutdown', etc.)
    m : PyNE Mesh
        Mesh object for which mats will be applied to.
    """
    if isinstance(lines, basestring):
        with open(lines) as f:
            lines = f.readlines()
    elif not isinstance(lines, collections.Sequence):
        raise TypeError("Lines argument not a file or sequence.")
    # Advance file to number density portion.
    header = 'Number Density [atoms/cm3]\n'
    line = ""
    while line != header:
        line = lines.pop(0)

    # Get decay time index from next line (the column the decay time answers
    # appear in.
    line_strs = lines.pop(0).replace('\t', '  ')
    time_index = [s.strip() for s in line_strs.split('  ') 
                  if s.strip()].index(time)

    # Create a dict of mats for the mesh.
    mats = {}
    count = 0
    # Read through file until enough material objects are create to fill mesh.
    while count != len(m):
        # Pop lines to the start of the next material.
        while lines.pop(0)[0] != '=':
            pass

        # Create a new material object and add to mats dict.
        line = lines.pop(0)
        nucvec = {}
        density = 0.0
        # Read lines until '=' delimiter at the end of a material.
        while line[0] != '=':
            nuc = line.split()[0]
            n = float(line.split()[time_index])
            if n != 0.0:
                nucvec[nuc] = n
                density += n * anum(nuc)/N_A

            line = lines.pop(0)
        mat = from_atom_frac(nucvec, density=density, mass=0)
        mats[count] = mat
        count += 1

    m.mats = mats


def irradiation_blocks(material_lib, element_lib, data_library, cooling, 
                       flux_file, irr_time, output = "number_density",
                       truncation=1E-12, impurity = (5E-6, 1E-3), 
                       dump_file = "dump_file"):
    """irradiation_blocks(material_lib, element_lib, data_library, cooling, 
                       flux_file, irr_time, output = "number_density",
                       truncation=1E-12, impurity = (5E-6, 1E-3), 
                       dump_file = "dump_file")

    This function returns a string of the irradation-related input blocks. This 
    function is meant to be used with files created by the mesh_to_geom 
    function, in order to append the remaining input blocks to form a complete
    ALARA input file. Only the simplest irradiation schedule is supported: a 
    single pulse of time <irr_time>. The notation in this function is consistent 
    with the ALARA users' guide, found at:

    http://alara.engr.wisc.edu/users.guide.html/

    Parameters
    ----------
    material_lib : str
        Path to material library.
    element_lib : str
        Path to element library.
    data_library : str
        The data_library card (see ALARA user's guide).
    cooling : str or iterable of str
        Cooling times for which output is requested. Given in ALARA form (e.g.
        "1 h", "0.5 y"). Note that "shutdown" is always implicitly included.
    flux_file : str
        Path to the "fluxin" file.
    irr_time : str
        The duration of the single pulse irradiation. Given in the ALARA form
        (e.g. "1 h", "0.5 y").
    output : str or iterable of str, optional.
        The requested output blocks (see ALARA users' guide).
    truncation : float, optional
        The chain truncation value (see ALARA users' guide).
    impurity : tuple of two floats, optional
       The impurity parameters (see ALARA users' guide).
    dump_file: str, optional
       Path to the dump file.
  
    Returns
    -------
    s : str
        Irradition-related ALARA input blocks.
    """

    s = ""

    # Material, element, and data_library blocks
    s += "material_lib {0}\n".format(material_lib)
    s += "element_lib {0}\n".format(element_lib)
    s += "data_library {0}\n\n".format(data_library)

    # Cooling times
    s += "cooling\n"
    if isinstance(cooling, collections.Iterable) and not isinstance(cooling, basestring):
        for c in cooling:
            s += "    {0}\n".format(c)
    else:
        s += "    {0}\n".format(cooling)

    s += "end\n\n"

    # Flux block
    s += "flux flux_1 {0} 1.0 0 default\n".format(flux_file)

    # Flux schedule
    s += ("schedule simple_schedule\n"
         "    {0} flux_1 pulse_once 0 s\nend\n\n".format(irr_time))

    s += "pulsehistory pulse_once\n    1 0.0 s\nend\n\n"
 
    # Output block
    s += "output zone\n    units Ci cm3\n"
    if isinstance(output, collections.Iterable) and not isinstance(output, basestring):
        for out in output:
            s += "    {0}\n".format(out)
    else:
        s += "    {0}\n".format(output)

    s += "end\n\n"

    # Other parameters
    s += "truncation {0}\n".format(truncation)
    s += "impurity {0} {1}\n".format(impurity[0], impurity[1])
    s += "dump_file {0}\n".format(dump_file)

    return s
