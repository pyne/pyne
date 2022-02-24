"""This module contains functions relevant to the ALARA activation code and the Chebyshev Rational Approximation Method
"""
from __future__ import print_function
from pyne.xs.data_source import SimpleDataSource
from pyne.data import N_A, decay_const, decay_children, branch_ratio
from pyne.nucname import serpent, alara, znum, anum
from pyne import nucname
from pyne.material import Material, from_atom_frac
from pyne.mesh import Mesh, MeshError, HAVE_PYMOAB
import os
import collections

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections
from warnings import warn
from pyne.utils import QA_warn, to_sec, str_to_unicode, IfBar
import numpy as np
import tables as tb
from io import open
import re

QA_warn(__name__)

try:
    basestring
except NameError:
    basestring = str

if HAVE_PYMOAB:
    from pyne.mesh import mesh_iterate
else:
    warn(
        "The PyMOAB optional dependency could not be imported. "
        "Some aspects of the mesh module may be incomplete.",
        ImportWarning,
    )

response_strings = {
    "decay_heat": "Total Decay Heat",
    "specific_activity": "Specific Activity",
    "alpha_heat": "Alpha Decay Heat",
    "beta_heat": "Beta Decay Heat",
    "gamma_heat": "Gamma Decay Heat",
    "wdr": "WDR/Clearance index",
    "photon_source": "Photon Source Distribution",
}


def mesh_to_fluxin(
    flux_mesh,
    flux_tag,
    fluxin="fluxin.out",
    reverse=False,
    sub_voxel=False,
    cell_fracs=None,
    cell_mats=None,
    print_progress=True,
):
    """This function creates an ALARA fluxin file from fluxes tagged on a PyNE
    Mesh object. Fluxes are printed in the order of the flux_mesh.__iter__().

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
    sub_voxel: bool, optional
        If true, sub-voxel r2s work flow will be sued. Flux of a voxel will
        be duplicated c times. Where c is the cell numbers of that voxel.
    cell_fracs : structured array, optional
        The output from dagmc.discretize_geom(). A sorted, one dimensional
        array, each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.

        The array must be sorted with respect to both idx and cell, with
        cell changing fastest.
    cell_mats : dict, optional
        Maps geometry cell numbers to PyNE Material objects.
        The cell_fracs and cell_mats are used only when sub_voxel=True.
        If sub_voxel=False, neither cell_fracs nor cell_mats will be used.
    print_progress : bool
        Print the progress with progress bar. Set to False to turn off progress.
    """

    tag_flux = flux_mesh.get_tag(flux_tag)

    # find number of e_groups
    e_groups = tag_flux[list(mesh_iterate(flux_mesh.mesh))[0]]
    e_groups = np.atleast_1d(e_groups)
    num_e_groups = len(e_groups)

    # write flux data block by block
    with open(fluxin, "w") as f:
        if not sub_voxel:
            bar = IfBar(
                "Writing alara fluxin",
                max=len(flux_mesh),
                suffix="%(percent).1f%% - %(eta)ds",
                show=print_progress,
            )
            for i, mat, ve in flux_mesh:
                f.write(_output_flux_block(ve, tag_flux, reverse))
                bar.next()
            bar.finish()
        else:
            ves = list(flux_mesh.iter_ve())
            bar = IfBar(
                "Writing alara fluxin",
                max=len(cell_fracs),
                suffix="%(percent).1f%% - %(eta)ds",
                show=print_progress,
            )
            for i, row in enumerate(cell_fracs):
                if len(cell_mats[row["cell"]].comp) != 0:
                    f.write(_output_flux_block(ves[row["idx"]], tag_flux, reverse))
                bar.next()
            bar.finish()


def photon_source_to_hdf5(filename, nucs="all", chunkshape=(10000,)):
    """Converts a plaintext photon source file to an HDF5 version for
    quick later use.

    This function produces a single HDF5 file named <filename>.h5 containing the
    table headings:

        idx : int
            The volume element index assuming the volume elements appear in xyz
            order (z changing fastest) within the photon source file in the case of
            a structured mesh or mesh.mesh_iterate() order for an unstructured mesh.
        nuc : str
            The nuclide name as it appears in the photon source file.
        time : str
            The decay time as it appears in the photon source file.
        phtn_src : 1D array of floats
            Contains the photon source density for each energy group.

    Parameters
    ----------
    filename : str
        The path to the file for read
    nucs : str
        Nuclides need to write into h5 file. For example:
            - 'all': default value. Write the information of all nuclides to h5.
            - 'total': used for r2s. Only write TOTAL value to h5.
    chunkshape : tuple of int
        A 1D tuple of the HDF5 chunkshape.

    """
    f = open(filename, "r")
    header = f.readline().strip().split("\t")
    f.seek(0)
    G = len(header) - 2

    phtn_dtype = _make_response_dtype("phtn_src", data_length=G)

    filters = tb.Filters(complevel=1, complib="zlib")
    # set the default output h5_filename
    h5f = tb.open_file(filename + ".h5", "w", filters=filters)
    tab = h5f.create_table("/", "data", phtn_dtype, chunkshape=chunkshape)

    chunksize = chunkshape[0]
    rows = np.empty(chunksize, dtype=phtn_dtype)
    idx = 0
    old = ""
    row_count = 0
    for i, line in enumerate(f, 1):
        tokens = line.strip().split("\t")

        # Keep track of the idx by delimiting by the last TOTAL line in a
        # volume element.
        if tokens[0] != "TOTAL" and old == "TOTAL":
            idx += 1

        if nucs.lower() == "all":
            row_count += 1
            j = (row_count - 1) % chunksize
            rows[j] = (
                idx,
                tokens[0].strip(),
                tokens[1].strip(),
                np.array(tokens[2:], dtype=np.float64),
            )
        elif nucs.lower() == "total":
            if tokens[0] == "TOTAL":
                row_count += 1
                j = (row_count - 1) % chunksize
                rows[j] = (
                    idx,
                    tokens[0].strip(),
                    tokens[1].strip(),
                    np.array(tokens[2:], dtype=np.float64),
                )
        else:
            h5f.close()
            f.close()
            raise ValueError("Nucs option {0} not support!".format(nucs))

        # Save the nuclide in order to keep track of idx
        old = tokens[0]

        if (row_count > 0) and (row_count % chunksize == 0):
            tab.append(rows)
            rows = np.empty(chunksize, dtype=phtn_dtype)
            row_count = 0

    if row_count % chunksize != 0:
        tab.append(rows[: j + 1])

    h5f.close()
    f.close()


def response_to_hdf5(filename, response, chunkshape=(10000,)):
    """Converts a plaintext output.txt file to an HDF5 version for
    quick later use.

    This function produces a single HDF5 file named <filename>.h5 containing the
    table headings:

        idx : int
            The volume element index assuming the volume elements appear in xyz
            order (z changing fastest) within the photon source file in the case of
            a structured mesh or mesh.mesh_iterate() order for an unstructured mesh.
        nuc : str
            The nuclide name as it appears in the output file.
        time : str
            The decay time as it appears in the output file.
        response : float
            Supported response:
                - decay_heat [W/cm3]
                - specific_activity [Bq/cm3]
                - alpha_heat [W/cm3]
                - beta_heat [W/cm3]
                - gamma_heat [W/cm3]
                - wdr
                - photon_source

    Parameters
    ----------
    filename : str
        The path to the file output.txt
    response : str
        Key word of the response
    chunkshape : tuple of int
        A 1D tuple of the HDF5 chunkshape.
    """
    f = open(filename, "r")
    f.seek(0)

    response_dtype = _make_response_dtype(response)

    filters = tb.Filters(complevel=1, complib="zlib")
    h5_filename = os.path.join(os.path.dirname(filename), "".join([response, ".h5"]))
    h5f = tb.open_file(h5_filename, "w", filters=filters)
    tab = h5f.create_table("/", "data", response_dtype, chunkshape=chunkshape)

    chunksize = chunkshape[0]
    rows = np.empty(chunksize, dtype=response_dtype)
    zone_idx = 0
    count = 1
    decay_times = []
    zone_start_state = 1
    response_start_state = 2
    state = None
    for line in f:
        # terminate condition
        if state == response_start_state and ("Totals for all zones." in line):
            break
        # get response string
        if state == zone_start_state and (response_strings[response] in line):
            state = response_start_state
        # get decay times
        elif (
            state == response_start_state
            and len(decay_times) == 0
            and ("isotope\t shutdown" in line)
        ):
            decay_times = read_decay_times(line)
        # get zone idx
        elif "Zone #" in line:
            zone_idx = _get_zone_idx(line)
            if zone_idx == 0:
                state = zone_start_state
        # skip lines if we haven't started the response or
        # the lines don't contain wanted data
        elif state == response_start_state and _is_data(line):

            tokens = line.strip().split()
            # put data into table
            # format of each row: zone_idx, nuc, time, decay_heat
            nuc = tokens[0].strip()
            if nuc.lower() == "total":
                nuc = nuc.upper()
            for dt, response_value in zip(decay_times, tokens[1:]):
                j = (count - 1) % chunksize
                rows[j] = (zone_idx, nuc, dt, response_value)
                if count % chunksize == 0:
                    tab.append(rows)
                    rows = np.empty(chunksize, dtype=response_dtype)
                count += 1

    if count % chunksize != 0:
        tab.append(rows[: j + 1])

    # close the file
    h5f.close()
    f.close()


def photon_source_hdf5_to_mesh(mesh, filename, tags, sub_voxel=False, cell_mats=None):
    """This function reads in an hdf5 file produced by photon_source_to_hdf5
    and tags the requested data to the mesh of a PyNE Mesh object. Any
    combinations of nuclides and decay times are allowed. The photon source
    file is assumed to be in mesh.__iter__() order

    Parameters
    ----------
    mesh : PyNE Mesh
       The object containing the PyMOAB instance to be tagged.
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
    sub_voxel: bool, optional
        If the sub_voxel is True, then the sub-voxel r2s will be used.
        Then the photon_source will be interpreted as sub-voxel photon source.
    cell_mats : dict, optional
        cell_mats is required when sub_voxel is True.
        Maps geometry cell numbers to PyNE Material objects.
    """

    # find number of energy groups
    with tb.open_file(filename) as h5f:
        num_e_groups = len(h5f.root.data[0][3])
    max_num_cells = 1
    ve0 = next(mesh.iter_ve())
    if sub_voxel:
        num_vol_elements = len(mesh)
        subvoxel_array = _get_subvoxel_array(mesh, cell_mats)
        # get max_num_cells
        max_num_cells = len(np.atleast_1d(mesh.cell_number[ve0]))

    # create a dict of tag handles for all keys of the tags dict
    tag_handles = {}
    tag_size = num_e_groups * max_num_cells
    for tag_name in tags.values():

        mesh.tag(
            tag_name,
            np.zeros(tag_size, dtype=float),
            "nat_mesh",
            size=tag_size,
            dtype=float,
        )
        tag_handles[tag_name] = mesh.get_tag(tag_name)

    decay_times = _read_h5_dt(filename)

    # iterate through each requested nuclide/dectay time
    for cond in tags.keys():
        with tb.open_file(filename) as h5f:
            # Convert nuclide to the form found in the ALARA phtn_src
            # file, which is similar to the Serpent form. Note this form is
            # different from the ALARA input nuclide form found in nucname.
            if cond[0] != "TOTAL":
                nuc = serpent(cond[0]).lower()
            else:
                nuc = "TOTAL"

            # time match, convert string mathch to float mathch
            dt = _find_dt(cond[1], decay_times)
            # create of array of rows that match the nuclide/decay criteria
            matched_data = h5f.root.data.read_where(
                "(nuc == '{0}') & (time == '{1}')".format(nuc, dt)
            )

        if not sub_voxel:
            idx = 0
            for i, _, ve in mesh:
                if matched_data[idx][0] == i:
                    tag_handles[tags[cond]][ve] = matched_data[idx][3]
                    idx += 1
                else:
                    tag_handles[tags[cond]][ve] = [0] * num_e_groups
        else:
            temp_mesh_data = np.empty(
                shape=(num_vol_elements, max_num_cells, num_e_groups), dtype=float
            )
            temp_mesh_data.fill(0.0)
            for sve, subvoxel in enumerate(subvoxel_array):
                temp_mesh_data[subvoxel["idx"], subvoxel["scid"], :] = matched_data[
                    sve
                ][3][:]
            for i, _, ve in mesh:
                tag_handles[tags[cond]][ve] = temp_mesh_data[i, :].reshape(
                    max_num_cells * num_e_groups
                )


def response_hdf5_to_mesh(mesh, filename, tags, response):
    """This function reads in an hdf5 file produced by response_to_hdf5
    and tags the requested data to the mesh of a PyNE Mesh object. Any
    combinations of nuclides and decay times are allowed. The photon source
    file is assumed to be in mesh.__iter__() order

    Parameters
    ----------
    mesh : PyNE Mesh
       The object containing the PyMOAB instance to be tagged.
    filename : str
        The path of the hdf5 version of the response file.
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
    response : str
        The keyword of the response. Supported responses:
            - decay_heat
            - specific_activity
            - alpha_heat
            - beta_heat
            - gamma_heat
            - wdr
            - photon_source
    """

    # create a dict of tag handles for all keys of the tags dict
    tag_handles = {}
    for tag_name in tags.values():

        mesh.tag(tag_name, np.zeros(1, dtype=float), "nat_mesh", size=1, dtype=float)
        tag_handles[tag_name] = mesh.get_tag(tag_name)

    decay_times = _read_h5_dt(filename)

    # iterate through each requested nuclide/dectay time
    for cond in tags.keys():
        with tb.open_file(filename) as h5f:
            # Convert nuclide to the form found in the ALARA response file
            # file, which is similar to the Serpent form. Note this form is
            # different from the ALARA input nuclide form found in nucname.
            if cond[0] != "TOTAL":
                nuc = serpent(cond[0]).lower()
            else:
                nuc = "TOTAL"

            # time match, convert string mathch to float mathch
            dt = _find_dt(cond[1], decay_times)
            # create of array of rows that match the nuclide/decay criteria
            matched_data = h5f.root.data.read_where(
                "(nuc == '{0}') & (time == '{1}')".format(nuc, dt)
            )

        idx = 0
        # index, mat, volume element
        for i, _, ve in mesh:
            if matched_data[idx][0] == i:
                tag_handles[tags[cond]][ve] = matched_data[idx][3]
                idx += 1
            else:
                tag_handles[tags[cond]][ve] = [0]


def record_to_geom(
    mesh, cell_fracs, cell_mats, geom_file, matlib_file, sig_figs=6, sub_voxel=False
):
    """This function preforms the same task as alara.mesh_to_geom, except the
    geometry is on the basis of the stuctured array output of
    dagmc.discretize_geom rather than a PyNE material object with materials.
    This allows for more efficient ALARA runs by minimizing the number of
    materials in the ALARA matlib. This is done by treating mixtures that are
    equal up to <sig_figs> digits to be the same mixture within ALARA.

    Parameters
    ----------
    mesh : PyNE Mesh object
        The Mesh object for which the geometry is discretized.
    cell_fracs : structured array
        The output from dagmc.discretize_geom(). A sorted, one dimensional
        array, each entry containing the following fields:

            :idx: int
                The volume element index.
            :cell: int
                The geometry cell number.
            :vol_frac: float
                The volume fraction of the cell withing the mesh ve.
            :rel_error: float
                The relative error associated with the volume fraction.

    cell_mats : dict
        Maps geometry cell numbers to PyNE Material objects. Each PyNE material
        object must have 'name' specified in Material.metadata.
    geom_file : str
        The name of the file to print the geometry and material blocks.
    matlib_file : str
        The name of the file to print the matlib.
    sig_figs : int
        The number of significant figures that two mixtures must have in common
        to be treated as the same mixture within ALARA.
    sub_voxel : bool
        If sub_voxel is True, the sub-voxel r2s will be used.
    """
    # Create geometry information header. Note that the shape of the geometry
    # (rectangular) is actually inconsequential to the ALARA calculation so
    # unstructured meshes are not adversely affected.
    geometry = "geometry rectangular\n\n"

    # Create three strings in order to create all ALARA input blocks in a
    # single mesh iteration.
    volume = "volume\n"  # volume input block
    mat_loading = "mat_loading\n"  # material loading input block
    mixture = ""  # mixture blocks

    unique_mixtures = []
    if not sub_voxel:
        for i, mat, ve in mesh:
            volume += "    {0: 1.6E}    zone_{1}\n".format(mesh.elem_volume(ve), i)

            ve_mixture = {}
            for row in cell_fracs[cell_fracs["idx"] == i]:
                cell_mat = cell_mats[row["cell"]]
                name = cell_mat.metadata["name"]
                if _is_void(name):
                    name = "mat_void"
                if name not in ve_mixture.keys():
                    ve_mixture[name] = np.round(row["vol_frac"], sig_figs)
                else:
                    ve_mixture[name] += np.round(row["vol_frac"], sig_figs)

            if ve_mixture not in unique_mixtures:
                unique_mixtures.append(ve_mixture)
                mixture += "mixture mix_{0}\n".format(unique_mixtures.index(ve_mixture))
                for key, value in ve_mixture.items():
                    mixture += "    material {0} 1 {1}\n".format(key, value)

                mixture += "end\n\n"

            mat_loading += "    zone_{0}    mix_{1}\n".format(
                i, unique_mixtures.index(ve_mixture)
            )
    else:
        ves = list(mesh.iter_ve())
        sve_count = 0
        for row in cell_fracs:
            if len(cell_mats[row["cell"]].comp) != 0:
                volume += "    {0: 1.6E}    zone_{1}\n".format(
                    mesh.elem_volume(ves[row["idx"]]) * row["vol_frac"], sve_count
                )
                cell_mat = cell_mats[row["cell"]]
                name = cell_mat.metadata["name"]
                if name not in unique_mixtures:
                    unique_mixtures.append(name)
                    mixture += "mixture {0}\n".format(name)
                    mixture += "    material {0} 1 1\n".format(name)
                    mixture += "end\n\n"
                mat_loading += "    zone_{0}    {1}\n".format(sve_count, name)
                sve_count += 1

    volume += "end\n\n"
    mat_loading += "end\n\n"

    with open(geom_file, "w") as f:
        f.write(geometry + volume + mat_loading + mixture)

    matlib = ""  # ALARA material library string

    printed_mats = []
    print_void = False
    for mat in cell_mats.values():
        name = mat.metadata["name"]
        if _is_void(name):
            print_void = True
            continue
        if name not in printed_mats:
            printed_mats.append(name)
            matlib += "{0}    {1: 1.6E}    {2}\n".format(
                name, mat.density, len(mat.comp)
            )
            for nuc, comp in mat.comp.items():
                matlib += "{0}    {1: 1.6E}    {2}\n".format(
                    alara(nuc), comp * 100.0, znum(nuc)
                )
            matlib += "\n"

    if print_void:
        matlib += "# void material\nmat_void 0.0 1\nhe 1 2\n"

    with open(matlib_file, "w") as f:
        f.write(matlib)


def _is_void(name):
    """Private function for determining if a material name specifies void."""
    lname = name.lower()
    return "vacuum" in lname or "void" in lname or "graveyard" in lname


def mesh_to_geom(mesh, geom_file, matlib_file):
    """This function reads the materials of a PyNE mesh object and prints the
    geometry and materials portion of an ALARA input file, as well as a
    corresponding matlib file. If the mesh is structured, xyz ordering is used
    (z changing fastest). If the mesh is unstructured the mesh_iterate order is
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
    volume = "volume\n"  # volume input block
    mat_loading = "mat_loading\n"  # material loading input block
    mixture = ""  # mixture blocks
    matlib = ""  # ALARA material library string

    for i, mat, ve in mesh:
        volume += "    {0: 1.6E}    zone_{1}\n".format(mesh.elem_volume(ve), i)
        mat_loading += "    zone_{0}    mix_{0}\n".format(i)
        matlib += "mat_{0}    {1: 1.6E}    {2}\n".format(
            i, mesh.density[i], len(mesh.comp[i])
        )
        mixture += "mixture mix_{0}\n" "    material mat_{0} 1 1\nend\n\n".format(i)

        for nuc, comp in mesh.comp[i].items():
            matlib += "{0}    {1: 1.6E}    {2}\n".format(
                alara(nuc), comp * 100.0, znum(nuc)
            )
        matlib += "\n"

    volume += "end\n\n"
    mat_loading += "end\n\n"

    with open(geom_file, "w") as f:
        f.write(geometry + volume + mat_loading + mixture)

    with open(matlib_file, "w") as f:
        f.write(matlib)


def num_density_to_mesh(lines, time, m):
    """num_density_to_mesh(lines, time, m)
    This function reads ALARA output containing number density information and
    creates material objects which are then added to a supplied PyNE Mesh object.
    The volumes within ALARA are assummed to appear in the same order as the
    idx on the Mesh object.
    All the strings in this function is text (unicode).

    Parameters
    ----------
    lines : list or str
        ALARA output from ALARA run with 'number_density' in the 'output' block
        of the input file. Lines can either be a filename or the equivalent to
        calling readlines() on an ALARA output file. If reading in ALARA output
        from stdout, call split('\n') before passing it in as the lines parameter.
    time : str
        The decay time for which number densities are requested (e.g. '1 h',
        'shutdown', etc.)
    m : PyNE Mesh
        Mesh object for which mats will be applied to.
    """
    if isinstance(lines, basestring):
        with open(lines) as f:
            lines = f.readlines()
    elif not isinstance(lines, collectionsAbc.Sequence):
        raise TypeError("Lines argument not a file or sequence.")
    # Advance file to number density portion.
    header = "Number Density [atoms/cm3]"
    line = ""
    while line.rstrip() != header:
        line = lines.pop(0)

    # Get decay time index from next line (the column the decay time answers
    # appear in.
    line_strs = lines.pop(0).replace("\t", "  ")
    time_index = [s.strip() for s in line_strs.split("  ") if s.strip()].index(time)

    # Create a dict of mats for the mesh.
    mats = {}
    count = 0
    # Read through file until enough material objects are create to fill mesh.
    while count != len(m):
        # Pop lines to the start of the next material.
        while (lines.pop(0) + " ")[0] != "=":
            pass

        # Create a new material object and add to mats dict.
        line = lines.pop(0)
        nucvec = {}
        density = 0.0
        # Read lines until '=' delimiter at the end of a material.
        while line[0] != "=":
            nuc = line.split()[0]
            n = float(line.split()[time_index])
            if n != 0.0:
                nucvec[nuc] = n
                density += n * anum(nuc) / N_A

            line = lines.pop(0)
        mat = from_atom_frac(nucvec, density=density, mass=0)
        mats[count] = mat
        count += 1

    m.mats = mats


def irradiation_blocks(
    material_lib,
    element_lib,
    data_library,
    cooling,
    flux_file,
    irr_time,
    output="number_density",
    truncation=1e-12,
    impurity=(5e-6, 1e-3),
    dump_file="dump_file",
):
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

    http://svalinn.github.io/ALARA/usersguide/index.html

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
    if isinstance(cooling, collectionsAbc.Iterable) and not isinstance(
        cooling, basestring
    ):
        for c in cooling:
            s += "    {0}\n".format(c)
    else:
        s += "    {0}\n".format(cooling)

    s += "end\n\n"

    # Flux block
    s += "flux flux_1 {0} 1.0 0 default\n".format(flux_file)

    # Flux schedule
    s += "schedule simple_schedule\n" "    {0} flux_1 pulse_once 0 s\nend\n\n".format(
        irr_time
    )

    s += "pulsehistory pulse_once\n    1 0.0 s\nend\n\n"

    # Output block
    s += "output zone\n    units Ci cm3\n"
    if isinstance(output, collectionsAbc.Iterable) and not isinstance(
        output, basestring
    ):
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


def phtn_src_energy_bounds(input_file):
    """Reads an ALARA input file and extracts the energy bounds from the
    photon_source block.

    Parameters
    ----------
    input_file : str
         The ALARA input file name, which must contain a photon_source block.

    Returns
    -------
    e_bounds : list of floats
    The lower and upper energy bounds for the photon_source discretization. Unit: eV.
    """
    phtn_src_lines = ""
    with open(input_file, "r") as f:
        line = f.readline()
        while not (" photon_source " in line and line.strip()[0] != "#"):
            line = f.readline()
        num_groups = float(line.split()[3])
        upper_bounds = [float(x) for x in line.split()[4:]]
        while len(upper_bounds) < num_groups:
            line = f.readline()
            upper_bounds += [
                float(x) for x in line.split("#")[0].split("end")[0].split()
            ]
    e_bounds = [0.0] + upper_bounds
    return e_bounds


def _build_matrix(N):
    """This function  builds burnup matrix, A. Decay only."""

    A = np.zeros((len(N), len(N)))

    # convert N to id form
    N_id = []
    for i in range(len(N)):
        if isinstance(N[i], str):
            ID = nucname.id(N[i])
        else:
            ID = N[i]
        N_id.append(ID)

    sds = SimpleDataSource()

    # Decay
    for i in range(len(N)):
        A[i, i] -= decay_const(N_id[i])

        # Find decay parents
        for k in range(len(N)):
            if N_id[i] in decay_children(N_id[k]):
                A[i, k] += branch_ratio(N_id[k], N_id[i]) * decay_const(N_id[k])
    return A


def _rat_apprx_14(A, t, n_0):
    """CRAM of order 14

    Parameters
    ---------
    A : numpy array
        Burnup matrix
    t : float
        Time step
    n_0: numpy array
        Inital composition vector
    """

    theta = np.array(
        [
            -8.8977731864688888199 + 16.630982619902085304j,
            -3.7032750494234480603 + 13.656371871483268171j,
            -0.2087586382501301251 + 10.991260561901260913j,
            3.9933697105785685194 + 6.0048316422350373178j,
            5.0893450605806245066 + 3.5888240290270065102j,
            5.6231425727459771248 + 1.1940690463439669766j,
            2.2697838292311127097 + 8.4617379730402214019j,
        ]
    )

    alpha = np.array(
        [
            -0.000071542880635890672853 + 0.00014361043349541300111j,
            0.0094390253107361688779 - 0.01784791958483017511j,
            -0.37636003878226968717 + 0.33518347029450104214j,
            -23.498232091082701191 - 5.8083591297142074004j,
            46.933274488831293047 + 45.643649768827760791j,
            -27.875161940145646468 - 102.14733999056451434j,
            4.8071120988325088907 - 1.3209793837428723881j,
        ]
    )

    alpha_0 = np.array([1.8321743782540412751e-14])

    s = 7
    A = A * t
    n = 0 * n_0

    for j in range(7):
        n = n + np.linalg.solve(
            A - theta[j] * np.identity(np.shape(A)[0]), alpha[j] * n_0
        )

    n = 2 * n.real
    n = n + alpha_0 * n_0

    return n


def _rat_apprx_16(A, t, n_0):
    """CRAM of order 16

    Parameters
    ---------
    A : numpy array
        Burnup matrix
    t : float
        Time step
    n_0: numpy array
        Inital composition vector
    """
    theta = np.array(
        [
            -10.843917078696988026 + 19.277446167181652284j,
            -5.2649713434426468895 + 16.220221473167927305j,
            5.9481522689511774808 + 3.5874573620183222829j,
            3.5091036084149180974 + 8.4361989858843750826j,
            6.4161776990994341923 + 1.1941223933701386874j,
            1.4193758971856659786 + 10.925363484496722585j,
            4.9931747377179963991 + 5.9968817136039422260j,
            -1.4139284624888862114 + 13.497725698892745389j,
        ]
    )

    alpha = np.array(
        [
            -0.0000005090152186522491565 - 0.00002422001765285228797j,
            0.00021151742182466030907 + 0.0043892969647380673918j,
            113.39775178483930527 + 101.9472170421585645j,
            15.059585270023467528 - 5.7514052776421819979j,
            -64.500878025539646595 - 224.59440762652096056j,
            -1.4793007113557999718 + 1.7686588323782937906j,
            -62.518392463207918892 - 11.19039109428322848j,
            0.041023136835410021273 - 0.15743466173455468191j,
        ]
    )

    alpha_0 = np.array([2.1248537104952237488e-16])

    s = 8
    A = A * t
    n = 0 * n_0

    for j in range(8):
        n = n + np.linalg.solve(
            A - theta[j] * np.identity(np.shape(A)[0]), alpha[j] * n_0
        )

    n = 2 * n.real
    n = n + alpha_0 * n_0
    return n


def cram(N, t, n_0, order):
    """This function returns matrix exponential solution n using CRAM14 or CRAM16

    Parameters
    ----------
    N : list or array
        Array of nuclides under consideration
    t : float
        Time step
    n_0 : list or array
        Nuclide concentration vector
    order : int
        Order of method. Only 14 and 16 are supported.
    """

    n_0 = np.array(n_0)
    A = _build_matrix(N)

    if order == 14:
        return _rat_apprx_14(A, t, n_0)

    elif order == 16:
        return _rat_apprx_16(A, t, n_0)

    else:
        msg = "Rational approximation of degree {0} is not supported.".format(order)
        raise ValueError(msg)


def _output_flux_block(ve, tag_flux, reverse):
    """
    This function is used to get neutron flux for fluxin

    Parameters
    ----------
    ve : entity, a mesh sub-voxel
    tag_flux : array, neutron flux of the sub-voxel
    reverse : bool, whether to reverse the flux
    """

    flux_data = np.atleast_1d(tag_flux[ve])
    if reverse:
        flux_data = np.flip(flux_data)
    outs = "".join(["{:.6E} ".format(value) for value in flux_data])
    outs = "\n".join(re.findall(".{1,78}", outs))
    outs += "\n\n"
    return outs


def _get_subvoxel_array(mesh, cell_mats):
    """
    This function returns an array of subvoxels.
    Parameters
    ----------
    mesh : PyNE Mesh object
        The Mesh object for which the geometry is discretized.

    return : subvoxel_array: structured array
        A sorted, one dimensional array, each entry containing the following
        fields:

            :svid: int
                The index of non-void subvoxel id
            :idx: int
                The idx of the voxel
            :scid: int
                The cell index of the cell in that voxel

    """
    cell_number_tag = mesh.cell_number
    subvoxel_array = np.zeros(
        0, dtype=[("svid", np.int64), ("idx", np.int64), ("scid", np.int64)]
    )
    temp_subvoxel = np.zeros(
        1, dtype=[("svid", np.int64), ("idx", np.int64), ("scid", np.int64)]
    )
    # calculate the total number of non-void sub-voxel
    non_void_sv_num = 0
    for i, _, ve in mesh:
        for c, cell in enumerate(np.atleast_1d(cell_number_tag[ve])):
            if cell > 0 and len(cell_mats[cell].comp):  # non-void cell
                temp_subvoxel[0] = (non_void_sv_num, i, c)
                subvoxel_array = np.append(subvoxel_array, temp_subvoxel)
                non_void_sv_num += 1

    return subvoxel_array


def _make_response_dtype(response_name, data_length=1):

    return np.dtype(
        [
            ("idx", np.int64),
            ("nuc", "S6"),
            ("time", "S20"),
            (response_name, np.float64, data_length),
        ]
    )


def _convert_unit_to_s(dt):
    """
    This function return a float number represent a time in unit of s.
    Parameters
    ----------
    dt : string.
        Decay time. Contain a num and an unit.

    Returns
    -------
    a float number
    """
    dt = str_to_unicode(dt)
    # get num and unit
    if dt == "shutdown":
        num, unit = "0.0", "s"
    else:
        num, unit = dt.split()
    return to_sec(float(num), unit)


def _find_dt(idt, decay_times):
    """
    This function returns a string representing a time in decay times.

    Parameters
    ----------
    idt : string
        Represents a time, input decay time
    decay_times : list of strings
        Decay times.

    Returns
    -------
    string from decay_times list that mathches idt
    """
    # Check the existence of idt in decay_times list.
    if idt in decay_times:
        return idt
    # Direct matching cannot be found. Convert units to [s] and compare.
    else:
        # convert idt to [s]
        idt_s = _convert_unit_to_s(idt)
        # Loop over decay times in decay_times list and compare to idt_s.
        for dt in decay_times:
            # Skip "shutdown" string in list.
            if str_to_unicode(dt) == "shutdown":
                continue
            # Convert to [s].
            dt_s = _convert_unit_to_s(dt)
            if idt_s == dt_s:
                # idt_s matches dt_s. return original string, dt.
                return dt
            elif dt_s != 0.0 and (abs(idt_s - dt_s) / dt_s) < 1e-6:
                return dt
        # if idt doesn't match any string in decay_times list, raise an error.
        raise ValueError("Decay time {0} not found in decay_times".format(idt))


def responses_output_zone(responses=None, wdr_file=None, alara_params=None):
    """
    This function returns a string that is an output block for an alara input
    file, configured for zone resolution.

    Parameters
    ----------
    responses : list of string
        A keyword represent the alara output zone. The following responses are
        supported:
            - decay_heat
            - specific_activity
            - alpha_heat
            - beta_heat
            - gamma_heat
            - wdr
            - photon_source
    wdr_file : string
        File name of the standard used to calculate wdr.
    alara_params: string
        Alara parameters

    Returns
    -------
    String represent the output zone block corresponding to the response.
    """

    # set default value for functions do not need response
    if responses == None:
        return ""
    # input check
    for response in responses:
        if response not in response_strings.keys():
            raise ValueError("response {0} not supported.".format(response))

    output_strings = {
        "decay_heat": "      total_heat",
        "specific_activity": "      specific_activity",
        "alpha_heat": "      alpha_heat",
        "beta_heat": "      beta_heat",
        "gamma_heat": "      gamma_heat",
    }
    if "wdr" in responses:
        output_strings["wdr"] = "".join(["       wdr ", wdr_file])
    if "photon_source" in responses:
        alara_lib = get_alara_lib(alara_params)
        output_strings["photon_source"] = "".join(
            ["      photon_source ", alara_lib, " phtn_src 1 2e7"]
        )
    output_zone = ["output zone"]
    for response in responses:
        output_zone.append(output_strings[response])
    output_zone.append("end")
    return "\n".join(output_zone)


def _is_data(line):
    """
    This function is used to check whether a line of alara output file contains
    wanted data. The line contains data is conposed of:
        - nuc : nuc name (total or total)
        - data for each decay time (including 'shutdown': floats

    Parameters
    ----------
    line : string
        A line from ALARA output.txt

    Returns
    -------
    True : if this line contains results data
    False : if this line doesn't contain results data
    """
    # check the list from the second value, if they are float, then return True
    tokens = line.strip().split()
    if len(tokens) < 2:
        return False
    # first block should be a valid nucname or 'total'
    if not (nucname.isnuclide(tokens[0]) or tokens[0] == "TOTAL".lower()):
        return False
    try:
        np.array(tokens[1:]).astype(float)
        return True
    except:
        return False


def read_decay_times(line):
    """
    This function reads a line contian decay times information from alara
    output file and return the decay times list.

    Parameters
    ----------
    line : string
        A line from ALARA output.txt

    Returns
    -------
    decay_times : array of string
        Array of decay times.
    """
    tokens = line.strip().split()
    decay_times = ["shutdown"]
    for i in range(2, len(tokens), 2):
        decay_times.append("".join([tokens[i], " ", tokens[i + 1]]))
    return decay_times


def _get_zone_idx(line):
    """
    This function is used to get the zone idx from a line of ALARA output.txt.

    Parameters
    ----------
    line : string
        A line from ALARA output.txt

    Returns
    -------
    int, zone index
    """

    last_word = line.strip().split()[-1]
    return int(last_word.split("_")[-1])


def get_alara_lib(alara_params):
    """
    This function is used to get the alara_library from alara_params.

    Parameters
    ----------
    alara_params: string
        ALARA parameters.

    Returns
    -------
    alara_lib: string
        Path to ALARA library.
    """
    lines = alara_params.split("\n")
    for line in lines:
        if "data_library" in line:
            alara_lib = line.strip().split()[-1]
            return alara_lib
    raise ValueError("alara_lib not found!")


def _read_h5_dt(filename):
    """
    This function reads decay times in h5 photon source and responses file.

    Parameters
    ----------
    filename : string
        Filename of the photon source file or response file.

    Returns
    -------
    dt : list
        List of the decay times, in unicode.
    """
    # creat a list of decay times (strings) in the source file
    dt = []
    with tb.open_file(filename) as h5f:
        for row in h5f.root.data:
            if row[2].decode() not in dt:
                dt.append(row[2].decode())
            else:
                break
    return dt
