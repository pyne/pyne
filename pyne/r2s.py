from os.path import isfile
from warnings import warn
from pyne.utils import QAWarning
import numpy as np

from pyne.mesh import Mesh
from pyne.mcnp import Meshtal
from pyne.alara import mesh_to_fluxin, record_to_geom, photon_source_to_hdf5, \
    photon_source_hdf5_to_mesh

warn(__name__ + " is not yet QA compliant.", QAWarning)


def resolve_mesh(mesh_reference, tally_num=None, flux_tag="n_flux",
                 output_material=False):
    """This function creates a method that will consume many mesh-like objects
       (e.g. mesh, an h5m file, a meshtal file, etc) and returns a robust PyNE
       mesh object accordingly.

    Parameters
    ----------
    mesh_reference : Mesh object, unstructured mesh file, Meshtal, meshtal
         file, or PyNE Meshtal object.
        The source of the neutron flux information. This can be a PyNE Meshtal
        object, a pyne Mesh object, or the filename an MCNP meshtal file, or
        the filename of an unstructured mesh tagged with fluxes.
    tally_num : int
        The MCNP FMESH4 tally number of the neutron flux tally within the
        meshtal file.
    flux_tag : str, optional
        The tag name for the neutron flux.
    output_material : bool, optional
        If true, output mesh will have materials as determined by
        dagmc.discretize_geom().

    Returns
    -------
    m : PyNE mesh object
        The PyNE mesh object of the flux data.
    """

    # mesh_reference is Mesh object
    if isinstance(mesh_reference, Mesh):
        m = mesh_reference
    #  mesh_reference is unstructured mesh file
    elif isinstance(mesh_reference, str) and isfile(mesh_reference) \
            and mesh_reference.endswith(".h5m"):
        m = Mesh(structured=False, mesh=mesh_reference)
    #  mesh_reference is Meshtal or meshtal file
    elif tally_num is not None:
        #  mesh_reference is meshtal file
        if isinstance(mesh_reference, str) and isfile(mesh_reference):
            mesh_reference = Meshtal(mesh_reference,
                                     {tally_num: (flux_tag, flux_tag + "_err",
                                                  flux_tag + "_total",
                                                  flux_tag + "_err_total")},
                                     meshes_have_mats=output_material)
            m = mesh_reference.tally[tally_num]
        #  mesh_reference is Meshtal object
        elif isinstance(mesh_reference, Meshtal):
            m = mesh_reference.tally[tally_num]
        else:
            raise ValueError("meshtal argument not a Mesh object, Meshtal"
                             " object, MCNP meshtal file or meshtal.h5m file.")
    # mesh_references is a Meshtal file but no tally_num provided
    else:
        raise ValueError(
            "Need to provide a tally number when reading a Meshtal file")

    return m


def irradiation_setup(flux_mesh, cell_mats, cell_fracs, alara_params,
                      tally_num=4, num_rays=10, grid=False, flux_tag="n_flux",
                      fluxin="alara_fluxin", reverse=False,
                      alara_inp="alara_inp", alara_matlib="alara_matlib",
                      output_mesh="r2s_step1.h5m", output_material=False,
                      decay_times=None, sub_voxel=False):
    """This function is used to setup the irradiation inputs after the first
    R2S transport step.

    Parameters
    ----------
    flux_mesh : PyNE Meshtal object, Mesh object, or str
        The source of the neutron flux information. This can be a PyNE Meshtal
        object, a pyne Mesh object, or the filename an MCNP meshtal file, or
        the filename of an unstructured mesh tagged with fluxes.
    cell_mats : dict
        Maps geometry cell numbers to PyNE Material objects.
    cell_fracs : record array
        The output of dagmc.discretize_geom().
    alara_params : str
        The ALARA input blocks specifying everything except the geometry
        and materials. This can either be passed as string or as a file name.
    tally_num : int
        The MCNP FMESH4 tally number of the neutron flux tally within the
        meshtal file.
    num_rays : int, optional
        The number of rays to fire down a mesh row for geometry discretization.
        This number must be a perfect square if grid=True.
    grid : bool, optional
        The if False, geometry discretization will be done with randomly fired
        rays. If true, a grid of sqrt(num_rays) x sqrt(num_rays) rays is used
        for each mesh row.
    flux_tag : str, optional
        The mesh tag for the neutron flux.
    fluxin : str, optional
        The name of the ALARA fluxin file to be created.
    reverse : bool, optional
        If True the fluxes in the fluxin file will be printed in the reverse
        order of how they appear within the flux vector tag. Since MCNP and
        the Meshtal class order fluxes from low energy to high energy, this
        option should only be true if the transmutation data being used is
        ordered from high energy to low energy.
    alara_inp : str, optional
        The name of the ALARA input file to be created.
    alara_matlib : str, optional
        The name of the alara_matlib file to be created.
    output_mesh : str, optional
        A mesh containing all the fluxes and materials used for irradiation
        setup.
    output_material : bool, optional
        If true, output mesh will have materials as determined by
        dagmc.discretize_geom()
    decay_times: list
        List of the decay times. If no decay times given, use '1 s'.
    sub_voxel : bool, optional
        If true, sub-voxel r2s work flow  will be used.
    """

    m = resolve_mesh(flux_mesh, tally_num, flux_tag, output_material)

    if output_material:
        m.cell_fracs_to_mats(cell_fracs, cell_mats)

    mesh_to_fluxin(m, flux_tag, fluxin, reverse,
                   sub_voxel, cell_fracs, cell_mats)
    record_to_geom(m, cell_fracs, cell_mats, alara_inp, alara_matlib,
                   sub_voxel=sub_voxel)

    # write decay times into alara_inp
    if decay_times is None:
        decay_times = ['1 s']
    decay_str = 'cooling\n'
    for dc in decay_times:
        decay_str = ''.join([decay_str, '    ', dc, '\n'])
    decay_str = ''.join([decay_str, 'end\n'])
    with open(alara_inp, 'a') as f:
        f.write(decay_str)

    if isfile(alara_params):
        with open(alara_params, 'r') as f:
            alara_params = f.read()

    with open(alara_inp, 'a') as f:
        f.write("\n" + alara_params)

    m.write_hdf5(output_mesh)


def photon_sampling_setup(mesh, phtn_src, tags):
    """This function reads in an ALARA photon source file and creates and tags
    photon source densities onto a Mesh object for the second R2S transport
    step.

    Parameters
    ----------
    mesh : PyNE Mesh
       The object containing the mesh instance to be tagged.
    phtn_src : str
        The path of the ALARA phtn_file.
    tags: dict
        A dictionary were the keys are tuples with two values. The first is a
        string denoting an nuclide in any form that is understood by
        pyne.nucname (e.g. '1001', 'U-235', '242Am') or 'TOTAL' for all
        nuclides. The second is a string denoting the decay time as it appears
        in the phtn_src file (e.g. 'shutdown', '1 h' '3 d'). The values of the
        dictionary are the requested tag names for the combination of nuclide
        and decay time. These tag names should be the tag names that are read
        by the sampling subroutine. For example:

        tags = {('U-235', 'shutdown'): 'tag1', ('TOTAL', '1 h'): 'tag2'}
    """
    photon_source_to_hdf5(phtn_src)
    h5_file = phtn_src + ".h5"
    photon_source_hdf5_to_mesh(mesh, h5_file, tags)


def total_photon_source_intensity(m, tag_name, sub_voxel=False):
    """This function reads mesh tagged with photon source densities and returns
    the total photon emission desinty.

    Parameters
    ----------
    m : PyNE Mesh
       The mesh-based photon emission density distribution in p/cm3/s.
    tag_name : str
       The name of the tag on the mesh with the photon emission density
       information.
    sub_voxel: bool, optional
        If true, sub-voxel r2s work flow will be used.

    Returns
    -------
    intensity : float
        The total photon emission density across the entire mesh (p/s).
    """

    sd_tag = m.get_tag(tag_name)
    intensity = 0.0
    if sub_voxel:
        cell_fracs = m.cell_fracs[:]
    else:
        # create a cell_fracs
        cell_fracs = np.ones(shape=(len(m), 1), dtype=float)
    max_num_cells = len(cell_fracs[0])
    num_e_groups = len(sd_tag[list(m.iter_ve())[0]]) / max_num_cells
    for idx, _, ve in m:
        ve_data = sd_tag[ve]
        for svid in range(max_num_cells):
            vol = m.elem_volume(ve) * cell_fracs[idx][svid]
            sv_data = ve_data[num_e_groups*svid:num_e_groups*(svid+1)]
            intensity += vol * np.sum(sv_data)
    return intensity
