from os.path import isfile

from pyne.mesh import Mesh
from pyne.mcnp import Meshtal
from pyne.alara import mesh_to_fluxin, mesh_to_geom
from pyne.dagmc import load, discretize_geom

def irradiation_setup(meshtal, tally_num, cell_mats, alara_params, geom=None,
                      num_rays=10, grid=False, flux_tag="n_flux",
                      fluxin="alara_fluxin", reverse=False, 
                      alara_inp="alara_geom", alara_matlib="alara_matlib",  
                      output_mesh="r2s_step1.h5m"):
    """def irradiation_setup(meshtal, tally_num, cell_mats, alara_params, geom=None,
                          num_rays=10, grid=False, flux_tag="n_flux",
                          fluxin="alara_fluxin", reverse=False, 
                          alara_inp="alara_geom", alara_matlib="alara_matlib",  
                          output_mesh="r2s_step1.h5m")

    This function is used to setup the irradiation inputs after the first
    transport step.

    Parameters
    ----------
    meshtal : PyNE Meshtal object or str
        The meshtal file contain the neutron flux data.
    tally_num : int
        The MCNP FMESH4 tally number of the neutron flux tally within the
        meshtal file.
    cell_mats : dict
        Maps geometry cell numbers to PyNE Material objects.
    alara_params : str
        The ALARA input blocks specifying everything except the geometry
        and materials. This can either be passed a string or as a file name.
    geom : str, optional
        The file name of a DAGMC-loadable faceted geometry. This is only
        necessary if the geometry is not already loaded into memory.
    num_rays : int, optional
        The number of rays to fire down a mesh row for geometry discretization.
        This number must be a perfect square if grid=True.
    grid : bool, optional
        The if False, geometry discretization will be done with randomly fired
        rays. If true, a grid of sqrt(num_rays) x sqrt(num_rays) rays is used
        for each mesh row.
    flux_tag : str, optional
        The iMesh tag for the neutron flux.
    fluxin : str, optional
        The name of the ALARA fluxin file to be created.
    reverse : bool, optional
        If True the fluxes in the fluxin file will be printed in the reverse
        order of how they appear within the flux vector tag. Since MCNP and
        the meshtal class order fluxes from low energy to high energy, this
        option should only be true if the transmutation data being used is
        ordered from high energy to low energy.
    alara_inp : str, optional
        The name of the ALARA input file to be created.
    alara_matlib : str, optional
        The name of the alara_matlib file to be created.
    output_mesh : str, optional
        A mesh containing all the fluxes and materials used for irradiation
        setup.
    """
    # Acquire fluxes
    if not isinstance(Meshtal, meshtal) and isfile(meshtal):
        meshtal = Meshtal(meshtal, flux_tag)

    if geom is not None and isfile(geom):
        load(geom)
    m = meshtal.tally[tally_num]
    vol_fracs = discretize_geom(m, num_rays, grid=grid)
    mesh.cell_fracs_to_mats(vol_fracs, mats)


    mesh_to_fluxin(m, flux_tag, fluxin, reverse)
    mesh_to_geom(m, inp, alara_matlib)

    if isfile(alara_params):
        with open(alara_params, 'r') as f:
            alara_params = f.read(alara_params)

    with open(alara_inp, 'a') as f:
        f.write("\n" + alara_params)
    
    m.mesh.save(output_mesh)

def photon_sampling_setup(phtn_src, mesh, tags):
    """This function reads in an hdf5 file produced by photon_source_to_hdf5

    Parameters
    ----------
    mesh : PyNE Mesh
       The object containing the imesh instance to be tagged.
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
    h5_file = "{0}.h5".format(phtn_src)
    photon_source_hdf5_to_mesh(mesh, h5_file, tags)

