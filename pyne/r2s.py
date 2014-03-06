from os.path import isfile

from pyne.mesh import Mesh
from pyne.mcnp import Meshtal
from pyne.alara import mesh_to_fluxin, mesh_to_geom
from pyne.dagmc import load, discretize_geom

def irradiation_setup(meshtal, cell_mats, alara_params, geom=None,
                      num_rays=10, grid=False, flux_tag="n_flux",
                      fluxin="alara_fluxin", alara_inp="alara_geom",
                      alara_matlib="alara_matlib",  output_mesh="r2s_step1.h5m")
    """def irradiation_setup(meshtal, cell_mats, alara_params, geom=None,
                          num_rays =10, grid=False, flux_tag="n_flux",
                          fluxin="alara_fluxin", alara_inp="alara_inp",
                          alara_matlib="alara_matlib", output_mesh="r2s_step1.h5m")

    This function is used to setup the irradiation inputs after the first
    transport step.

    Parameters
    ----------
    meshtal : PyNE Meshtal object or str
        The meshtal file contain the neutron flux data.
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
    alara_inp : str, optional
        The name of the ALARA input file to be created.
    alara_matlib : str, optional
        The name of the alara_matlib file to be created.
    output_mesh : str, optional
        A mesh containing all the fluxes and materials used for irradiation
        setup.
    """
    # Aquire fluxes
    is not isinstance(Meshtal, meshtal) and isfile(meshtal):
        meshtal = Meshtal(meshtal, flux_tag)

    if geom is not None and isfile(geom):
        load(geom)

    if isfile(alara_params):
        with open(alara_params, 'r') as f:
            alara_params = f.read(alara_params)

    vol_fracs = discretize_geom(meshtal.mesh, num_rays, grid=grid)
    mesh.cell_fracs_to_mats(vol_fracs, mats)
    mesh_to_fluxin(meshtal.mesh, flux_tag, fluxin, reverse)
    mesh_to_geom(meshtal.mesh, alara_geom, alara_matlib)
    
    mesh.mesh.save(output_mesh)
