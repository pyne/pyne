"""This module contains functions relevant to the ALARA activation code.
"""
import numpy as np

try:
    from itaps import iMesh, iBase, iMeshExtensions
except ImportError:
    warnings.warn("the PyTAPS optional dependency could not be imported. "
         "Some aspects of the alara module may be incomplete.", ImportWarning)

from mesh import Mesh, MeshError

def flux_mesh_to_fluxin(flux_mesh, flux_tag, fluxin=fluxin.out, reverse=False):
    """This function creates an ALARA fluxin file from fluxes tagged on a PyNE
    Mesh object.

    Parameters:
    ----------
    flux_mesh : PyNE Mesh object
        Contains the mesh with fluxes tagged on each volume element.
    flux_tag : string
        The name of the tag of the flux mesh. Flux values for different energy
        groups are assumed to be represented as vector tags.
    fluxin : string
        The name of the ALARA fluxin file to be output.
    reverse: bool
        If true, fluxes will be printed in the reverse order as they appear in
        the flux vector tagged on the mesh.
    """
    tag_flux = flux_mesh.mesh.getTagHandle(flux_tag)

    if flux_mesh.structured:
        ves = flux_mesh.structured_iterate_hex("zyx")
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
