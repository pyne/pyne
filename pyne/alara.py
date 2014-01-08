"""This module contains functions relevant to the ALARA activation code.
"""
try:
    from itaps import iMesh, iBase, iMeshExtensions
except ImportError:
    warnings.warn("the PyTAPS optional dependency could not be imported. "
         "Some aspects of the ALARA module may be incomplete.", ImportWarning)

from pyne.Mesh import Mesh, MeshError

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
    
    flux_tag = flux_mesh.getTagHandle(flux_tag)

    if flux_mesh.structured:
        ves = flux_mesh.structured_iterateHex("zyx")
    else:
        ves = flux.mesh.mesh.iterate(iBase.Type.region, iMesh.Toplogy.all)

    #Establish for loop bounds based on if forward or backward printing
    #is requested
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
        #Print flux data to file
        count = 0
        flux_data = tag_flux[ve]
        for i in range(start, stop, direction):
            output += "{0} ".format(flux_data[i])
            #fluxin formatting: create a new line after every 8th entry
            count += 1
            if count % 8 == 0:
                output += "\n"
    
        output += "\n\n"

    with open(fluxin.out) as f:
        f.write(output)
