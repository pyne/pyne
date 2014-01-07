#!/usr/bin/env python
"""
This module contains functions for mesh-based Monte Carlo variance reduction.
"""
from itaps import iMesh
from itaps import iBase
from itaps import iMeshExtensions
import numpy as np
from itertools import izip

from pyne.mesh import Mesh
from pyne.mesh import MeshError

def cadis(adj_flux_mesh, adj_flux_tag, q_mesh, q_tag, 
          ww_mesh, ww_tag, q_bias_mesh, q_bias_tag, beta=5):
    """This function reads PyNE meshes tagged with adjoint fluxes and unbiased
    source densities and outputs meshes of weight window lower bounds and
    biased source densities as computed by the Consistant Adjoint-Driven
    Importance Sampling (CADIS) method. Note values can be stored be can be
    stored on the same Mesh object, all difference Mesh objects, or any 
    combination in between. Meshes can be structured or unstructured. Note that
    this function is suitable for Forward Weighted (FW) CADIS as well, the only
    difference being the adjoint source used for the estimation of the adjoint 
    flux.

    Parameters
    -----------
    adj_flux_mesh : PyNE Mesh object
        The mesh containing the adjoint fluxes.
    adj_flux_tag : string
        The name of the adjoing flux tag on adj_mesh.
    q_mesh : PyNE Mesh object
        The mesh containing the unbiased source density.
    q_tag : string
        The source daensity tag on q_mesh.
    ww_mesh : PyNE Mesh object
        Name of the mesh to store the output weight window mesh.
    ww_tag : string
        Name of the tag to store output weight window values on ww_mesh.
    q_bias_mesh : PyNE Mesh object
        Name of the mesh to store the output bias source mesh.
    q_bias_tag : PyNE Mesh object 
        Name of the tag to store output weight window values on ww_mesh.
    beta : float
        The ratio of the weight window upper bound to the weight window lower
        bound. The default value is 5: the value used in MCNP. 
    """
    
    #find number of energy groups

    e_groups = adj_flux_mesh.mesh.getTagHandle(adj_flux_tag)[list(
                   adj_flux_mesh.mesh.iterate(iBase.Type.region, 
                                              iMesh.Topology.all))[0]]
    #if there is one e_group, convert e_groups to list
    try:
        len(e_groups)
    except TypeError:
        e_groups = [e_groups]

    num_e_groups = len(e_groups)

    #verify source (q) mesh has the same number of energy groups
    q_e_groups = q_mesh.mesh.getTagHandle(q_tag)[list(
                      q_mesh.mesh.iterate(iBase.Type.region, 
                                          iMesh.Topology.all))[0]]
    try:
        len(q_e_groups)
    except TypeError:
        q_e_groups = [q_e_groups]

    num_q_e_groups = len(q_e_groups)

    if num_q_e_groups != num_e_groups:
        raise TypeError("{0} on {1} and {2} on {3} "
              "must be of the same dimension".format(adj_flux_mesh, adj_flux_tag, 
                                                     q_mesh, q_tag))

    #calculate total response (R)
    R = [0]*num_e_groups
    #create volume element (ve) iterators
    adj_ves = adj_flux_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)
    q_ves = q_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)

    for adj_ve, q_ve in zip(list(adj_ves), list(q_ves)):
        adj_flux = adj_flux_mesh.mesh.getTagHandle(adj_flux_tag)[adj_ve]
        q = q_mesh.mesh.getTagHandle(q_tag)[q_ve]

        try:
            len(adj_flux)
        except TypeError:
            adj_flux = [adj_flux]
        try:
            len(q)
        except TypeError:
            q = [q]

        for i in range(0, num_e_groups):
            R[i] += adj_flux[i]*q[i]

    #generate weight windows and biased source densities using total response
    tag_ww = ww_mesh.mesh.createTag(ww_tag, num_e_groups, float)
    ww_ves = ww_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)
    tag_q_bias = q_bias_mesh.mesh.createTag(q_bias_tag, num_e_groups, float) 
    q_bias_ves = q_bias_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)
    q_ves.reset()
    adj_ves.reset()

    for adj_ve, q_ve, ww_ve, q_bias_ve in izip(list(adj_ves), list(q_ves), 
                                               list(ww_ves), list(q_bias_ves)):
        adj_flux = adj_flux_mesh.mesh.getTagHandle(adj_flux_tag)[adj_ve]
        q = q_mesh.mesh.getTagHandle(q_tag)[q_ve]

        try:
            len(adj_flux)
        except TypeError:
            adj_flux = [adj_flux]
        try:
            len(q)
        except TypeError:
            q = [q]

        tag_ww[ww_ve] = [R[i]/(adj_flux[i]*q[i]*(beta + 1)/2) 
                         for i in range(0, num_e_groups)]

        tag_q_bias[q_bias_ve] = [adj_flux[i]*q[i]/R[i] 
                                 for i in range(0, num_e_groups)]

#    tag_ww = ww_mesh.mesh.createTag(ww_tag, num_e_groups, float)
#    ww_ves = ww_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)
#    q_ves.reset()
#    adj_ves.reset()
#
#    for ww_ve, adj_ve, q_ve in izip(list(ww_ves), list(adj_ves), list(q_ves)):
#        adj_flux = adj_flux_mesh.mesh.getTagHandle(adj_flux_tag)[adj_ve]
#        q = q_mesh.mesh.getTagHandle(q_tag)[q_ve]
#
#        try:
#            len(adj_flux)
#        except TypeError:
#            adj_flux = [adj_flux]
#        try:
#            len(q)
#        except TypeError:
#            q = [q]
#
#        tag_ww[ww_ve] = [R[i]/(adj_flux[i]*q[i]*(beta + 1)/2) 
#                         for i in range(0, num_e_groups)]
#
#    tag_q_bias = q_bias_mesh.mesh.createTag(q_bias_tag, num_e_groups, float) 
#    q_bias_ves = q_bias_mesh.mesh.iterate(iBase.Type.region, iMesh.Topology.all)
#    q_ves.reset()
#    adj_ves.reset()
#
#    for q_bias_ve, adj_ve, q_ve in izip(list(q_bias_ves), list(adj_ves), list(q_ves)):
#        adj_flux = adj_flux_mesh.mesh.getTagHandle(adj_flux_tag)[adj_ve]
#        q = q_mesh.mesh.getTagHandle(q_tag)[q_ve]
# 
#        try:
#            len(adj_flux)
#        except TypeError:
#            adj_flux = [adj_flux]
#        try:
#            len(q)
#        except TypeError:
#            q = [q]
#
#        tag_q_bias[q_bias_ve] = [adj_flux[i]*q[i]/R[i] 
#                                 for i in range(0, num_e_groups)]
