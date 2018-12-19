"""Purpose:

This module provides a modified ray traverse method, which is modifed based on
http://www.cse.yorku.ca/~amana/research/grid.pdf, to calculate which voxels
intersect with a given ray (line segment).
Several rays are used to calculate the voxels that intersect with the facet.
A series of rays are constructed based on both the facet and the mesh grid to
make sure all the wanted voxels will be found. Combining all the found voxels
for each ray, we get the voxels intersecting with the facet.

These voxels will be used to construct the bounding volume for source sampling.
But currently these functions are separate from MOAB. Interface of MOAB with
PyNE is still in need.

.. moduleauthor:: X. Zhang
"""
from __future__ import division
import math
import numpy as np


##################################
### Facet voxel traverse tools ###
##################################


def _is_point_in_mesh(bounds, point):
    """
    Check whether a point is in the mesh.

    Parameters:
    -----------
    bounds : numpy array
        Boundaries. A 2D array of the boundaries along x, y, and z directions.
        Expected format: [[x_min, ..., x_max],
                          [y_min, ..., y_max],
                          [z_min, ..., z_max]]
    point : numpy array
        An array with size of 3, represents the point to be checked.

    Return:
    -------
    True : when point in the mesh.
    False : when point outside the mesh.
    """
    if not (bounds[0][0] <= point[0] <= bounds[0][-1]) or \
       not (bounds[1][0] <= point[1] <= bounds[1][-1]) or \
       not (bounds[2][0] <= point[2] <= bounds[2][-1]):
        return False
    else:
        return True

def _find_voxel_idx_1d(bounds_1d, cor, vec_1d):
    """
    Find the voxel id (in a given dimension).

    Parameters:
    -----------
    bounds_1d : numpy array
        Boundary of the mesh in a given direction.
    cor : float
        Coordinate of the point in the given direction.
    vec_1d : float
        Vector attribute in the given direction. Used for determining the
        search direction.

    Return:
    -------
    idx : int
        Mesh index of the given direction. This value is -1 when the point not
        in the mesh.
    ----------
    """
    idx = -1
#    bounds_1d = np.array(bounds_1d)
    if vec_1d > 0:
        idx = np.searchsorted(bounds_1d, cor, side='right')
    else:
        idx = np.searchsorted(bounds_1d, cor, side='left')
    if idx == 0 or idx == len(bounds_1d):
    # idx could be -1, that means point outside the mesh
        return -1
    else:
        return idx - 1

def _calc_vec_dir(start, end):
    """
    Calculate the direction from start to end (start -> end).
    
    Parameters:
    -----------
    start : numpy array
        An array of size 3. The start point.
    end :
        An array of size 3. The end point.

    Return:
    -------
    vec : numpy array
        Direction from start to end, an unit vector.
    """
    return np.divide(end-start, np.linalg.norm(end-start))

def _find_next_grid_1d(cor, vec_1d, bounds_1d):
    """
    Find the next grid coordinate of a direction.

    Parameters:
    -----------
    cor : float
        Coordinate of current point in the direction.
    vec_1d : float
        Direction value.
    bounds_1d : numpy array
        Boundary values of current direction.

    Return:
    -------
    n_grid : float
        The coordinate value of next mesh grid.
    """
    voxel_idx = _find_voxel_idx_1d(bounds_1d, cor, vec_1d)
    if voxel_idx == -1:
        return vec_1d * np.inf
    else:
        if vec_1d > 0:
            return bounds_1d[voxel_idx+1]
        else:
            return bounds_1d[voxel_idx]

def _calc_max_travel_length_in_current_voxel(point, end, vec, bounds):
    """
    Calculate the maximum travel length in current voxel before reach next
    boundary.
    
    Parameters:
    -----------
    point : numpy array
        An array of size 3. Temporary point position.
    end : numpy array
        An array of size 3. End point
    vec : numpy array
        Direction vector from start to end: [v_x, v_y, v_z].
    bounds : numpy array
        Boundaries of the mesh grid.
        Expected format: [[x_min, ..., x_max],
                          [y_min, ..., y_max],
                          [z_min, ..., z_max]]

    Return:
    -------
    dist_maxs : numpy array
        The maximum travel lenght in current voxel for 3 directions. Format: 
        [dist_max_x, dist_max_y, dist_max_z]. These value could be both
        positive and negtive. And could be 'inf' or '-inf'.
    """
    n_x_grid = _find_next_grid_1d(point[0], vec[0], bounds[0])
    n_y_grid = _find_next_grid_1d(point[1], vec[1], bounds[1])
    n_z_grid = _find_next_grid_1d(point[2], vec[2], bounds[2])
    next_grid = np.array([n_x_grid, n_y_grid, n_z_grid])
    # Divide by direction cosines to calculate the distance along vec
    dist_maxs = np.divide(next_grid - point, vec)
    return dist_maxs

def _move_to_next_voxel(idxs, point, dist_temp, dist_maxs, vec):
    """
    Move the point to next voxel boundary.

    Parameters:
    -----------
    idxs: list
        List of indexes, format: [x_idx, y_idx, z_idx]
    point : numpy array
        An array of size 3. Current point position.
    dist_temp: float
        Temporary travel length from the start.
    dist_maxs: numpy array
        The maximum travel length in current voxel for 3 directions.
    vec : numpy array
        Direction vector from start to end.

    Returns:
    --------
    idxs: list
        Updated indexes.
    updated_point : numpy array
        An array of size 3. Updated point positon.
    dist_temp : float
        Updated travel length since the beginning of the ray."""

    dist_min = min(dist_maxs[dist_maxs>0])
    update_dir = list(dist_maxs).index(dist_min)
    # vec[update_dir] will not be 0
    if vec[update_dir] > 0:
        idxs[update_dir] += 1
    elif vec[update_dir] < 0:
        idxs[update_dir] -= 1 
    updated_point = np.add(point, dist_min * vec)
    dist_temp += dist_min
    return idxs, updated_point, dist_temp

def _is_ray_on_boundary(start, vec, bounds):
    """
    Check whether the ray on a boundary surface.

    Parameters:
    -----------
    start : numpy array
        Start point coordinate.
    vec : numpy array
        Direction vector.
    bounds : numpy array
        Boundaries of the mesh grid.

    Return:
    -------
    True : bool
        Ray is on a boundary surface
    False : bool
        Ray is not on a boundary surface
    """
    if 1.0 in vec:
        for idx in np.argwhere(vec!=1.0).flatten():
            if start[idx] in bounds[idx]:
                return True
    return False

def _ray_voxel_traverse(bounds, start, end):
    """
    Ray (line segment) traversal in structured mesh.
    Return the voxels that intersect with the line segment.

    Parameters:
    -----------
    bounds : numpy array
        Boundaries of the mesh grid.
    start : numpy array
        Start point.
    end : numpy array
        End point.

    Return:
    -------
    voxles : list of tuples
        List of voxel indices. These voxels intersect with the line segment.
    """
    # Initialization
    # Find the voxel that start point in
    point = np.copy(start)
    idxs = [-1, -1, -1]
    voxels = set()
    vec = _calc_vec_dir(start, end)
    # if the ray on a boundary, then return empty set
    if _is_ray_on_boundary(start, vec, bounds):
        return voxels
    dist_end = np.linalg.norm(end-start)
    dist_temp = 0.0
    while dist_temp < dist_end:
        if _is_point_in_mesh(bounds, point):
            # the point is in the mesh
            # calculate the voxel id
            for dr in range(len(idxs)):
                idxs[dr] = _find_voxel_idx_1d(bounds[dr], point[dr], vec[dr])
            if -1 in idxs:
                # point on mesh boundary with direction toward outside the mesh
                return voxels
            # add current voxel
            voxels.add((idxs[0], idxs[1], idxs[2]))
            dist_maxs = _calc_max_travel_length_in_current_voxel(
                point, end, vec, bounds)
            idxs, point, dist_temp = _move_to_next_voxel(
                idxs, point, dist_temp, dist_maxs, vec)
        else:
            # Check whether the ray intersecs the mesh. Calculate the maximum travel
            # length of the ray to the outside boundaries of the mesh. If all the
            # lengths are negtive, then there is no intersection, return. Ohterwise,
            # choose the minimum length and move the point to that position.
            # left, right, back, front, down, up
            dist_maxs = np.zeros(6)
            for dr in range(len(idxs)):
                dist_maxs[dr*2] = np.divide(bounds[dr][0] - point[dr],
                                            vec[dr])
                dist_maxs[dr*2+1] = np.divide(bounds[dr][-1] - point[dr],
                                              vec[dr])
            if any((x > 0).all() for x in dist_maxs):
                # current point outside the mesh
                # there maybe intersects to the mesh
                dist_min = min(dist_maxs[dist_maxs > 0])
                # Move the point to the boundary, calculate the new coordinate
                point = np.add(point, dist_min * vec)
            else:
                # current point outside the mesh, not any intersection
                return voxels
    return voxels
 
def _calc_min_grid_step(bounds):
    """
    Calculate the minimum grid step of the given mesh.

    Parameters:
    -----------
    bounds : numpy array
        Boundaries of mesh grid.

    Return:
    -------
    min_step: float
        The minimum grid step.
    """
    min_step = np.inf
    for bound_list in bounds:
        min_step = min(min_step, min(np.array(bound_list[1:]) -
                                     np.array(bound_list[:-1])))
    return min_step

def _divide_tri_edge(start, end, step):
    """
    Divide the shortest triangel edge into several parts uniformly. Return the
    insert points.

    Parameters:
    -----------
    start : numpy array
        An array of size 3, represents the start point.
    end : numpy array
        An array of size 3, represents the start point.
    step : float
        The minimum grid step. Used to create the divisor of the shortest edge.

    Return:
    -------
    intert_points : list
        List of the ray end points. These points are located on the line
        segment from start
        to end.
    """
    dist = np.linalg.norm(end-start)
    vec = _calc_vec_dir(start, end)
    dist_temp = 0.0
    insert_points = []
    # +2 here to make sure it's a conservative result
    num_segment = math.floor(dist/step) + 2
    x_list = np.linspace(start[0], end[0], num_segment)
    y_list = np.linspace(start[1], end[1], num_segment)
    z_list = np.linspace(start[2], end[2], num_segment)
    for (x, y, z) in zip(x_list, y_list, z_list):
        insert_points.append(np.array([x, y, z]))
    return insert_points

def _create_rays_from_triangle_facet(A, B, C, min_step):
    """
    Create rays from triangle facet for traversing.
    Divide the shortest edge of the triangle with the min_step.

    Parameters:
    -----------
    A : numpy array
        An array of size 3, represents a vertex of the triangle
    B : numpy array
        An array of size 3, represents a vertex of the triangle
    C : numpy array
        An array of size 3, represents a vertex of the triangle
    min_step: float
        The minimum step of the mesh grid. Used to divide the shortest edge
        of the facet triangle.

    Return:
    -------
    rays : list of rays, [start, end] pairs
        List of the created rays. Used to represent the facet.
    """
    a = np.linalg.norm(B-C)
    b = np.linalg.norm(C-A)
    c = np.linalg.norm(A-B)
    min_edge = min([a, b, c])
    if a == min_edge:
        target = A
        source = (B, C)
    elif b == min_edge:
        target = B
        source = (A, C)
    else:
        target = C
        source = (A, B)

    insert_points = _divide_tri_edge(source[0], source[1], min_step)
    rays = [(target, point) for point in insert_points ]
    return rays


def _facet_voxel_traverse(A, B, C, bounds):
    """
    Calculate all the voxles that intersect with the
    facet. Facet is represented with 3 triangle vertices.

    Parameters:
    -----------
    A : numpy array
        An array of size 3, represents a vertex of the triangle.
        Vertex of the triangle facet
    B : numpy array
        An array of size 3, represents a vertex of the triangle.
    C : numpy array
        An array of size 3, represents a vertex of the triangle.
    bounds : numpy array
        Boundaries of the mesh grid.

    Return:
    -------
    voxel_list: set
        Set of voxels that intersect with the facet.
    """
    min_step = _calc_min_grid_step(bounds)
    rays = _create_rays_from_triangle_facet(A, B, C, min_step)
    voxels = set()
    for ray in rays:
        temp_voxels = _ray_voxel_traverse(bounds, ray[0], ray[1])
        # merge two lists
        voxels = voxels.union(temp_voxels)
    return voxels

def _is_tri_intersects_box(triangle, box_center, box_extents):
    """
    Check whether the triangle intersects with the axis aligned bounding box.

    Based on the "Fast 3D Triangle-Box Overlap Test" by Tomas Akenine-Moller.
    Related paper could see:
    https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tribox.pdf.


    Background about separate axis theorem (SAT) could see:
    http://www.jkh.me/files/tutorials/Separating%20Axis%20Theorem%20for%20Oriented%20Bounding%20Boxes.pdf

    This function is basically converted from the two codes below:
    http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox3.txt
    https://gist.github.com/zvonicek/fe73ba9903f49d57314cf7e8e0f05dcf

    Parameters:
    -----------
    triangle : numpy array
        A numpy array with shape of (3, 3), represents the vertices of the
        triangle. Format: [[v0x, v0y, v0z], [v1x, v1y, v1z], [v2x, v2y, v2z]].
    box_center : numpy array
        An array of size 3, represents center of the box.
    box_extents : numpy array
        An array of size 3, represents format: [width, length, height]

    Return:
    -------
    True : bool
        The triangle intersects with the box.
    False : bool
        The triangle does not intersect with the box.
    """

    # variables names used in the code
    # e, unit vector of axis. e[0]=(1, 0, 0), e[1]=(0, 1, 0), e[2]=(0, 0, 1).
    # v, vertices of the triangles. v[0], v[1], and v[2].
    # f, edge vector of the triangle. f[0]=v[1]-v[0],
    #                                 f[1]=v[2]-v[1],
    #                                 f[2]=v[0]-v[2].
    # a, axis for test in bullet 3. a[i][j] = e[i] x f[j], 
    #                         i.e., a[0][0] = e[0] x f[0]

    # define basis axis
    e = np.eye(3, 3)
    # Translate triangle as conceptually moving AABB to origin
    v = triangle - box_center
    # Compute edge vectors for triangle
    f = np.zeros(shape=(3, 3))
    f[0] = triangle[1] - triangle[0]
    f[1] = triangle[2] - triangle[1]
    f[2] = triangle[0] - triangle[2]

    ## region Test axes a00..a22 (category 3)
    for i in range(3):
        for j in range(0, 3):
            # calculate the mormal to the edge-axis plane
            a = np.cross(e[i], f[j])
            # project the triangle vertices onto axis a
            p = np.dot(v, a)
            # project the box onto the axis a.
            r = np.dot(box_extents, np.absolute(a))
            # check whether the box there is overlap of projection on the axis a
            if min(p) > r or max(p) < -r:
                return False
    ## endregion category 3

    ## region Test the three axes corresponding to the face normals of AABB b (category 1)
    # Exit if...
    for i in range(3):
        if max(v[:, i]) < -box_extents[i] or min(v[:, i]) > box_extents[i]:
            return False
    ## endregion category 1

    ## region Test separating axis corresponding to triangle face normal (category 2)
    plane_normal = np.cross(f[0], f[1])
    plane_distance = np.dot(plane_normal, v[0])
    # Compute the projection interval radius of b onto L(t) = b.c + t * p.n
    r = np.dot(box_extents, np.absolute(plane_normal))
    # Intersection occurs when plane distance falls within [-r,+r] interval
    if plane_distance > r:
        return False
    ## endregion category 2

    return True

