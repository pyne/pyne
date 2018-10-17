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


def _distance(start, end):
    """
    Calculate the distance between two points.

    Parameters:
    -----------
    start : numpy array 
        An array of size 3, represents the start point of  the line segment
    end : numpy array
        An array of size 3, represents the end point of  the line segment

    Return:
    -------
    dist : float
        The distance between the start point and the end point.
    """
    dist = math.sqrt((end[0] - start[0]) * (end[0] - start[0]) +
                     (end[1] - start[1]) * (end[1] - start[1]) +
                     (end[2] - start[2]) * (end[2] - start[2]))
    return dist

def _is_point_in_mesh(bounds, point):
    """
    Check whether a point is in the mesh.

    Parameters:
    -----------
    bounds : list
        Boundaries. Nested list of the boundaries along x, y, and z directions.
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
    if bounds[0][0] <= point[0] <= bounds[0][-1] and \
       bounds[1][0] <= point[1] <= bounds[1][-1] and \
       bounds[2][0] <= point[2] <= bounds[2][-1]:
        return True
    else:
        return False

def _find_voxel_idx_1d(bounds_1d, cor, vec_1d):
    """
    Find the voxel id (in a given dimension).

    Parameters:
    -----------
    bounds_1d : list
        Boundary of the mesh in a given direction.
    cor : float
        Coordinate of the point in the given direction.
    vec_1d : float
        Vector attribute in the given direction.

    Return:
    -------
    idx : int
        Mesh index of the given direction. This value is -1 when the point not
        in the mesh.
    ----------
    """
    idx = -1
    if vec_1d > 0:
        for i in range(len(bounds_1d)-1):
            if bounds_1d[i] <= cor < bounds_1d[i+1]:
                idx = i
                break
    else:
        for i in range(len(bounds_1d)-1, -1, -1):
            if bounds_1d[i] < cor <= bounds_1d[i+1]:
                idx = i
                break
    # idx could be -1, that means point outside the mesh
    return idx

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
    vec = end - start
    dist = math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])
    vec = np.divide(vec, dist)
    return vec

def _find_next_grid_1d(cor, vec_1d, bounds_1d):
    """
    Find the next grid coordinate of a direction.

    Parameters:
    -----------
    cor : float
        Coordinate of current point in the direction.
    vec_1d : float
        Direction value.
    bounds_1d : list of float
        Boundary values of current direction.

    Return:
    -------
    n_grid : float
        The coordinate value of next mesh grid.
    """
    n_grid = 0.0
    flag = False
    if vec_1d > 0:
        for bnd_pos in bounds_1d:
            if cor < bnd_pos:
                n_grid = bnd_pos
                flag = True
                break
    else:
        for i in range(len(bounds_1d)-1, -1, -1):
            if cor > bounds_1d[i]:
                n_grid = bounds_1d[i]
                flag = True
                break
    if flag:
        return n_grid
    else:
        # point goes out of the mesh
        return vec_1d * float('inf')


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
    bounds : list
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
    dist_max_x = np.divide(n_x_grid - point[0], vec[0])
    dist_max_y = np.divide(n_y_grid - point[1], vec[1])
    dist_max_z = np.divide(n_z_grid - point[2], vec[2])
    return np.array([dist_max_x, dist_max_y, dist_max_z])

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

    dist_min = min(i for i in dist_maxs if i > 0)
    update_dir = list(dist_maxs).index(dist_min)
    update_idx = [0, 0, 0]
    # vec[update_dir] will not be 0
    if vec[update_dir] > 0:
        update_idx[update_dir] = 1
    elif vec[update_dir] < 0:
        update_idx[update_dir] = -1
    idxs = [x + y for x, y in zip(idxs, update_idx)] 
    updated_point = np.add(point, dist_min * vec)
    dist_temp += dist_min
    return idxs, updated_point, dist_temp

def _move_to_boundary(point, dist_min, vec):
    """
    Move the point to the nearest boundary when initializing the problem.

    Parameters:
    -----------
    point : numpy array
        An array of size 3. Point.
    dist_min : float
        Distance to the nearest (positive) boundary. For a 3D mesh, there are 6
        outside boundaries. Therefore there are 6 value of the distance from the
        start to the boundaries. This dist_min is the smallest one among the
        positive distance. (The distance is negtive when the direction is
        negtive.)
    vec : numpy array
        Direction vector.

    Return:
    -------
    updated_point : numpy array
        Modified point. This point lies on the nearest boundary. Therefore, we
        moved point from outside the mesh to the mesh boundary (or the extension
        of the boundary, which still outside the mesh).
    """

    updated_point = np.add(point, dist_min * vec)
    return updated_point

def _is_ray_on_boundary(start, vec, bounds):
    """
    Check whether the ray on a boundary surface.

    Parameters:
    -----------
    start : numpy array
        Start point coordinate.
    vec : numpy array
        Direction vector.
    bounds : list
        Boundaries of the mesh grid.

    Return:
    -------
    True : bool
        Ray is on a boundary surface
    False : bool
        Ray is not on a boundary surface
    """
    if any((x == 1.0).all() for x in vec):
        # this is a axis aligned ray
        for i in range(len(vec)):
            if vec[i] == 1.0:
                continue
            else:
               if any((x == start[i]).all() for x in bounds[i]):
                   return True
        return False
    else:
        return False

def _ray_voxel_traverse(bounds, start, end):
    """
    Ray (line segment) traversal in structured mesh.
    Return the voxels that intersect with the line segment.

    Parameters:
    -----------
    bounds : list
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
        return set()
    dist_end = _distance(start, end)
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
            dist_maxs = [0.0] * 6
            for dr in [0, 1, 2]:
                dist_maxs[dr*2] = np.divide(bounds[dr][0] - point[dr],
                                            vec[dr])
                dist_maxs[dr*2+1] = np.divide(bounds[dr][-1] - point[dr],
                                              vec[dr])
            if any((x > 0).all() for x in dist_maxs):
                # current point outside the mesh
                # there maybe intersects to the mesh
                dist_min = min(i for i in dist_maxs if i > 0)
                # Move the point to the boundary, calculate the new coordinate
                point = _move_to_boundary(point, dist_min, vec)
            else:
                # current point outside the mesh, not any intersection
                return voxels
    return voxels
 
def _calc_min_grid_step(bounds):
    """
    Calculate the minimum grid step of the given mesh.

    Parameters:
    -----------
    bounds : list
        Boundaries of mesh grid.

    Return:
    -------
    min_step: float
        The minimum grid step.
    """
    min_step = float('inf')
    for bound_list in bounds:
        for bound_idx in range(len(bound_list)-1):
            min_step = min(min_step,
                           bound_list[bound_idx+1] - bound_list[bound_idx])
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
    dist = _distance(start, end)
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

def _create_rays_from_points(start, insert_points):
    """
    Create rays (a pair of two points) from given point list. In which, 'start'
    is the start point of all the rays, and the point in the `insert_points'
    are the end points.

    Parameters:
    -----------
    start : numpy array
        An array of size 3, represents the start point.
    insert_points : list of arrays
        List of points.

    Return:
    -------
    rays: list
        List of rays (line segments).
    """
    rays = []
    for point in insert_points:
        rays.append([start, point])
    return rays

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
    a = _distance(B, C)
    b = _distance(A, C)
    c = _distance(A, B)
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
    rays = _create_rays_from_points(target, insert_points)
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
    bounds : list
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

