"""This module provides physical cross-section models and helper functions."""

cimport numpy as np
import numpy as np

from pyne cimport nucname
from pyne import nucname


##############################
### Partial group collapse ###
##############################

def partial_energy_matrix_mono(np.ndarray[np.float64_t, ndim=1] E_g, np.ndarray[np.float64_t, ndim=1] E_n, int slope=-1):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    Here, both of the energy arrays must be monotonic. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : 1d numpy float array 
        Lower resolution energy group structure [MeV] that is of length G+1. 
        Ordered based on slope.
    E_n : 1d numpy float array 
        Higher resolution energy group structure [MeV] that is of length N+1. 
        Ordered based on slope.
    slope : int, optional
        Gives the monotonicity of E_g and E_n.  If positive, then they are 
        monotonicly increasing (lowest-to-highest).  If negative, they are
        monotonicly decreasing (highest-to-lowest).

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    # Some convienence parameters
    cdef Py_ssize_t G, N, g, lig, uig
    cdef np.ndarray[np.float64_t, ndim=2] pem

    G = E_g.shape[0] - 1
    N = E_n.shape[0] - 1

    index_E_n = np.arange(N+1)

    # Get the interior points for each gth group in n-space
    if slope < 0:
        inner_mask = np.array([(E_n <= E_g[g]) & (E_g[g+1] <= E_n) for g in range(G)])
    elif 0 < slope:
        inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])
    else:
        raise ValueError("slope must be positive or negative.")

    # Get the upper and lower nth index for every gth group
    lower_index = np.array([index_E_n[inner_mask[g]][0] for g in range(G)])
    upper_index = np.array([index_E_n[inner_mask[g]][-1] for g in range(G)])

    # Convert the mask to initialize the partial enery matrix
    # Hack off the last index of the mask to make the right size
    pem = np.array(inner_mask[:, :-1], dtype=float)

    # Check for partial contibutions at the edges
    for g in range(G):
        # Lower bound
        lig = lower_index[g]
        if lig != 0:
            pem[g,lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g,uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem



def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    The group structures must have the same monotonicity. This is useful for performing group collapses.

    Parameters
    ----------
    E_g : sequence of floats
        Lower resolution energy group structure [MeV] that is of length G+1. 
    E_n : sequence of floats
        Higher resolution energy group structure [MeV] that is of length N+1. 

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    cdef np.ndarray[np.float64_t, ndim=2] pem

    E_g = np.asarray(E_g, dtype=float)
    E_n = np.asarray(E_n, dtype=float)

    if (E_g[:-1] > E_g[1:]).all() and (E_n[:-1] > E_n[1:]).all():
        # Both energy arrays are monotonically decreasing
        assert E_g[0] <= E_n[0]
        assert E_n[-1] <= E_g[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, -1)
    elif (E_g[:-1] < E_g[1:]).all() and (E_n[:-1] < E_n[1:]).all():
        # Both energy arrays are monotonically increasing
        assert E_n[0] <= E_g[0]
        assert E_g[-1] <= E_n[-1]
        pem = partial_energy_matrix_mono(E_g, E_n, 1)
    else:
        raise ValueError("E_g and E_n are not both monotonic in the same direction.")

    return pem

