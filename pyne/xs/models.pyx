"""This module provides physical cross-section models and helper functions."""

cimport numpy as np
import numpy as np

from pyne cimport nucname
from pyne import nucname


##############################
### Partial group collapse ###
##############################

def partial_energy_matrix_inc(np.ndarray[np.float64_t, ndim=1] E_g, np.ndarray[np.float64_t, ndim=1] E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    Here, both of the energy arrays must be monotonically increasing.

    Parameters
    ----------
    E_g : 1d numpy float array 
        Lower resolution energy group structure [MeV] that is of length G+1. 
        Ordered from lowest-to-highest energy.
    E_n : 1d numpy float array 
        Higher resolution energy group structure [MeV] that is of length N+1. 
        Ordered from lowest-to-highest energy.

    Returns
    -------
    pem : 2d numpy float array of fractions 
        This is a GxN sized matrix that when dotted with a high-resolution 
        flux (or cross section) produces a low-resolution flux (or cross section).
    """
    # Some convienence parameters
    cdef Py_ssize_t G, N
    G = E_g.shape[0] - 1
    N = E_n.shape[0] - 1

    index_E_n = np.arange(N+1)

    # Get the interior points for each gth group in n-space
    inner_mask = np.array([(E_g[g] <= E_n) & (E_n <= E_g[g+1]) for g in range(G)])

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
            pem[g][lig-1] = (E_n[lig] - E_g[g]) / (E_n[lig] - E_n[lig-1])

        # Upper bound
        uig = upper_index[g]
        if uig < N:
            pem[g][uig] = (E_g[g+1] - E_n[uig]) / (E_n[uig+1] - E_n[uig])

    return pem


def partial_energy_matrix(E_g, E_n):
    """Gerenates a matrix of fractional values that may be used to converts a high-resolution 
    flux array with group structure E_n to a low-resolution flux array with group-structure E_g.
    The group structures must have the same monotonicity.

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
    cdef np.ndarray[np.float64_t, ndim=1] pem

    E_g = np.asarray(E_g, dtype=float)
    E_n = np.asarray(E_n, dtype=float)

    if (E_g[:-1] < E_g[1:]).all() and (E_n[:-1] < E_n[1:]).all():
        # Both energy arrays are monotonically increasing
        assert E_n[0] <= E_g[0]
        assert E_g[-1] <= E_n[-1]
        pem = partial_energy_matrix_inc(E_g, E_n)
    else:
        raise ValueError("E_g and E_n are not both monotonic in the same direction.")

    return pem

