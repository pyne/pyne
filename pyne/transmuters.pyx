# Cython Imports
from libcpp.map cimport map as cpp_map
from libcpp.vector cimport vector as cpp_vector

# Local imports
cimport cpp_transmuters
cimport pyne.stlcontainers as conv
import pyne.stlcontainers as conv
from pyne cimport nucname
from pyne import nucname

# startup numpy
cimport numpy as np
import numpy as np
np.import_array()
np.import_ufunc()


def cram(A, n0, int order=14):
    """Basic CRAM solver that takes a (flat) A matrix and an inital nuclide
    atom fraction composition map and returns the value.

    Parameters
    ----------
    A : 1D array-like
        The transmutation matrix [unitless]
    n0 : Mapping
        The initial compositions [atom fraction]
    order : int, optional
        The order of approximation, default 14.

    Returns
    -------
    n1 : Mapping
        The result of the transmutation [atom fraction]
    """
    # convert A
    A = np.asarray(A, dtype=np.float64)
    cdef int Alen = len(A)
    cdef double* Aptr = <double*> np.PyArray_DATA(A)
    cdef cpp_vector[double] cpp_A = cpp_vector[double]()
    cpp_A.reserve(Alen)
    cpp_A.assign(Aptr, Aptr + Alen)
    # Convert n0
    cdef cpp_map[int, double] cpp_n0 = cpp_map[int, double]()
    for key, value in n0.items():
        cpp_n0[nucname.id(key)] = value
    # call cram and return
    cdef cpp_map[int, double] cpp_n1 = cpp_transmuters.cram(cpp_A, cpp_n0, order)
    cdef conv._MapIntDouble n1 = conv.MapIntDouble()
    n1.map_ptr = new cpp_map[int, double](cpp_n1)
    return n1
