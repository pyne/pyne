"""Cython-based utils to be imported into utils."""

from __future__ import division, unicode_literals
from libc.stdlib cimport malloc, free
from libc.stdlib cimport atof
from libc.string cimport strtok, strcpy, strncpy

cimport numpy as np
cimport pyne.cpp_utils
from cython.operator cimport dereference as deref
import numpy as np

def fromstring_split(s, sep=None, dtype=float):
    """A replacement for numpy.fromstring() using the Python str.split()
    and np.array().

    Parameters
    ----------
    s : str
        String of data.
    sep : str or None
        String of separator characters, has the same meaning as in
        str.split().
    dtype : np.dtype
        Numpy dtype to cast elements enough.

    Returns
    -------
    data : ndarray, 1d
        Will always return a 1d array of dtype.  You must reshape to the
        appropriate shape.

    See Also
    --------
    fromstring_token : May faster depending on the data.

    """
    cdef list rawdata
    rawdata = s.split(sep)
    return np.array(rawdata, dtype=dtype)


def fromstring_token(s, sep=" ", bint inplace=False, int maxsize=-1):
    """A replacement for numpy.fromstring() using the C standard
    library atof() and strtok() functions.

    Parameters
    ----------
    s : str
        String of data.
    sep : str
        String of separator characters.  Unlike numpy.fromstring(),
        all characters are separated on independently.
    inplace : bool
        Whether s should tokenized in-place or whether a copy should
        be made.  If done in-place, the first instance of sep between
        any tokens will replaced with the NULL character.
    maxsize : int
        Specifies the size of the array to pre-allocate.  If negative,
        this will be set to the maximum possible number of elements,
        ie len(s)/2 + 1.

    Returns
    -------
    data : ndarray, 1d, float64
        Will always return a 1d float64 array.  You must cast and reshape
        to the appropriate type and shape.

    See Also
    --------
    fromstring_split : May faster depending on the data.

    """
    cdef char* cstring
    cdef char* cs
    cdef char* csep
    cdef int i, I
    cdef np.ndarray[np.float64_t, ndim=1] cdata

    s_bytes = s.encode()
    I = len(s_bytes)
    sep_bytes = sep.encode()
    csep = sep_bytes

    if inplace:
        cs = s_bytes
    else:
        cs = <char *> malloc(I * sizeof(char))
        strcpy(cs, s_bytes)

    if maxsize < 0:
        maxsize = (I // 2) + 1

    data = np.empty(maxsize, dtype=np.float64)
    cdata = data

    i = 0
    cstring = strtok(cs, csep)
    while cstring != NULL:
        cdata[i] = atof(cstring)
        cstring = strtok(NULL, csep)
        i += 1

    if not inplace:
        free(cs)

    data = data[:i].copy()
    return data


def endftod(s):
    """Converts a string from ENDF number format to float64.

    Parameters
    ----------
    s : char *
        Plain string to convert.

    Returns
    -------
    float64
    """
    cdef char * cs
    if isinstance(s, str):
        s = s.encode()
    cs = s
    return pyne.cpp_utils.endftod(cs)


def use_fast_endftod():
    """ Switches to fast ENDF string parser"""
    pyne.cpp_utils.use_fast_endftod()


def fromendf_tok(s):
    """A replacement for numpy.fromstring().

    Parameters:
    -----------
    s : str
        String of data, consisting of complete lines of ENDF data.

    Returns:
    --------
    data : ndarray, 1d, float64
        Will always return a 1d float64 array.  You must reshape to the
        appropriate shape.
    """
    cdef char * cs
    if isinstance(s, str):
        s = s.encode()
    cs = s
    cdef int i, num_entries
    cdef char entry[12]
    cdef long pos = 0
    cdef np.ndarray[np.float64_t, ndim=1] cdata
    i = 0
    num_entries = len(cs)//81 * 6
    cdata = np.empty(num_entries, dtype=np.float64)
    while i < num_entries:
        pos = i*11 + i//6 * 15
        strncpy(entry, cs+pos, 11)
        cdata[i] = pyne.cpp_utils.endftod(entry)
        i += 1
    return cdata
