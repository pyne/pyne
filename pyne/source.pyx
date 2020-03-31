from pyne cimport stlcontainers
from libc.stdlib cimport free
from libcpp.string cimport string as std_string
import numpy as np
cimport numpy as np

import collections
from pyne import stlcontainers

from warnings import warn
from pyne.utils import QAWarning

warn(__name__ + " is not yet QA compliant.", QAWarning)

cdef class PointSource:
    """

    Attributes
    ----------


    Methods
    -------

    Notes
    -----
    This class was defined in source.h

    The class is found in the "pyne" namespace"""

    # constuctors
    def __cinit__(self, x=0.0, y=0.0, z=0.0, i=0.0, j=0.0, k=0.0, E=14., p="n", weight=1):
        p_bytes = p.encode()
        self._inst = new cpp_source.PointSource(x, y, z, i, j, k, E, std_string(<char *> p_bytes), weight)
        self._free_inst = True

    def mcnp(self):
        cdef std_string card
        card = self._inst.mcnp()
        return card.decode()
    
    
    def __dealloc__(self):
        """PointSource C++ destructor."""
        if self._free_inst:
            del self._inst
