"""Python wrapper for C++ standard library functionailty."""

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()

