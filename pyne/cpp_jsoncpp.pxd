"""C++ wrapper for jsoncpp."""

cdef extern from "json/json.h" namespace "Json":

    cdef cppclass Value:
        # Constuctors
        Value() except +
