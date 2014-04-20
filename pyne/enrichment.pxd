"""Python header for enrichment library."""

cimport cpp_enrichment
cimport pyne.material

from warnings import warn

warn(__name__ + " is not yet V&V compliant.", ImportWarning)

cdef class Cascade:
    cdef cpp_enrichment.Cascade * _inst
    cdef public pyne.material._Material _mat_feed
    cdef public pyne.material._Material _mat_prod
    cdef public pyne.material._Material _mat_tail
