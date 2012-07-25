"""Python header for enrichment library."""
cimport cpp_enrichment

cdef class Cascade:
    cdef cpp_enrichment.Cascade * ptr
