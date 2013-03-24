cimport extra_types

cdef extra_types.complex_t py2c_complex(object pyv):
    cdef extra_types.complex_t cv = extra_types.complex_t()
    pyv = complex(pyv)
    cv.re = pyv.real
    cv.im = pyv.imag
    return cv

