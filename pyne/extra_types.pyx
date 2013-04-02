cdef complex_t py2c_complex(object pyv):
    cdef complex_t cv = complex_t()
    pyv = complex(pyv)
    cv.re = pyv.real
    cv.im = pyv.imag
    return cv

