"""C++ wrapper for extra types header."""
from libc.stdio cimport FILE

# Dirty ifdef, else, else preprocessor hack
# see http://comments.gmane.org/gmane.comp.python.cython.user/4080
cdef extern from *:
    cdef void emit_ifc "#if defined(__STDC__) //" ()
    cdef void emit_ifcpp "#if defined(__cplusplus) //" ()
    cdef void emit_elifc "#elif defined(__STDC__) //" ()
    cdef void emit_elifcpp "#elif defined(__cplusplus) //" ()
    cdef void emit_else "#else //" ()
    cdef void emit_endif "#endif //" ()

ctypedef unsigned char uchar
ctypedef long long int64
ctypedef unsigned short uint16
ctypedef unsigned int uint32
ctypedef unsigned long long uint64
ctypedef long double float128

cdef extern from "extra_types.h":

    ctypedef struct complex_t "xd_complex_t":
        double re
        double im

cdef complex_t py2c_complex(object pyv)

cdef extern from "Python.h":

    object PyFile_FromFile(FILE *fp, char *name, char *mode, int (*close)(FILE*))
    FILE* PyFile_AsFile(object p)


#emit_ifcpp()
#cdef extern from "<exception>" namespace "std":

#    cdef cppclass exception:
#        exception()
#        exception(const exception&)
#        exception& operator= (const exception&)
#        ~exception()
#        const char * what()

#emit_endif()
