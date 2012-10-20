"""Cython header for C/C++ standard library functionailty."""

cdef extern from 'stdlib.h':
    double atof(char*)


cdef extern from 'string.h':
    ctypedef char const_char "const char"
    char* strtok(char*, char*)
    char* strcpy(char*, char*)
    void* memcpy(void*, void*, size_t)

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char * c_str()

