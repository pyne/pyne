################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
# With some edits from BaM241                  #
#                    Come on, guys. I mean it! #
################################################


from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector
from pyne cimport cpp_tally
from pyne cimport stlcontainers

cdef vector[double] to_vector_double(value)
cdef vector[int] to_vector_int(value)
cdef vector[std_string] to_vector_string(value)


cdef class Tally:
    cdef void * _inst
    cdef public bint _free_inst
    cdef public stlcontainers._MapStrStr _rx2fluka
    cdef public stlcontainers._MapStrStr _rx2mcnp5
    cdef public stlcontainers._MapStrStr _rx2mcnp6
    pass


