################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


cimport dtypes
cimport endf_mt
cimport stlcontainers
from libcpp.map cimport map as cpp_map
from pyne cimport cpp_endf2
from pyne cimport cpp_endf_mt



cdef class endf_id:
    cdef void * _inst
    cdef public bint _free_inst
    pass





cdef class library:
    cdef void * _inst
    cdef public bint _free_inst
    cdef public stlcontainers._Mapendf_idmt_base _contents
    pass




{'cpppxd_footer': '', 'pyx_header': '', 'pxd_header': '', 'pxd_footer': '', 'cpppxd_header': '', 'pyx_footer': ''}