################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as cpp_vector
from pyne cimport cpp_endf2
from pyne cimport cpp_endf_mt

cdef extern from "endf.h" namespace "pyne::endf":

    cdef cppclass endf_id:
        # constructors
        endf_id() except +

        # attributes
        int mat
        int mf
        int mt

        # methods

        pass



cdef extern from "endf.h" namespace "endf":

    cdef cppclass library:
        # constructors
        library() except +

        # attributes
        cpp_map[endf_id, cpp_endf_mt.mt_base *] contents

        # methods
        cpp_vector[cpp_vector[int]] get_content_list() except +
        void read_endf() except +
        void read_endf(std_string) except +
        pass




{'cpppxd_footer': '', 'pyx_header': '', 'pxd_header': '', 'pxd_footer': '', 'cpppxd_header': '', 'pyx_footer': ''}