################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
"""
"""
cimport dtypes
cimport numpy as np
from libc.stdlib cimport free
from libcpp cimport bool as cpp_bool
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as cpp_vector

import numpy as np

np.import_array()



cdef class AliasTable:
    """Constructor
    
    Attributes
    ----------
    n (int) :
    prob (std::vector< double >) : Number of bins in the PDF.
    alias (std::vector< int >) : Probabilities.
    
    
    Methods
    -------
    AliasTable
    ~AliasTable
    sample_pdf
    
    Notes
    -----
    This class was defined in source_sampling.h
    
    The class is found in the "pyne" namespace"""



    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults
        self._alias = None
        self._prob = None

    def __init__(self, p):
        """AliasTable(self, p)
        Constructor
        
        Parameters
        ----------
        p : std::vector< double >
        
        Returns
        -------
        None
        
        """
        cdef cpp_vector[double] p_proxy
        cdef int ip
        cdef int p_size
        cdef double * p_data
        # p is a ('vector', 'float64', 0)
        p_size = len(p)
        if isinstance(p, np.ndarray) and (<np.ndarray> p).descr.type_num == np.NPY_FLOAT64:
            p_data = <double *> np.PyArray_DATA(<np.ndarray> p)
            p_proxy = cpp_vector[double](<size_t> p_size)
            for ip in range(p_size):
                p_proxy[ip] = p_data[ip]
        else:
            p_proxy = cpp_vector[double](<size_t> p_size)
            for ip in range(p_size):
                p_proxy[ip] = <double> p[ip]
        self._inst = new cpp_source_sampling.AliasTable(p_proxy)
    
    
    def __dealloc__(self):
        if self._free_inst and self._inst is not NULL:
            free(self._inst)

    # attributes
    property alias:
        """no docstring for alias, please file a bug report!"""
        def __get__(self):
            cdef np.ndarray alias_proxy
            cdef np.npy_intp alias_proxy_shape[1]
            if self._alias is None:
                alias_proxy_shape[0] = <np.npy_intp> (<cpp_source_sampling.AliasTable *> self._inst).alias.size()
                alias_proxy = np.PyArray_SimpleNewFromData(1, alias_proxy_shape, np.NPY_INT32, &(<cpp_source_sampling.AliasTable *> self._inst).alias[0])
                self._alias = alias_proxy
            return self._alias
    
        def __set__(self, value):
            cdef cpp_vector[int] value_proxy
            cdef int ivalue
            cdef int value_size
            cdef int * value_data
            # value is a ('vector', 'int32', 0)
            value_size = len(value)
            if isinstance(value, np.ndarray) and (<np.ndarray> value).descr.type_num == np.NPY_INT32:
                value_data = <int *> np.PyArray_DATA(<np.ndarray> value)
                value_proxy = cpp_vector[int](<size_t> value_size)
                for ivalue in range(value_size):
                    value_proxy[ivalue] = value_data[ivalue]
            else:
                value_proxy = cpp_vector[int](<size_t> value_size)
                for ivalue in range(value_size):
                    value_proxy[ivalue] = <int> value[ivalue]
            (<cpp_source_sampling.AliasTable *> self._inst).alias = value_proxy
            self._alias = None
    
    
    property n:
        """no docstring for n, please file a bug report!"""
        def __get__(self):
            return int((<cpp_source_sampling.AliasTable *> self._inst).n)
    
        def __set__(self, value):
            (<cpp_source_sampling.AliasTable *> self._inst).n = <int> value
    
    
    property prob:
        """no docstring for prob, please file a bug report!"""
        def __get__(self):
            cdef np.ndarray prob_proxy
            cdef np.npy_intp prob_proxy_shape[1]
            if self._prob is None:
                prob_proxy_shape[0] = <np.npy_intp> (<cpp_source_sampling.AliasTable *> self._inst).prob.size()
                prob_proxy = np.PyArray_SimpleNewFromData(1, prob_proxy_shape, np.NPY_FLOAT64, &(<cpp_source_sampling.AliasTable *> self._inst).prob[0])
                self._prob = prob_proxy
            return self._prob
    
        def __set__(self, value):
            cdef cpp_vector[double] value_proxy
            cdef int ivalue
            cdef int value_size
            cdef double * value_data
            # value is a ('vector', 'float64', 0)
            value_size = len(value)
            if isinstance(value, np.ndarray) and (<np.ndarray> value).descr.type_num == np.NPY_FLOAT64:
                value_data = <double *> np.PyArray_DATA(<np.ndarray> value)
                value_proxy = cpp_vector[double](<size_t> value_size)
                for ivalue in range(value_size):
                    value_proxy[ivalue] = value_data[ivalue]
            else:
                value_proxy = cpp_vector[double](<size_t> value_size)
                for ivalue in range(value_size):
                    value_proxy[ivalue] = <double> value[ivalue]
            (<cpp_source_sampling.AliasTable *> self._inst).prob = value_proxy
            self._prob = None
    
    
    # methods
    def sample_pdf(self, rand1, rand2):
        """sample_pdf(self, rand1, rand2)
        Samples the alias table
        
        Parameters
        ----------
        rand2 : double
        
        rand1 : double
        
        Returns
        -------
        res1 : int
        
        """
        cdef int rtnval
        rtnval = (<cpp_source_sampling.AliasTable *> self._inst).sample_pdf(<double> rand1, <double> rand2)
        return int(rtnval)
    
    
    

    pass





cdef class Sampler:
    """Constuctor for analog and uniform sampling
    
    Attributes
    ----------
    filename (std::string) : MOAB mesh file path.
    src_tag_name (std::string) : Unbiased source density
        distribution.
    bias_tag_name (std::string) : Biased source density
        distribution.
    e_bounds (std::vector< double >) : Energy boundaries.
    num_e_groups (int) : Number of groups in tag
    num_bias_groups (int) : Number of groups tag
    mode (Mode) : Problem mode: analog, uniform, user.
    mesh (MBInterface *) : MOAB mesh.
    num_ves (int) : Number of mesh volume elements on
    ve_type (MBEntityType) : Type of mesh volume: MBTET or MBHEX.
    verts_per_ve (int) : Number of verticles per mesh volume
        element.
    all_edge_points (std::vector< ) : Four connected points on a VE.
    biased_weights (std::vector< double >) : Birth weights for
        biased sampling.
    at (None) : Alias table used for sampling.
    
    
    Methods
    -------
    Sampler
    ~Sampler
    mesh_geom_data
    mesh_tag_data
    normalize_pdf
    num_groups
    particle_birth
    read_bias_pdf
    sample_e
    sample_w
    sample_xyz
    setup
    
    Notes
    -----
    This class was defined in source_sampling.h
    
    The class is found in the "pyne" namespace"""



    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults


    def _sampler_sampler_0(self, filename, src_tag_name, e_bounds, bias_tag_name):
        """Sampler(self, filename, src_tag_name, e_bounds, bias_tag_name)
         This method was overloaded in the C-based source. To overcome
        this we ill put the relevant docstring for each version below.
        Each version will begin with a line of # characters.
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        bias_tag_name : std::string
        
        filename : std::string
        
        Returns
        -------
        None
        
        ################################################################
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        uniform : bool
        
        filename : std::string
        
        Returns
        -------
        None
        
        """
        cdef char * filename_proxy
        cdef char * src_tag_name_proxy
        cdef cpp_vector[double] e_bounds_proxy
        cdef int ie_bounds
        cdef int e_bounds_size
        cdef double * e_bounds_data
        cdef char * bias_tag_name_proxy
        filename_bytes = filename.encode()
        src_tag_name_bytes = src_tag_name.encode()
        # e_bounds is a ('vector', 'float64', 0)
        e_bounds_size = len(e_bounds)
        if isinstance(e_bounds, np.ndarray) and (<np.ndarray> e_bounds).descr.type_num == np.NPY_FLOAT64:
            e_bounds_data = <double *> np.PyArray_DATA(<np.ndarray> e_bounds)
            e_bounds_proxy = cpp_vector[double](<size_t> e_bounds_size)
            for ie_bounds in range(e_bounds_size):
                e_bounds_proxy[ie_bounds] = e_bounds_data[ie_bounds]
        else:
            e_bounds_proxy = cpp_vector[double](<size_t> e_bounds_size)
            for ie_bounds in range(e_bounds_size):
                e_bounds_proxy[ie_bounds] = <double> e_bounds[ie_bounds]
        bias_tag_name_bytes = bias_tag_name.encode()
        self._inst = new cpp_source_sampling.Sampler(std_string(<char *> filename_bytes), std_string(<char *> src_tag_name_bytes), e_bounds_proxy, std_string(<char *> bias_tag_name_bytes))
    
    
    def _sampler_sampler_1(self, filename, src_tag_name, e_bounds, uniform, threshold=0):
        """Sampler(self, filename, src_tag_name, e_bounds, uniform)
         This method was overloaded in the C-based source. To overcome
        this we ill put the relevant docstring for each version below.
        Each version will begin with a line of # characters.
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        bias_tag_name : std::string
        
        filename : std::string
        
        Returns
        -------
        None
        
        ################################################################
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        uniform : bool

        threshold : double
        
        filename : std::string
        
        Returns
        -------
        None
        
        """
        cdef char * filename_proxy
        cdef char * src_tag_name_proxy
        cdef cpp_vector[double] e_bounds_proxy
        cdef int ie_bounds
        cdef int e_bounds_size
        cdef double * e_bounds_data
        filename_bytes = filename.encode()
        src_tag_name_bytes = src_tag_name.encode()
        # e_bounds is a ('vector', 'float64', 0)
        e_bounds_size = len(e_bounds)
        if isinstance(e_bounds, np.ndarray) and (<np.ndarray> e_bounds).descr.type_num == np.NPY_FLOAT64:
            e_bounds_data = <double *> np.PyArray_DATA(<np.ndarray> e_bounds)
            e_bounds_proxy = cpp_vector[double](<size_t> e_bounds_size)
            for ie_bounds in range(e_bounds_size):
                e_bounds_proxy[ie_bounds] = e_bounds_data[ie_bounds]
        else:
            e_bounds_proxy = cpp_vector[double](<size_t> e_bounds_size)
            for ie_bounds in range(e_bounds_size):
                e_bounds_proxy[ie_bounds] = <double> e_bounds[ie_bounds]
        self._inst = new cpp_source_sampling.Sampler(std_string(<char *> filename_bytes), std_string(<char *> src_tag_name_bytes), e_bounds_proxy, <bint> uniform, threshold)
    
    
    _sampler_sampler_0_argtypes = frozenset(((0, str), (1, str), (2, np.ndarray), (3, str), ("filename", str), ("src_tag_name", str), ("e_bounds", np.ndarray), ("bias_tag_name", str)))
    _sampler_sampler_1_argtypes = frozenset(((0, str), (1, str), (2, np.ndarray), (3, bool), ("filename", str), ("src_tag_name", str), ("e_bounds", np.ndarray), ("uniform", bool)))
    
    def __init__(self, *args, **kwargs):
        """Sampler(self, filename, src_tag_name, e_bounds, uniform)
         This method was overloaded in the C-based source. To overcome
        this we ill put the relevant docstring for each version below.
        Each version will begin with a line of # characters.
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        bias_tag_name : std::string
        
        filename : std::string
        
        Returns
        -------
        None
        
        ################################################################
        
        Constuctor for analog and uniform sampling
        
        Parameters
        ----------
        e_bounds : std::vector< double >
        
        src_tag_name : std::string
        
        uniform : bool
        
        filename : std::string
        
        Returns
        -------
        None
        
        """
        types = set([(i, type(a)) for i, a in enumerate(args)])
        types.update([(k, type(v)) for k, v in kwargs.items()])
        # vtable-like dispatch for exactly matching types
        if types <= self._sampler_sampler_0_argtypes:
            self._sampler_sampler_0(*args, **kwargs)
            return
        if types <= self._sampler_sampler_1_argtypes:
            self._sampler_sampler_1(*args, **kwargs)
            return
        # duck-typed dispatch based on whatever works!
        try:
            self._sampler_sampler_0(*args, **kwargs)
            return
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            self._sampler_sampler_1(*args, **kwargs)
            return
        except (RuntimeError, TypeError, NameError):
            pass
        raise RuntimeError('method __init__() could not be dispatched')
    
    def __dealloc__(self):
        if self._free_inst and self._inst is not NULL:
            free(self._inst)

    # attributes

    # methods
    def particle_birth(self, rands):
        """particle_birth(self, rands)
        Samples particle birth parameters
        
        Parameters
        ----------
        rands : std::vector< double >
        
        Returns
        -------
        res1 : std::vector< double >
        
        """
        cdef cpp_vector[double] rands_proxy
        cdef int irands
        cdef int rands_size
        cdef double * rands_data
        cdef cpp_vector[double] rtnval
        
        cdef np.npy_intp rtnval_proxy_shape[1]
        # rands is a ('vector', 'float64', 0)
        rands_size = len(rands)
        if isinstance(rands, np.ndarray) and (<np.ndarray> rands).descr.type_num == np.NPY_FLOAT64:
            rands_data = <double *> np.PyArray_DATA(<np.ndarray> rands)
            rands_proxy = cpp_vector[double](<size_t> rands_size)
            for irands in range(rands_size):
                rands_proxy[irands] = rands_data[irands]
        else:
            rands_proxy = cpp_vector[double](<size_t> rands_size)
            for irands in range(rands_size):
                rands_proxy[irands] = <double> rands[irands]
        rtnval = (<cpp_source_sampling.Sampler *> self._inst).particle_birth(rands_proxy)
        rtnval_proxy_shape[0] = <np.npy_intp> rtnval.size()
        rtnval_proxy = np.PyArray_SimpleNewFromData(1, rtnval_proxy_shape, np.NPY_FLOAT64, &rtnval[0])
        rtnval_proxy = np.PyArray_Copy(rtnval_proxy)
        return rtnval_proxy
    
    
    

    pass






{'cpppxd_footer': '', 'pyx_header': '', 'pxd_header': '', 'pxd_footer': '', 'cpppxd_header': '', 'pyx_footer': ''}
