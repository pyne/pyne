import ctypes
from contextlib import contextmanager
from numpy import ndarray, float64
from numpy.ctypeslib import as_ctypes, ndpointer

from pydagmc import daglib, EntityHandle, _ErrorCode, DagmcError


lib = daglib

lib.dag_version.restype = ctypes.c_float
lib.dag_version.argtypes = []

lib.dag_rev_version.restype = ctypes.c_uint
lib.dag_rev_version.argtypes = []

def _geom_dim_check( result, func, arguments ):
    if arguments[0] not in (2,3): raise DagmcError('Incorrect geometric dimension: '+str(arguments[0]))
    return result

lib.geom_id_list.restype = ctypes.POINTER( ctypes.c_int )
lib.geom_id_list.argtypes = [ ctypes.c_int, ctypes.POINTER( ctypes.c_int ) ]
lib.geom_id_list.errcheck = _geom_dim_check

lib.handle_from_id.restype = EntityHandle
lib.handle_from_id.argtypes = [ ctypes.c_int, ctypes.c_int ]
lib.handle_from_id.errcheck = _geom_dim_check

lib.id_from_handle.restype = ctypes.c_int
lib.id_from_handle.argtypes = [ EntityHandle ]

def _moab_error_check( result, func, arguments ):
    if result.value != 0:
        raise DagmcError( "Error code " + str(result.value) + " returned from " + func.__name__ )
    return result

def _returns_moab_errors( function ):
    function.restype = _ErrorCode
    function.errcheck = _moab_error_check

_returns_moab_errors( lib.dag_load )
lib.dag_load.argtypes = [ ctypes.c_char_p ]

#lib.dag_ray_fire



class listPOINTER:
    def __init__(self, etype, length):
        self.etype = etype
        self.length = length

    def from_param(self, param):
        if isinstance( param, ndarray ) and len(param) == self.length:
            return as_ctypes( param )
        elif isinstance( param, (list,tuple) ) and len(param) == self.length:
            return (self.etype * len(param))(*param)
        else:
            raise DagmcError( "Problem with listPOINTER" )


#_vec3 = listPOINTER( ctypes.c_double, 3 )
_vec3 = ndpointer( dtype=float64, shape=(3,), flags='CONTIGUOUS' )

_returns_moab_errors( lib.dag_pt_in_vol )
lib.dag_pt_in_vol.argtypes = [ EntityHandle, _vec3, ctypes.POINTER(ctypes.c_int),
                               _vec3, ctypes.c_void_p ]

lib.dag_alloc_ray_history.restype = ctypes.c_void_p
lib.dag_alloc_ray_history.argtypes = []

lib.dag_dealloc_ray_history.restype = None
lib.dag_dealloc_ray_history.argtypes = [ ctypes.c_void_p ]

@contextmanager
def _ray_history():
    history = lib.dag_alloc_ray_history()
    yield history
    lib.dag_dealloc_ray_history( history )

_returns_moab_errors( lib.dag_ray_fire )
lib.dag_ray_fire.argtypes = [ EntityHandle, _vec3, _vec3,
                              ctypes.POINTER( EntityHandle ), ctypes.POINTER( ctypes.c_double ), 
                              ctypes.c_void_p, ctypes.c_double ]


_returns_moab_errors( lib.dag_next_vol )
lib.dag_next_vol.argtypes = [ EntityHandle, EntityHandle, ctypes.POINTER( EntityHandle ) ]

lib.vol_is_graveyard.restype = ctypes.c_int
lib.vol_is_graveyard.argtypes = [EntityHandle]

lib.vol_is_implicit_complement.restype = ctypes.c_int
lib.vol_is_implicit_complement.argtypes = [EntityHandle]

_returns_moab_errors( lib.get_volume_metadata )
lib.get_volume_metadata.argtypes = [ EntityHandle, ctypes.POINTER( ctypes.c_int ),
                                     ctypes.POINTER( ctypes.c_double ), ctypes.POINTER( ctypes.c_double ) ]

_returns_moab_errors( lib.get_volume_boundary )
lib.get_volume_boundary.argtypes = [ EntityHandle, _vec3, _vec3 ]
