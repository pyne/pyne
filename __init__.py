import ctypes, ctypes.util
import sys
import os.path

# Attempt to find the dagmc_bridge shared library
# ctypes.util.find_library is (sometimes) useful, but does not work the
# same way on all platforms.
_dirname = os.path.dirname( __file__ )
if sys.platform == 'darwin':
    _bridgepath = ctypes.util.find_library(os.path.join(_dirname,'dagmc_bridge','dagmc_bridge'))
elif sys.platform.startswith('linux'):
    _bridgepath = os.path.join(_dirname,'dagmc_bridge','dagmc_bridge.so')
else:
    print(sys.stderr, 'Warning: unknown platform, lib import may fail')
    _bridgepath = ctypes.util.find_library('dagmc_bridge')

daglib = ctypes.CDLL(_bridgepath)

# Set the entity handle type; I don't know a cleaner way than to ask the daglib
# to return the byte width of the type
def get_EH_type():
    EH_size = daglib.dag_ent_handle_size()
    if( EH_size == 8):
        EH_t = ctypes.c_uint64
    elif( EH_size == 4):
        EH_t = ctypes.c_uint32
    else:
        raise TypeError("Unrecognized entity handle size in dagmc library: " + int(EH_size) )
    return type("EntityHandle",(EH_t,),{})

EntityHandle = get_EH_type()
_ErrorCode   = type( "ErrorCode", (ctypes.c_int,),{})

class DagmcError( Exception ):
    
    pass

